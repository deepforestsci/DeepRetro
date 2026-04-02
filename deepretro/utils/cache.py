"""
Explicit in-memory cache primitives for expensive library operations.

Provides a small process-local ``CacheManager`` with tag-based invalidation,
TTL support, and lightweight statistics. Callers must instantiate and pass
cache objects explicitly; this module does not create shared caches implicitly.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import sys
import time
from dataclasses import dataclass
from typing import Any

import structlog

logger = structlog.get_logger(__name__)

_MISS = object()


@dataclass
class CacheStats:
    """Statistics for cache hits, misses, approximate size, and entry count."""

    hits: int
    misses: int
    size_bytes: int
    num_entries: int


@dataclass
class _CacheEntry:
    """Single in-memory cache entry."""

    value: Any
    expires_at: float | None
    tag: str | None


def _make_args_hash(*args: Any, **kwargs: Any) -> str:
    """
    Generate a deterministic hash of arguments for cache keying.

    Tries JSON first for common types; falls back to pickle for complex objects.
    """
    payload = {"args": args, "kwargs": kwargs}
    try:
        raw = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()
    except (TypeError, ValueError):
        raw = pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)
        return hashlib.sha256(raw).hexdigest()


def make_cache_key(namespace: str, *args: Any, version: int = 1, **kwargs: Any) -> str:
    """
    Build a deterministic cache key for a namespaced operation.

    Parameters
    ----------
    namespace : str
        Stable operation name, such as ``"run_az"``.
    *args : Any
        Positional arguments that affect the cached result.
    version : int, optional
        Cache version. Bump when behavior changes and old entries should be
        invalidated, by default 1.
    **kwargs : Any
        Keyword arguments that affect the cached result.

    Returns
    -------
    str
        A deterministic key suitable for ``CacheManager.get`` and ``set``.
    """
    if not namespace:
        raise ValueError("namespace must be a non-empty string")
    args_hash = _make_args_hash(*args, **kwargs)
    return f"v{version}:{namespace}:{args_hash}"


class CacheManager:
    """
    Process-local in-memory cache manager with tag support and TTL.

    Each ``CacheManager`` instance owns its own dictionary, so cache state is
    isolated unless the caller explicitly reuses the same object.
    """

    def __init__(self) -> None:
        """Initialize an empty in-memory cache."""
        self._entries: dict[str, _CacheEntry] = {}
        self._tags: dict[str, set[str]] = {}
        self._hits = 0
        self._misses = 0
        self._log = logger.bind(component="cache")

    def _purge_if_expired(self, key: str) -> bool:
        """Remove the key if it has expired."""
        entry = self._entries.get(key)
        if entry is None:
            return False
        if entry.expires_at is None or entry.expires_at > time.monotonic():
            return False
        self._delete_key(key)
        self._log.debug("cache.expired", key=key)
        return True

    def _delete_key(self, key: str) -> bool:
        """Remove a key from the cache and tag index if present."""
        entry = self._entries.pop(key, None)
        if entry is None:
            return False
        if entry.tag is not None:
            tagged_keys = self._tags.get(entry.tag)
            if tagged_keys is not None:
                tagged_keys.discard(key)
                if not tagged_keys:
                    self._tags.pop(entry.tag, None)
        return True

    def _purge_expired_entries(self) -> None:
        """Drop expired keys so stats reflect live entries only."""
        for key in list(self._entries):
            self._purge_if_expired(key)

    def _estimate_size_bytes(self) -> int:
        """Return a shallow approximation of in-memory cache size."""
        size = sys.getsizeof(self._entries) + sys.getsizeof(self._tags)
        for key, entry in self._entries.items():
            size += sys.getsizeof(key)
            size += sys.getsizeof(entry)
            size += sys.getsizeof(entry.value)
        for tag, keys in self._tags.items():
            size += sys.getsizeof(tag)
            size += sys.getsizeof(keys)
        return size

    def get(self, key: str, default: Any = _MISS) -> Any:
        """
        Retrieve a value by key.

        Parameters
        ----------
        key : str
            Cache key.
        default : Any, optional
            Value returned when the key is not cached, by default an internal
            sentinel object.

        Returns
        -------
        Any
            Cached value, or ``default`` if not found.
        """
        if self._purge_if_expired(key):
            self._misses += 1
            self._log.debug("cache.miss", key=key)
            return default

        entry = self._entries.get(key)
        if entry is None:
            self._misses += 1
            self._log.debug("cache.miss", key=key)
            return default

        self._hits += 1
        self._log.debug("cache.hit", key=key)
        return entry.value

    def set(
        self,
        key: str,
        value: Any,
        *,
        expire: float | None = None,
        tag: str | None = None,
    ) -> None:
        """
        Store a value with optional TTL and tag.

        Parameters
        ----------
        key : str
            Cache key.
        value : Any
            Value to store.
        expire : float | None, optional
            Time-to-live in seconds. None means no expiry.
        tag : str | None, optional
            Tag for later eviction via ``evict_tag``.
        """
        self._delete_key(key)
        expires_at = None if expire is None else time.monotonic() + expire
        self._entries[key] = _CacheEntry(value=value, expires_at=expires_at, tag=tag)
        if tag is not None:
            self._tags.setdefault(tag, set()).add(key)
        self._log.debug("cache.set", key=key, expire=expire, tag=tag)

    def evict_tag(self, tag: str) -> int:
        """
        Remove all live entries with the given tag.

        Parameters
        ----------
        tag : str
            Tag identifying entries to remove.

        Returns
        -------
        int
            Number of entries evicted.
        """
        keys = list(self._tags.get(tag, set()))
        removed = 0
        for key in keys:
            removed += int(self._delete_key(key))
        self._log.info("cache.evict_tag", tag=tag, removed=removed)
        return removed

    def clear(self) -> None:
        """Remove all entries from this cache instance."""
        self._entries.clear()
        self._tags.clear()
        self._log.info("cache.clear")

    def stats(self) -> CacheStats:
        """
        Return cache statistics.

        Returns
        -------
        CacheStats
            Hits, misses, approximate size in bytes, and live entry count.
        """
        self._purge_expired_entries()
        return CacheStats(
            hits=self._hits,
            misses=self._misses,
            size_bytes=self._estimate_size_bytes(),
            num_entries=len(self._entries),
        )

    def close(self) -> None:
        """Release cache contents for compatibility with previous callers."""
        self.clear()
