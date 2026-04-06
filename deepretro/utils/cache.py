"""
Explicit in-memory cache primitives for expensive library operations.

This module provides a small process-local ``CacheManager`` backed by in-memory
Python dictionaries. Each cache entry stores a value, an optional expiry
deadline measured with ``time.monotonic()``, and an optional tag. A secondary
``tag -> set[key]`` index makes group eviction efficient.

The cache is explicit rather than global: callers instantiate ``CacheManager``
objects and pass them where needed. The implementation does not use locks, so a
single ``CacheManager`` instance is not thread-safe.
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

__all__ = [
    "CacheEntry",
    "CacheManager",
    "CacheStats",
    "make_args_hash",
    "make_cache_key",
]


@dataclass
class CacheStats:
    """
    Snapshot of live cache statistics returned by :meth:`CacheManager.stats`.

    The reported values describe the cache after expired entries have been
    purged. They are intended for diagnostics and monitoring rather than exact
    process-memory accounting.

    Attributes
    ----------
    hits : int
        Number of successful ``CacheManager.get`` lookups.
    misses : int
        Number of failed ``CacheManager.get`` lookups, including expired keys.
    size_bytes : int
        Shallow approximation of the live cache footprint in bytes. The
        estimate includes the top-level entry and tag dictionaries, their keys,
        and the immediate cached values, but does not traverse referenced
        objects recursively.
    num_entries : int
        Number of live entries remaining after expired values are purged. This
        reflects the keys that still participate in lookups and tag eviction.
    """

    hits: int
    misses: int
    size_bytes: int
    num_entries: int


@dataclass
class CacheEntry:
    """
    Single in-memory cache entry.

    Each entry stores a cached payload, an optional expiry deadline measured
    with ``time.monotonic()``, and an optional tag used for group invalidation.
    Tags let callers associate multiple cache keys with the same logical input
    such as one molecule, model configuration, or request family.

    Attributes
    ----------
    value : Any
        Cached payload returned by ``CacheManager.get``.
    expires_at : float | None
        ``time.monotonic()`` deadline when the key becomes stale. ``None`` means
        the entry does not expire automatically.
    tag : str | None
        Optional group label attached when calling ``cache.set(..., tag=...)``.
        All keys written with the same tag can be removed together with
        ``CacheManager.evict_tag``, which is useful when multiple cached values
        should be invalidated as one group.
    """

    value: Any
    expires_at: float | None
    tag: str | None


def make_args_hash(*args: Any, **kwargs: Any) -> str:
    """
    Generate a deterministic hash of arguments for cache keying.

    Tries JSON first for common types; falls back to pickle for complex objects.

    Examples
    --------
    >>> make_args_hash("CCO", az_model="USPTO")
    '6ad01e27a3a319962ad084787e060ab0fa0e661cc7d3e018e96747b06f7bacf7'
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

    Examples
    --------
    >>> make_cache_key("run_az", "CCO", az_model="USPTO", version=1)
    'v1:run_az:6ad01e27a3a319962ad084787e060ab0fa0e661cc7d3e018e96747b06f7bacf7'
    """
    if not namespace:
        raise ValueError("namespace must be a non-empty string")
    args_hash = make_args_hash(*args, **kwargs)
    return f"v{version}:{namespace}:{args_hash}"


class CacheManager:
    """
    Process-local in-memory cache manager with tag support and TTL.

    Each instance owns two in-memory indexes: ``_entries`` maps cache keys to
    :class:`CacheEntry` objects, and ``_tags`` maps each tag to the set of keys
    currently carrying that tag. ``get`` removes expired keys lazily,
    ``evict_tag`` removes every key associated with a tag, and ``stats`` first
    purges expired values so the reported counts reflect live entries only.

    The cache is process-local and not thread-safe. Reuse the same
    ``CacheManager`` instance only when callers intentionally want to share
    state.

    Examples
    --------
    >>> cache = CacheManager()
    >>> key = make_cache_key("call_llm", "CCO", model="gpt-5.4", version=1)
    >>> miss = object()
    >>> cache.get(key, default=miss) is miss
    True
    >>> cache.set(key, {"molecule": "CCO"}, expire=300, tag="molecule:CCO")
    >>> cache.get(key)
    {'molecule': 'CCO'}
    >>> cache.evict_tag("molecule:CCO")
    1
    """

    def __init__(self) -> None:
        """Initialize an empty in-memory cache."""
        self._entries: dict[str, CacheEntry] = {}
        self._tags: dict[str, set[str]] = {}
        self._hits = 0
        self._misses = 0
        self._log = logger.bind(component="cache")

    def purge_if_expired(self, key: str) -> bool:
        """
        Remove a key if its expiry deadline has passed.

        Parameters
        ----------
        key : str
            Cache key to inspect.

        Returns
        -------
        bool
            ``True`` when the key existed and was removed because it had
            expired, otherwise ``False``.
        """
        entry = self._entries.get(key)
        if entry is None:
            return False
        if entry.expires_at is None or entry.expires_at > time.monotonic():
            return False
        self.delete_key(key)
        self._log.debug("cache.expired", key=key)
        return True

    def delete_key(self, key: str) -> bool:
        """
        Remove a key from the cache and tag index if present.

        Parameters
        ----------
        key : str
            Cache key to remove.

        Returns
        -------
        bool
            ``True`` when an entry was removed, otherwise ``False``.
        """
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

    def purge_expired_entries(self) -> None:
        """
        Remove every expired entry currently stored in the cache.

        This is useful before inspecting cache size or exporting diagnostics.
        """
        for key in list(self._entries):
            self.purge_if_expired(key)

    def estimate_size_bytes(self) -> int:
        """
        Return a shallow approximation of the current in-memory cache size.

        The estimate includes the top-level dictionaries, keys, tag sets, and
        the immediate cached values. Referenced objects are not traversed
        recursively.
        """
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

        Examples
        --------
        >>> cache = CacheManager()
        >>> miss = object()
        >>> cache.get("missing", default=miss) is miss
        True
        """
        if self.purge_if_expired(key):
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
            Optional group label for later eviction via ``evict_tag``. Multiple
            keys may share the same tag.

        Examples
        --------
        >>> cache = CacheManager()
        >>> cache.set("demo", {"smiles": "CCO"}, expire=60, tag="molecule:CCO")
        """
        self.delete_key(key)
        expires_at = None if expire is None else time.monotonic() + expire
        self._entries[key] = CacheEntry(value=value, expires_at=expires_at, tag=tag)
        if tag is not None:
            self._tags.setdefault(tag, set()).add(key)
        self._log.debug("cache.set", key=key, expire=expire, tag=tag)

    def evict_tag(self, tag: str) -> int:
        """
        Remove all live entries with the given tag.

        Parameters
        ----------
        tag : str
            Group label identifying entries to remove. A single tag may be
            attached to multiple cache keys.

        Returns
        -------
        int
            Number of entries evicted.

        Examples
        --------
        >>> cache = CacheManager()
        >>> cache.set("a", 1, tag="batch:1")
        >>> cache.set("b", 2, tag="batch:1")
        >>> cache.evict_tag("batch:1")
        2
        """
        keys = list(self._tags.get(tag, set()))
        removed = 0
        for key in keys:
            removed += int(self.delete_key(key))
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
            A snapshot containing hit count, miss count, shallow byte estimate,
            and live entry count.

        Examples
        --------
        >>> cache = CacheManager()
        >>> cache.set("demo", 1)
        >>> stats = cache.stats()
        >>> (stats.hits, stats.misses, stats.num_entries)
        (0, 0, 1)
        """
        self.purge_expired_entries()
        return CacheStats(
            hits=self._hits,
            misses=self._misses,
            size_bytes=self.estimate_size_bytes(),
            num_entries=len(self._entries),
        )

    def close(self) -> None:
        """Release cache contents for compatibility with previous callers."""
        self.clear()
