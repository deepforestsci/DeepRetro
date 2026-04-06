deepretro.utils.cache
=====================

Explicit in-memory caching utilities for expensive library operations.

Overview
--------

The cache module provides:

- **CacheEntry**: public record describing a cached value, its expiry deadline,
  and optional eviction tag
- **CacheManager**: process-local cache with tag support, TTL, and statistics
- **make_args_hash()**: deterministic argument hashing used by cache keys
- **make_cache_key()**: deterministic cache keys for explicit cache lookups

The cache only lives for the current Python process. If the process exits,
the cached values are discarded. ``CacheManager`` uses plain Python data
structures without locking, so one instance should not be shared across threads
unless the caller adds synchronization.

Key Format
----------

``make_cache_key()`` returns keys in the form
``v<version>:<namespace>:<64-char-sha256>``. For example:

.. code-block:: python

   from deepretro.utils.cache import make_args_hash, make_cache_key

   print(make_args_hash("CCO", az_model="USPTO"))
   # 6ad01e27a3a319962ad084787e060ab0fa0e661cc7d3e018e96747b06f7bacf7

   print(make_cache_key("run_az", "CCO", az_model="USPTO", version=1))
   # v1:run_az:6ad01e27a3a319962ad084787e060ab0fa0e661cc7d3e018e96747b06f7bacf7

Algorithm Notes
---------------

``make_args_hash()`` serializes ``args`` and ``kwargs`` into a stable payload,
tries JSON first for common serializable values, and falls back to pickle for
complex Python objects. ``CacheManager`` stores each key in ``_entries`` and
maintains a secondary ``tag -> set[key]`` mapping so one tag can evict multiple
keys with ``evict_tag()``.

Usage
-----

Create and pass cache objects explicitly:

.. code-block:: python

   from deepretro.utils.cache import CacheManager, make_cache_key

   cache = CacheManager()
   key = make_cache_key("call_llm", "CCO", model="claude-opus-4-6", version=1)
   cache_miss = object()
   result = cache.get(key, default=cache_miss)

   if result is cache_miss:
       result = {"molecule": "CCO", "model": "claude-opus-4-6"}
       cache.set(key, result, expire=3600, tag="CCO")

   # Evict all entries for a molecule from this cache instance
   cache.evict_tag("CCO")

   # Inspect helper methods directly during tests or diagnostics
   cache.purge_expired_entries()
   print(cache.estimate_size_bytes())

   # Clear this cache instance
   cache.clear()

   # Inspect statistics
   stats = cache.stats()
   print(stats.hits, stats.misses, stats.num_entries)

Tag Semantics
-------------

A tag is an arbitrary group label attached to one or more cache keys when
calling ``cache.set(..., tag="...")``. Tags are useful when many cached values
should be invalidated together, such as all results derived from one molecule or
all outputs from one model configuration.

API
---

.. automodule:: deepretro.utils.cache
   :members:
