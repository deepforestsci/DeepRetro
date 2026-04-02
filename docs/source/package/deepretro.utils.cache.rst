deepretro.utils.cache
=====================

Explicit in-memory caching utilities for expensive library operations.

Overview
--------

The cache module provides:

- **CacheManager**: process-local cache with tag support, TTL, and statistics
- **make_cache_key()**: deterministic cache keys for explicit cache lookups

The cache only lives for the current Python process. If the process exits,
the cached values are discarded.

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

   # Clear this cache instance
   cache.clear()

   # Inspect statistics
   stats = cache.stats()
   print(stats.hits, stats.misses, stats.num_entries)

API
---

.. automodule:: deepretro.utils.cache
   :members:
   :undoc-members:
