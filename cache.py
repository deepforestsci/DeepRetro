import diskcache as dc

# Create a disk cache
cache = dc.Cache('cache_api')


def clear_cache():
    """Clear the cache"""
    cache.clear()


def cache_results(func):
    """Decorator to cache results using diskcache"""

    def wrapper(*args, **kwargs):
        cache_key = func.__name__ + "_" + str(args) + str(
            kwargs)  # Unique key with function name
        if cache_key in cache:
            return cache[cache_key]
        else:
            result = func(*args, **kwargs)
            cache[cache_key] = result
            return result

    return wrapper
