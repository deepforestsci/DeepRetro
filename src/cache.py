"""Caching utilities for the DeepRetro application.

This module provides a decorator and functions to cache results of expensive
computations using the `diskcache` library. It helps in speeding up repeated
analyses by storing and retrieving previously computed results.

Attributes:
    cache (diskcache.Cache): The global cache object instance used for storing
        and retrieving cached items.
"""
import diskcache as dc
import hashlib
import json
import functools

cache = dc.Cache('cache_api_new')


def _generate_cache_key(func_name, *args, **kwargs):
    """Generates a unique cache key based on the function name and its arguments.

    This is an internal helper function. The key is formed by combining the
    function's name with an MD5 hash of its serialized arguments to ensure
    uniqueness.

    Args:
        func_name (str): The name of the function being called.
        *args: Positional arguments passed to the function.
        **kwargs: Keyword arguments passed to the function.

    Returns:
        str: A unique string cache key.
    """
    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
    args_hash = hashlib.md5(arg_string.encode('utf-8')).hexdigest()
    return f"{func_name}:{args_hash}"


def cache_results(func):
    """Decorator to cache the results of a function using diskcache.

    If the decorated function is called with the same arguments (as determined
    by the generated cache key), the cached result is returned. Otherwise,
    the function is executed, and its result is stored in the cache for
    future calls.

    Args:
        func (callable): The function to be decorated.

    Returns:
        callable: The wrapped function with caching capabilities.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """Wrapper that implements the caching logic.

        Checks if a result for the given function and arguments already exists
        in the cache. If so, it returns the cached result. Otherwise, it
        executes the original function, caches its result, and then returns
        the result.

        Args:
            *args: Positional arguments for the wrapped function.
            **kwargs: Keyword arguments for the wrapped function.

        Returns:
            Any: The result of the wrapped function, which could be retrieved
                 from the cache or by executing the function.
        """
        arg_string = json.dumps({
            'args': args,
            'kwargs': kwargs
        },
                                sort_keys=True)
        cache_key = _generate_cache_key(func.__name__, *args, **kwargs)
        if cache_key in cache:
            return cache[cache_key]['result']
        else:
            result = func(*args, **kwargs)
            cache[cache_key] = {'result': result, 'input_args': arg_string}
            return result

    return wrapper


def clear_entire_cache() -> None:
    """Clears all items from the global disk cache.

    This function will remove all persisted cache entries managed by the
    `cache` object defined in this module.
    """
    cache.clear()


def clear_cache_for_molecule(molecule):
    """Clears cache entries potentially related to a specific molecule string.

    This function iterates through the cache and removes entries where the
    provided `molecule` string (typically a SMILES string) is found within
    the string representation of the cached function's input arguments.

    Note:
        The matching is based on a substring check within the serialized
        arguments, which might not be perfectly precise for all use cases
        but serves as a targeted way to invalidate cache entries related
        to a particular molecule.

    Args:
        molecule (str): The molecule string (e.g., SMILES) to search for
                        in the `input_args` of cached items.
    """
    keys_to_delete = []
    for key in cache.iterkeys():
        data = cache.get(key)
        if data and 'input_args' in data:
            # Just check if the string representation of molecule is in the JSON
            if molecule in data['input_args']:
                keys_to_delete.append(key)
    for k in keys_to_delete:
        del cache[k]
