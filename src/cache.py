"""Caching utilities for DeepRetro API functions.

Provides decorators and functions to cache and manage results for retrosynthesis and related computations.
"""
import diskcache as dc
import hashlib
import json
import functools

cache = dc.Cache('cache_api_new')

def _generate_cache_key(func_name, *args, **kwargs):
    """
    Generates a unique cache key based on the function name and
    a hash of the arguments.

    Parameters
    ----------
    func_name : str
        Name of the function.
    *args : tuple
        Positional arguments to the function.
    **kwargs : dict
        Keyword arguments to the function.

    Returns
    -------
    str
        A unique cache key string.
    """
    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
    args_hash = hashlib.md5(arg_string.encode('utf-8')).hexdigest()
    return f"{func_name}:{args_hash}"

def cache_results(func):
    """
    Decorator to cache the results of a function.

    Parameters
    ----------
    func : callable
        Function to be decorated.

    Returns
    -------
    callable
        Decorated function with caching.

    Examples
    --------
    >>> @cache_results
    ... def add(a, b):
    ...     return a + b
    >>> add(2, 3)  # doctest: +SKIP
    5
    >>> add(2, 3)  # This call will use the cache  # doctest: +SKIP
    5
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """
        Wrapper function to cache the results of the function.

        Returns
        -------
        Any
            Result of the function.
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
    """
    Function to clear the entire cache.

    Examples
    --------
    >>> clear_entire_cache()  # doctest: +SKIP
    # All cache entries are removed.
    """
    cache.clear()

def clear_cache_for_molecule(molecule):
    """
    Clear cache entries specifically related to a given molecule string.

    Parameters
    ----------
    molecule : str
        SMILES string of the molecule

    Examples
    --------
    >>> @cache_results
    ... def foo(mol):
    ...     return mol[::-1]
    >>> foo('CCO')  # doctest: +SKIP
    'OCC'
    >>> clear_cache_for_molecule('CCO')  # doctest: +SKIP
    # All cache entries for 'CCO' are removed.
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
