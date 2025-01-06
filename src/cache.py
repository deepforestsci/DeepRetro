import diskcache as dc
import hashlib
import json

cache = dc.Cache('cache_api')


def _generate_cache_key(func_name, *args, **kwargs):
    """
    Generates a unique cache key based on the function name and
    a hash of the arguments.
    """
    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
    args_hash = hashlib.md5(arg_string.encode('utf-8')).hexdigest()
    return f"{func_name}:{args_hash}"


def cache_results(func):
    """Decorator to cache the results of a function.

    Parameters
    ----------
    func : _type_
        Function to be decorated
    """

    def wrapper(*args, **kwargs):
        """Wrapper function to cache the results of the function.

        Returns
        -------
        _type_
            Result of the function
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
    """Function to clear the entire cache.
    """
    cache.clear()


def clear_cache_for_molecule(molecule):
    """
    Clear cache entries specifically related to a given molecule string.

    Parameters
    ----------
    molecule : str
        SMILES string of the molecule
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
