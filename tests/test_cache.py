import pytest
from unittest import mock
import json
import hashlib

# Module to be tested
from src import cache as cache_module

@pytest.fixture
def mock_diskcache_cache():
    """Mocks the cache object in the cache_module."""
    # Create a dictionary-like mock for the cache
    # We can use a real dictionary to simulate the cache's behavior for get, set, in, del
    mock_cache_storage = {}

    mock_cache = mock.Mock()
    mock_cache.__contains__ = mock.Mock(side_effect=lambda key: key in mock_cache_storage)
    mock_cache.__getitem__ = mock.Mock(side_effect=lambda key: mock_cache_storage[key])
    mock_cache.__setitem__ = mock.Mock(side_effect=lambda key, value: mock_cache_storage.setdefault(key, value))
    mock_cache.__delitem__ = mock.Mock(side_effect=lambda key: mock_cache_storage.pop(key, None))
    mock_cache.clear = mock.Mock(side_effect=mock_cache_storage.clear)
    # For iterkeys, we need to return an iterator over the keys of our mock storage
    mock_cache.iterkeys = mock.Mock(side_effect=lambda: iter(list(mock_cache_storage.keys())))
    mock_cache.get = mock.Mock(side_effect=lambda key, default=None: mock_cache_storage.get(key, default))

    with mock.patch('src.cache.cache', mock_cache):
        yield mock_cache, mock_cache_storage # Yield both for assertions

# Tests will be added below
def test_generate_cache_key_format_and_consistency():
    """Test _generate_cache_key for correct format and consistent hashing."""
    func_name = "test_func"
    args1 = (1, "hello")
    kwargs1 = {"a": True, "b": None}

    key1 = cache_module._generate_cache_key(func_name, *args1, **kwargs1)
    
    assert key1.startswith(f"{func_name}:")
    hash_part1 = key1.split(":")[1]
    assert len(hash_part1) == 32 # MD5 hexdigest length

    # Calling again with same args should produce the same key
    key2 = cache_module._generate_cache_key(func_name, *args1, **kwargs1)
    assert key1 == key2

def test_generate_cache_key_arg_sensitivity():
    """Test that _generate_cache_key is sensitive to argument changes."""
    func_name = "test_func"
    
    key_base = cache_module._generate_cache_key(func_name, 1, x=10)
    
    # Change positional arg
    key_diff_arg = cache_module._generate_cache_key(func_name, 2, x=10)
    assert key_base != key_diff_arg

    # Change kwarg value
    key_diff_kwarg_val = cache_module._generate_cache_key(func_name, 1, x=20)
    assert key_base != key_diff_kwarg_val

    # Change kwarg key
    key_diff_kwarg_key = cache_module._generate_cache_key(func_name, 1, y=10)
    assert key_base != key_diff_kwarg_key

    # Change function name
    key_diff_func_name = cache_module._generate_cache_key("another_func", 1, x=10)
    assert key_base != key_diff_func_name

    # Order of kwargs should not matter due to sort_keys=True in json.dumps
    key_kwarg_order1 = cache_module._generate_cache_key(func_name, val=1, name="test")
    key_kwarg_order2 = cache_module._generate_cache_key(func_name, name="test", val=1)
    assert key_kwarg_order1 == key_kwarg_order2

# --- Tests for cache_results decorator ---

# A simple function to be decorated for testing purposes
@cache_module.cache_results
def example_function_to_cache(arg1, kwarg1=None):
    """A sample function whose results we want to cache."""
    # print(f"Executing example_function_to_cache with {arg1}, {kwarg1}") # For debugging if needed
    return f"result:{arg1}-{kwarg1}"

# Mock the original function for spying on calls
# We need to wrap it carefully because it's already decorated by the time we access it in the module
# A better way for spying is to mock the function *before* decoration if possible,
# or test the decorator with a function defined *inside* the test scope.
# For simplicity here, we'll assume example_function_to_cache is primarily for testing the decorator.

# We will use a callable mock to act as the "original" function inside the test
# to verify if it gets called or not.

def test_cache_results_decorator_cache_miss_and_hit(mock_diskcache_cache):
    """Test cache_results decorator for cache miss and subsequent cache hit."""
    mock_cache, mock_cache_storage = mock_diskcache_cache
    
    # This is our stand-in for the actual function logic, so we can count calls
    original_function_mock = mock.Mock(wraps=lambda arg1, kwarg1=None: f"result:{arg1}-{kwarg1}")

    # Apply the decorator to our mock *within the test*
    # This is a common pattern for testing decorators: apply it to a controlled mock.
    @cache_module.cache_results
    def decorated_test_func(arg1, kwarg1=None):
        return original_function_mock(arg1, kwarg1=kwarg1)

    # --- Cache Miss ---
    arg1_val = "data1"
    kwarg1_val = "optionA"
    expected_result = f"result:{arg1_val}-{kwarg1_val}"

    # First call: should be a cache miss
    result1 = decorated_test_func(arg1_val, kwarg1=kwarg1_val)
    
    assert result1 == expected_result
    original_function_mock.assert_called_once_with(arg1_val, kwarg1=kwarg1_val) # Original func was called

    # Verify item was added to cache
    # Construct the expected key (can use _generate_cache_key or replicate its logic if it's simple)
    # For robustness, let's use the actual _generate_cache_key from the module
    expected_key = cache_module._generate_cache_key(decorated_test_func.__name__, arg1_val, kwarg1=kwarg1_val)
    assert expected_key in mock_cache_storage
    assert mock_cache_storage[expected_key]['result'] == expected_result
    
    # Check that the cache mock's __setitem__ was called
    mock_cache.__setitem__.assert_called_once()
    args_set, _ = mock_cache.__setitem__.call_args
    assert args_set[0] == expected_key # key
    assert args_set[1]['result'] == expected_result # value

    # --- Cache Hit ---
    original_function_mock.reset_mock() # Reset call count for the next check
    mock_cache.__getitem__.reset_mock() # Reset mock for __getitem__ for precise assertion

    # Second call with same args: should be a cache hit
    result2 = decorated_test_func(arg1_val, kwarg1=kwarg1_val)
    
    assert result2 == expected_result
    original_function_mock.assert_not_called() # Original func NOT called again
    
    # Verify item was fetched from cache
    mock_cache.__getitem__.assert_called_once_with(expected_key)

def test_cache_results_decorator_different_args(mock_diskcache_cache):
    """Test that different arguments result in different cache entries and calls."""
    mock_cache, mock_cache_storage = mock_diskcache_cache
    original_function_mock = mock.Mock(side_effect=lambda arg, **kwargs: f"res:{arg}:{kwargs.get('k')}")

    @cache_module.cache_results
    def decorated_test_func(arg, k=None):
        return original_function_mock(arg, k=k)

    # Call 1
    res1 = decorated_test_func("val1", k="opt1")
    assert res1 == "res:val1:opt1"
    original_function_mock.assert_called_once_with("val1", k="opt1")
    key1 = cache_module._generate_cache_key(decorated_test_func.__name__, "val1", k="opt1")
    assert key1 in mock_cache_storage

    original_function_mock.reset_mock()

    # Call 2 (different args)
    res2 = decorated_test_func("val2", k="opt2")
    assert res2 == "res:val2:opt2"
    original_function_mock.assert_called_once_with("val2", k="opt2") # Called again due to different args
    key2 = cache_module._generate_cache_key(decorated_test_func.__name__, "val2", k="opt2")
    assert key2 in mock_cache_storage
    assert key1 != key2

    original_function_mock.reset_mock()

    # Call 3 (repeat of call 1 - should be a hit)
    res3 = decorated_test_func("val1", k="opt1")
    assert res3 == "res:val1:opt1"
    original_function_mock.assert_not_called()

    # Call 4 (repeat of call 2 - should be a hit)
    res4 = decorated_test_func("val2", k="opt2")
    assert res4 == "res:val2:opt2"
    original_function_mock.assert_not_called()

# --- Tests for clear_entire_cache --- 
def test_clear_entire_cache(mock_diskcache_cache):
    """Test that clear_entire_cache calls cache.clear()."""
    mock_cache, mock_cache_storage = mock_diskcache_cache
    
    # Add something to the mock storage to ensure clear does something
    mock_cache_storage["some_key"] = "some_value"
    assert len(mock_cache_storage) == 1

    cache_module.clear_entire_cache()

    mock_cache.clear.assert_called_once()
    assert len(mock_cache_storage) == 0 # Verify our mock storage was cleared

# --- Tests for clear_cache_for_molecule ---
def test_clear_cache_for_molecule(mock_diskcache_cache):
    """Test clearing cache entries related to a specific molecule."""
    mock_cache, mock_cache_storage = mock_diskcache_cache

    target_molecule = "CCO" # Ethanol
    other_molecule = "CCC"  # Propane

    # Populate mock_cache_storage with some entries
    # Key generation logic for testing (simplified from _generate_cache_key for clarity in test setup)
    key1_target = "func1:key_for_CCO_and_other_stuff"
    mock_cache_storage[key1_target] = {
        "result": "res1", 
        "input_args": json.dumps({'args': (target_molecule, "paramX"), 'kwargs': {'temp': 25}})
    }

    key2_other = "func2:key_for_CCC_only"
    mock_cache_storage[key2_other] = {
        "result": "res2", 
        "input_args": json.dumps({'args': (other_molecule,), 'kwargs': {}})
    }
    
    key3_target_and_other = "func3:key_for_CCO_and_CCC"
    mock_cache_storage[key3_target_and_other] = {
        "result": "res3", 
        "input_args": json.dumps({'args': (target_molecule, other_molecule), 'kwargs': {'time': 60}})
    }

    key4_no_input_args = "func4:no_input_args_field"
    mock_cache_storage[key4_no_input_args] = {"result": "res4"} # Missing 'input_args'
    
    key5_malformed_input_args = "func5:malformed_input_args"
    mock_cache_storage[key5_malformed_input_args] = {"result": "res5", "input_args": "not_json_but_contains_CCO"}

    # Call the function to be tested
    cache_module.clear_cache_for_molecule(target_molecule)

    # Assertions
    # Check which keys were deleted from our mock_cache_storage
    # __delitem__ on the mock_cache would have been called
    
    # key1_target should be deleted
    assert key1_target not in mock_cache_storage
    # key3_target_and_other should be deleted
    assert key3_target_and_other not in mock_cache_storage
    # key5_malformed_input_args should be deleted because "CCO" is in the string
    assert key5_malformed_input_args not in mock_cache_storage

    # key2_other should NOT be deleted
    assert key2_other in mock_cache_storage
    # key4_no_input_args should NOT be deleted (and shouldn't cause an error)
    assert key4_no_input_args in mock_cache_storage

    # Verify that mock_cache.__delitem__ was called for the correct keys
    # The order of deletion might not be guaranteed if iterkeys() doesn't guarantee order (though dict keys are ordered in modern Python)
    # So, check the set of called arguments for __delitem__
    expected_deleted_keys = {key1_target, key3_target_and_other, key5_malformed_input_args}
    actual_deleted_keys = set()
    for call_args in mock_cache.__delitem__.call_args_list:
        actual_deleted_keys.add(call_args[0][0]) # The first arg of the call is the key
    
    assert actual_deleted_keys == expected_deleted_keys
    assert mock_cache.__delitem__.call_count == len(expected_deleted_keys)

def test_clear_cache_for_molecule_no_matches(mock_diskcache_cache):
    """Test clear_cache_for_molecule when no entries match the molecule."""
    mock_cache, mock_cache_storage = mock_diskcache_cache
    
    mock_cache_storage["key1"] = {"input_args": json.dumps({"args": ("CCC",)}), "result": "data"}
    mock_cache_storage["key2"] = {"input_args": json.dumps({"args": ("OOO",)}), "result": "data2"}
    
    initial_keys = set(mock_cache_storage.keys())

    cache_module.clear_cache_for_molecule("XXX") # Molecule not in any entry
    
    assert not mock_cache.__delitem__.called # No deletions should occur
    assert set(mock_cache_storage.keys()) == initial_keys # Storage should be unchanged

def test_clear_cache_for_molecule_empty_cache(mock_diskcache_cache):
    """Test clear_cache_for_molecule on an empty cache."""
    mock_cache, _ = mock_diskcache_cache # We don't need mock_cache_storage here
    
    cache_module.clear_cache_for_molecule("CCO")
    
    # iterkeys will be called even on an empty cache, it will just return an empty iterator.
    mock_cache.iterkeys.assert_called_once() 
    
    # Since the cache is empty, get and __delitem__ should not have been called.
    mock_cache.get.assert_not_called()
    mock_cache.__delitem__.assert_not_called() 