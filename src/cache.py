import diskcache as dc
import hashlib
import json
import os
import logging

# Try importing mysql.connector, handle if not installed for MySQL caching
try:
    import mysql.connector
    MYSQL_AVAILABLE = True
except ImportError:
    MYSQL_AVAILABLE = False

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Configuration ---
CACHE_TYPE = os.getenv('CACHE_TYPE', 'diskcache').lower()
MYSQL_CACHE_TABLE = 'results_cache' # Keep table name fixed for simplicity

# Global instance for disk cache (if used)
disk_cache_instance = None
mysql_config = {}

if CACHE_TYPE == 'diskcache':
    disk_cache_instance = dc.Cache('cache_api_new')
    logger.info("Using DiskCache for caching.")
elif CACHE_TYPE == 'mysql':
    if not MYSQL_AVAILABLE:
        logger.error("CACHE_TYPE set to mysql, but mysql-connector-python is not installed. Falling back to diskcache.")
        CACHE_TYPE = 'diskcache' # Fallback
        disk_cache_instance = dc.Cache('cache_api_new')
    else:
        mysql_config = {
            'host': os.getenv('MYSQL_HOST'),
            'user': os.getenv('MYSQL_USER'),
            'password': os.getenv('MYSQL_PASSWORD'),
            'database': os.getenv('MYSQL_DATABASE')
        }
        # Basic check for required MySQL config
        if not all(mysql_config.values()):
            logger.error("MYSQL cache type selected, but required environment variables (MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE) are missing. Falling back to diskcache.")
            CACHE_TYPE = 'diskcache' # Fallback
            disk_cache_instance = dc.Cache('cache_api_new')
        else:
            logger.info(f"Using MySQL for caching (Host: {mysql_config['host']}, DB: {mysql_config['database']}).")
            # Basic connection test (optional but recommended)
            try:
                conn = mysql.connector.connect(**mysql_config)
                conn.close()
                logger.info("MySQL connection successful.")
            except mysql.connector.Error as err:
                logger.error(f"MySQL connection failed: {err}. Falling back to diskcache.")
                CACHE_TYPE = 'diskcache' # Fallback
                disk_cache_instance = dc.Cache('cache_api_new')
else:
    logger.warning(f"Unknown CACHE_TYPE '{CACHE_TYPE}'. Defaulting to diskcache.")
    CACHE_TYPE = 'diskcache'
    disk_cache_instance = dc.Cache('cache_api_new')


# --- Cache Key Generation ---
def _generate_cache_key(func_name, *args, **kwargs):
    """
    Generates a unique cache key based on the function name and
    a hash of the arguments.
    """
    # Ensure args are consistently serializable (handle potential non-serializable args if necessary)
    try:
        arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
    except TypeError as e:
        logger.error(f"Could not serialize arguments for caching function {func_name}: {e}")
        # Decide how to handle: raise error, return specific key, or skip caching?
        # For now, let's create a key indicating the error to avoid caching potentially bad state
        return f"{func_name}:serialization_error_{hash(str(args)+str(kwargs))}"
    args_hash = hashlib.md5(arg_string.encode('utf-8')).hexdigest()
    # Using MD5 hash directly as part of the key, ensure it fits DB constraints if needed
    return f"{func_name}:{args_hash}"


# --- MySQL Helper Functions (Basic) ---
def _get_mysql_connection():
    """Establishes a new MySQL connection."""
    try:
        return mysql.connector.connect(**mysql_config)
    except mysql.connector.Error as err:
        logger.error(f"MySQL connection error: {err}")
        return None

def _fetch_from_mysql(cache_key):
    """Fetches result from MySQL cache."""
    conn = _get_mysql_connection()
    if not conn: return None
    cursor = None
    try:
        cursor = conn.cursor()
        query = f"SELECT result FROM {MYSQL_CACHE_TABLE} WHERE cache_key = %s"
        cursor.execute(query, (cache_key,))
        result = cursor.fetchone()
        if result:
            logger.debug(f"MySQL Cache HIT for key: {cache_key}")
            return json.loads(result[0]) # Deserialize JSON result
        else:
            logger.debug(f"MySQL Cache MISS for key: {cache_key}")
            return None
    except mysql.connector.Error as err:
        logger.error(f"MySQL fetch error for key {cache_key}: {err}")
        return None
    except json.JSONDecodeError as err:
         logger.error(f"MySQL JSON decode error for key {cache_key}: {err}")
         # Consider deleting the corrupt entry?
         return None # Treat as miss
    finally:
        if cursor: cursor.close()
        if conn: conn.close()


def _store_in_mysql(cache_key, input_args_json, result):
    """Stores result in MySQL cache."""
    conn = _get_mysql_connection()
    if not conn: return
    cursor = None
    try:
        serialized_result = json.dumps(result) # Serialize result to JSON
        cursor = conn.cursor()
        # Use INSERT IGNORE or INSERT ... ON DUPLICATE KEY UPDATE for robustness
        query = f"""
            INSERT INTO {MYSQL_CACHE_TABLE} (cache_key, input_args, result)
            VALUES (%s, %s, %s)
            ON DUPLICATE KEY UPDATE result = VALUES(result), input_args = VALUES(input_args)
        """
        cursor.execute(query, (cache_key, input_args_json, serialized_result))
        conn.commit()
        logger.debug(f"Stored result in MySQL for key: {cache_key}")
    except mysql.connector.Error as err:
        logger.error(f"MySQL store error for key {cache_key}: {err}")
        if conn: conn.rollback() # Rollback on error
    except TypeError as e:
        logger.error(f"Could not serialize result for MySQL caching (key {cache_key}): {e}")
    finally:
        if cursor: cursor.close()
        if conn: conn.close()

# --- Caching Decorator ---
def cache_results(func):
    """Decorator to cache the results of a function based on CACHE_TYPE."""

    def wrapper(*args, **kwargs):
        """Wrapper function to cache the results of the function."""
        cache_key = _generate_cache_key(func.__name__, *args, **kwargs)

        if CACHE_TYPE == 'mysql':
            # --- MySQL Caching Logic ---
            cached_result = _fetch_from_mysql(cache_key)
            if cached_result is not None:
                return cached_result
            else:
                # Cache miss, execute function
                result = func(*args, **kwargs)
                # Generate arg_string again for storing (could optimize)
                try:
                    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
                except TypeError:
                    arg_string = "{'error': 'args not serializable'}" # Placeholder
                _store_in_mysql(cache_key, arg_string, result)
                return result

        else:
            # --- DiskCache Caching Logic (Original) ---
            if cache_key in disk_cache_instance:
                logger.debug(f"DiskCache HIT for key: {cache_key}")
                # Ensure compatibility with how data might be stored (dict vs direct value)
                cached_data = disk_cache_instance[cache_key]
                if isinstance(cached_data, dict) and 'result' in cached_data:
                    return cached_data['result']
                else:
                    # If stored directly or old format, return it
                    # Might need adjustment based on how things were stored previously
                    logger.warning(f"DiskCache entry for {cache_key} has unexpected format. Returning raw value.")
                    return cached_data
            else:
                logger.debug(f"DiskCache MISS for key: {cache_key}")
                result = func(*args, **kwargs)
                # Store input args along with result for consistency and potential clearing
                try:
                    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
                except TypeError:
                    arg_string = "{'error': 'args not serializable'}" # Placeholder
                disk_cache_instance[cache_key] = {'result': result, 'input_args': arg_string}
                return result

    return wrapper


# --- Cache Clearing Functions ---
def clear_entire_cache() -> None:
    """Function to clear the entire cache based on CACHE_TYPE."""
    if CACHE_TYPE == 'mysql':
        conn = _get_mysql_connection()
        if not conn: return
        cursor = None
        try:
            cursor = conn.cursor()
            # TRUNCATE is faster but requires more privileges, DELETE is safer
            query = f"DELETE FROM {MYSQL_CACHE_TABLE}"
            cursor.execute(query)
            conn.commit()
            logger.info("Cleared entire MySQL cache table.")
        except mysql.connector.Error as err:
            logger.error(f"MySQL clear cache error: {err}")
            if conn: conn.rollback()
        finally:
            if cursor: cursor.close()
            if conn: conn.close()
    else:
        # DiskCache Logic
        if disk_cache_instance:
            disk_cache_instance.clear()
            logger.info("Cleared entire DiskCache.")


def clear_cache_for_molecule(molecule: str):
    """
    Clear cache entries specifically related to a given molecule string, based on CACHE_TYPE.
    Requires 'input_args' to be stored as a JSON string in the cache.
    """
    if CACHE_TYPE == 'mysql':
        conn = _get_mysql_connection()
        if not conn: return
        cursor = None
        try:
            cursor = conn.cursor()
            # Use LIKE to find the molecule within the JSON string
            # This is inefficient but simple for this prototype
            # Assumes molecule SMILES is stored like: '..."smiles": "...",...' or in args list
            like_pattern = f'%"{molecule}"%' # Basic check, might need refinement
            query = f"DELETE FROM {MYSQL_CACHE_TABLE} WHERE input_args LIKE %s"
            cursor.execute(query, (like_pattern,))
            deleted_count = cursor.rowcount
            conn.commit()
            logger.info(f"Cleared {deleted_count} MySQL cache entries matching molecule: {molecule}")
        except mysql.connector.Error as err:
            logger.error(f"MySQL clear_cache_for_molecule error: {err}")
            if conn: conn.rollback()
        finally:
            if cursor: cursor.close()
            if conn: conn.close()

    else:
        # DiskCache Logic
        if not disk_cache_instance: return
        keys_to_delete = []
        deleted_count = 0
        try:
            # Iterating keys can be slow on large caches
            for key in disk_cache_instance.iterkeys():
                data = disk_cache_instance.get(key)
                # Check if data is stored in the expected dictionary format
                if data and isinstance(data, dict) and 'input_args' in data:
                    # Check if the molecule string appears within the input_args JSON string
                    if molecule in data['input_args']:
                        keys_to_delete.append(key)
                elif isinstance(data, str) and molecule in data: # Fallback for potentially old string format?
                     keys_to_delete.append(key)

            for k in keys_to_delete:
                try:
                    del disk_cache_instance[k]
                    deleted_count += 1
                except KeyError:
                    pass # Key already deleted, ignore
            logger.info(f"Cleared {deleted_count} DiskCache entries matching molecule: {molecule}")
        except Exception as e:
             logger.error(f"Error during DiskCache clear_cache_for_molecule: {e}")

# Example of how to use (Optional - for testing)
# @cache_results
# def expensive_calculation(a, b):
#     print(f"Performing expensive calculation for {a}, {b}...")
#     import time
#     time.sleep(2)
#     return a + b

# if __name__ == "__main__":
#     # Set environment variables before running for MySQL:
#     # export CACHE_TYPE=mysql
#     # export MYSQL_HOST=localhost
#     # export MYSQL_USER=youruser
#     # export MYSQL_PASSWORD=yourpass
#     # export MYSQL_DATABASE=yourdb
#     # Ensure the table 'results_cache' exists in 'yourdb'
#     print("First call (expect delay):")
#     print(expensive_calculation(5, 3))
#     print("Second call (expect cache hit):")
#     print(expensive_calculation(5, 3))

#     print("Clearing cache for molecule 'test_mol'")
#     # Simulate storing something with a 'molecule' arg if needed for testing clear
#     # _store_in_mysql or use disk_cache_instance directly depending on CACHE_TYPE
#     clear_cache_for_molecule("test_mol")

#     print("Clearing entire cache")
#     clear_entire_cache()
