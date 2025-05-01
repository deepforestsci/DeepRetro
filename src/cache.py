import diskcache as dc
import hashlib
import json
import os
import logging
import sqlite3 # Use built-in SQLite

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Configuration ---
CACHE_TYPE = os.getenv('CACHE_TYPE', 'diskcache').lower()
SQLITE_DB_PATH = os.getenv('SQLITE_CACHE_PATH', 'cache_api_new.db') # Path to SQLite DB file
SQLITE_CACHE_TABLE = 'results_cache' # Keep table name fixed

# Global instance for disk cache (if used)
disk_cache_instance = None

if CACHE_TYPE == 'diskcache':
    disk_cache_instance = dc.Cache('cache_api_new')
    logger.info("Using DiskCache for caching.")
elif CACHE_TYPE == 'sqlite':
    logger.info(f"Using SQLite for caching (DB File: {SQLITE_DB_PATH}).")
    # Initialize SQLite DB and table if not exists
    try:
        conn = sqlite3.connect(SQLITE_DB_PATH)
        cursor = conn.cursor()
        # Create table with compatible types (TEXT for most things is fine)
        # PRIMARY KEY ensures uniqueness and implicit index for cache_key
        cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS {SQLITE_CACHE_TABLE} (
                cache_key TEXT PRIMARY KEY,
                input_args TEXT,
                result TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        """)
        conn.commit()
        conn.close()
        logger.info(f"SQLite database and table '{SQLITE_CACHE_TABLE}' initialized successfully.")
    except sqlite3.Error as err:
        logger.error(f"SQLite initialization failed: {err}. Falling back to diskcache.")
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
    try:
        arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
    except TypeError as e:
        logger.error(f"Could not serialize arguments for caching function {func_name}: {e}")
        return f"{func_name}:serialization_error_{hash(str(args)+str(kwargs))}"
    args_hash = hashlib.md5(arg_string.encode('utf-8')).hexdigest()
    return f"{func_name}:{args_hash}"


# --- SQLite Helper Functions ---
def _get_sqlite_connection():
    """Establishes a connection to the SQLite database file."""
    try:
        # Consider adding timeout parameter if contention is expected
        return sqlite3.connect(SQLITE_DB_PATH)
    except sqlite3.Error as err:
        logger.error(f"SQLite connection error: {err}")
        return None

def _fetch_from_sqlite(cache_key):
    """Fetches result from SQLite cache."""
    conn = _get_sqlite_connection()
    if not conn: return None
    cursor = None
    try:
        cursor = conn.cursor()
        # Use placeholder (?) for security
        query = f"SELECT result FROM {SQLITE_CACHE_TABLE} WHERE cache_key = ?"
        cursor.execute(query, (cache_key,))
        result = cursor.fetchone()
        if result:
            logger.debug(f"SQLite Cache HIT for key: {cache_key}")
            # Result is stored as TEXT (JSON string)
            return json.loads(result[0])
        else:
            logger.debug(f"SQLite Cache MISS for key: {cache_key}")
            return None
    except sqlite3.Error as err:
        logger.error(f"SQLite fetch error for key {cache_key}: {err}")
        return None
    except json.JSONDecodeError as err:
         logger.error(f"SQLite JSON decode error for key {cache_key}: {err}")
         # Consider deleting the corrupt entry?
         return None # Treat as miss
    finally:
        # No need to close cursor explicitly with sqlite3 in basic cases
        if conn: conn.close()


def _store_in_sqlite(cache_key, input_args_json, result):
    """Stores result in SQLite cache."""
    conn = _get_sqlite_connection()
    if not conn: return
    cursor = None
    try:
        serialized_result = json.dumps(result) # Serialize result to JSON string
        cursor = conn.cursor()
        # Use SQLite's UPSERT capability (INSERT ON CONFLICT)
        query = f"""
            INSERT INTO {SQLITE_CACHE_TABLE} (cache_key, input_args, result)
            VALUES (?, ?, ?)
            ON CONFLICT(cache_key) DO UPDATE SET
                result = excluded.result,
                input_args = excluded.input_args;
        """
        # Using excluded.column refers to the values that would have been inserted
        cursor.execute(query, (cache_key, input_args_json, serialized_result))
        conn.commit()
        logger.debug(f"Stored result in SQLite for key: {cache_key}")
    except sqlite3.Error as err:
        logger.error(f"SQLite store error for key {cache_key}: {err}")
        if conn: conn.rollback() # Rollback on error
    except TypeError as e:
        logger.error(f"Could not serialize result for SQLite caching (key {cache_key}): {e}")
    finally:
        if conn: conn.close()

# --- Caching Decorator ---
def cache_results(func):
    """Decorator to cache the results of a function based on CACHE_TYPE."""

    def wrapper(*args, **kwargs):
        """Wrapper function to cache the results of the function."""
        cache_key = _generate_cache_key(func.__name__, *args, **kwargs)

        if CACHE_TYPE == 'sqlite':
            # --- SQLite Caching Logic ---
            cached_result = _fetch_from_sqlite(cache_key)
            if cached_result is not None:
                return cached_result
            else:
                # Cache miss, execute function
                result = func(*args, **kwargs)
                # Generate arg_string again for storing
                try:
                    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
                except TypeError:
                    arg_string = "{'error': 'args not serializable'}" # Placeholder
                _store_in_sqlite(cache_key, arg_string, result)
                return result

        else:
            # --- DiskCache Caching Logic (Original) ---
            if disk_cache_instance and cache_key in disk_cache_instance:
                logger.debug(f"DiskCache HIT for key: {cache_key}")
                cached_data = disk_cache_instance[cache_key]
                if isinstance(cached_data, dict) and 'result' in cached_data:
                    return cached_data['result']
                else:
                    logger.warning(f"DiskCache entry for {cache_key} has unexpected format. Returning raw value.")
                    return cached_data
            else:
                logger.debug(f"DiskCache MISS for key: {cache_key}")
                result = func(*args, **kwargs)
                try:
                    arg_string = json.dumps({'args': args, 'kwargs': kwargs}, sort_keys=True)
                except TypeError:
                    arg_string = "{'error': 'args not serializable'}" # Placeholder
                if disk_cache_instance:
                    disk_cache_instance[cache_key] = {'result': result, 'input_args': arg_string}
                return result

    return wrapper


# --- Cache Clearing Functions ---
def clear_entire_cache() -> None:
    """Function to clear the entire cache based on CACHE_TYPE."""
    if CACHE_TYPE == 'sqlite':
        conn = _get_sqlite_connection()
        if not conn: return
        cursor = None
        try:
            cursor = conn.cursor()
            # DELETE is standard SQL and works fine here
            query = f"DELETE FROM {SQLITE_CACHE_TABLE}"
            cursor.execute(query)
            conn.commit()
            # Optional: Vacuum to reclaim disk space after large delete
            # conn.execute("VACUUM;")
            # conn.commit()
            logger.info("Cleared entire SQLite cache table.")
        except sqlite3.Error as err:
            logger.error(f"SQLite clear cache error: {err}")
            if conn: conn.rollback()
        finally:
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
    if CACHE_TYPE == 'sqlite':
        conn = _get_sqlite_connection()
        if not conn: return
        cursor = None
        try:
            cursor = conn.cursor()
            # Use LIKE to find the molecule within the JSON string (same logic as before)
            like_pattern = f'%"{molecule}"%'
            query = f"DELETE FROM {SQLITE_CACHE_TABLE} WHERE input_args LIKE ?"
            cursor.execute(query, (like_pattern,))
            deleted_count = cursor.rowcount # Get count of deleted rows
            conn.commit()
            logger.info(f"Cleared {deleted_count} SQLite cache entries matching molecule: {molecule}")
        except sqlite3.Error as err:
            logger.error(f"SQLite clear_cache_for_molecule error: {err}")
            if conn: conn.rollback()
        finally:
            if conn: conn.close()

    else:
        # DiskCache Logic (No changes needed here)
        if not disk_cache_instance: return
        keys_to_delete = []
        deleted_count = 0
        try:
            for key in disk_cache_instance.iterkeys():
                data = disk_cache_instance.get(key)
                if data and isinstance(data, dict) and 'input_args' in data:
                    if molecule in data['input_args']:
                        keys_to_delete.append(key)
                elif isinstance(data, str) and molecule in data:
                     keys_to_delete.append(key)

            for k in keys_to_delete:
                try:
                    del disk_cache_instance[k]
                    deleted_count += 1
                except KeyError:
                    pass
            logger.info(f"Cleared {deleted_count} DiskCache entries matching molecule: {molecule}")
        except Exception as e:
             logger.error(f"Error during DiskCache clear_cache_for_molecule: {e}")

# Example usage remains conceptually the same, just set CACHE_TYPE=sqlite
# and potentially SQLITE_CACHE_PATH if you don't want the default 'cache_api_new.db'
