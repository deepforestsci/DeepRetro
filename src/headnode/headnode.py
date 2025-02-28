from flask import Flask, request, jsonify
import requests
import yaml
import json
import sqlite3
import threading
import uuid
import time
from functools import wraps
import hashlib
import os
from datetime import datetime

app = Flask(__name__)


# Load configuration from YAML file
def load_config():
    with open('config.yml', 'r') as file:
        return yaml.safe_load(file)


# Initialize database
def init_db():
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Create logs table for persistent logging
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS logs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        request_id TEXT,
        timestamp TEXT,
        smiles TEXT,
        parameters TEXT,
        status TEXT,
        worker_node TEXT,
        execution_time REAL,
        error_message TEXT
    )
    ''')

    # Modified cache table to store individual parameter sets and their results
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS parameter_cache (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        smiles TEXT,
        parameters TEXT,
        output_json_1 TEXT,
        output_json_2 TEXT,
        output_json_3 TEXT,
        timestamp_1 TEXT,
        timestamp_2 TEXT,
        timestamp_3 TEXT,
        UNIQUE(smiles, parameters)
    )
    ''')

    # Keep the original cache table for compatibility
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS cache (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        smiles TEXT,
        parameters TEXT,
        output_json_1 TEXT,
        output_json_2 TEXT,
        output_json_3 TEXT,
        output_json_4 TEXT,
        output_json_5 TEXT,
        output_json_6 TEXT,
        output_json_7 TEXT,
        output_json_8 TEXT,
        timestamp TEXT,
        UNIQUE(smiles, parameters)
    )
    ''')

    # Create a table to track running jobs
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS running_jobs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        request_id TEXT,
        smiles TEXT,
        parameters TEXT,
        start_time TEXT,
        worker_node TEXT,
        worker_index INTEGER,
        UNIQUE(request_id, worker_index)
    )
    ''')

    conn.commit()
    conn.close()


# Create a hash of parameters for consistent comparison
def hash_parameters(params):
    # Sort the keys to ensure consistent order
    ordered_params = json.dumps(params, sort_keys=True)
    return hashlib.md5(ordered_params.encode()).hexdigest()


# Authentication decorator
def require_api_key(f):

    @wraps(f)
    def decorated_function(*args, **kwargs):
        config = load_config()
        api_key = request.headers.get('X-API-KEY')
        if api_key and api_key == config['api_key']:
            return f(*args, **kwargs)
        else:
            # Log failed authentication attempt
            log_event(None, None, None, 'auth_failed', None, 0,
                      'Invalid API key')
            return jsonify({"error": "Unauthorized"}), 401

    return decorated_function


# Log events to the database
def log_event(request_id,
              smiles,
              parameters,
              status,
              worker_node,
              execution_time,
              error_message=None,
              worker_index=None):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()
    timestamp = datetime.now().isoformat()

    # Convert parameters to JSON string if it's a dictionary or list
    if isinstance(parameters, (dict, list)):
        parameters = json.dumps(parameters)

    worker_info = worker_node
    if worker_index is not None:
        worker_info = f"{worker_node}[{worker_index}]" if worker_node else None

    cursor.execute(
        '''
    INSERT INTO logs (request_id, timestamp, smiles, parameters, status, worker_node, execution_time, error_message)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ''', (request_id, timestamp, smiles, parameters, status, worker_info,
          execution_time, error_message))

    conn.commit()
    conn.close()


# Check if a job is already running for a specific worker index and parameter set
def is_job_running(request_id, smiles, parameters, worker_index):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Specific parameter set for this worker
    params_hash = hash_parameters(parameters)

    cursor.execute(
        '''
    SELECT request_id, worker_node, start_time FROM running_jobs
    WHERE request_id = ? AND smiles = ? AND parameters = ? AND worker_index = ?
    ''', (request_id, smiles, params_hash, worker_index))

    result = cursor.fetchone()

    # If not found with this request ID, check if there's any job running with these parameters
    if not result:
        cursor.execute(
            '''
        SELECT request_id, worker_node, start_time FROM running_jobs
        WHERE smiles = ? AND parameters = ?
        ''', (smiles, params_hash))
        result = cursor.fetchone()

    conn.close()

    return result


# Register a running job for a specific worker
def register_running_job(request_id, smiles, parameters, worker_node,
                         worker_index):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)
    start_time = datetime.now().isoformat()

    try:
        cursor.execute(
            '''
        INSERT INTO running_jobs (request_id, smiles, parameters, start_time, worker_node, worker_index)
        VALUES (?, ?, ?, ?, ?, ?)
        ''', (request_id, smiles, params_hash, start_time, worker_node,
              worker_index))
        conn.commit()
        success = True
    except sqlite3.IntegrityError:
        # Job is already running
        success = False

    conn.close()
    return success


# Remove a job from running_jobs when complete
def complete_running_job(request_id, worker_index):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    cursor.execute(
        '''
    DELETE FROM running_jobs WHERE request_id = ? AND worker_index = ?
    ''', (request_id, worker_index))

    conn.commit()
    conn.close()


# Check parameter cache for a specific parameter set
def check_parameter_cache(smiles, parameters):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)

    cursor.execute(
        '''
    SELECT output_json_1, output_json_2, output_json_3, 
           timestamp_1, timestamp_2, timestamp_3 
    FROM parameter_cache
    WHERE smiles = ? AND parameters = ?
    ''', (smiles, params_hash))

    result = cursor.fetchone()
    conn.close()

    if result:
        # Extract outputs and timestamps
        outputs = result[:3]
        timestamps = result[3:]

        # Create a list of results with their timestamps
        all_results = []
        for i, (output, timestamp) in enumerate(zip(outputs, timestamps)):
            if output:
                try:
                    # Parse JSON and add metadata
                    parsed_result = json.loads(output)
                    parsed_result["cache_index"] = i + 1
                    parsed_result["timestamp"] = timestamp
                    all_results.append(parsed_result)
                except:
                    pass

        # Sort by timestamp (newest first) if available
        all_results.sort(key=lambda x: x.get("timestamp", ""), reverse=True)

        # Return newest result for backward compatibility
        return all_results[0] if all_results else None, all_results
    return None, []


# Check cache for all results in a parameters list
def check_all_parameter_cache(smiles, parameters_list):
    results = []
    all_results = []
    cache_hits = True

    for params in parameters_list:
        result, all_param_results = check_parameter_cache(smiles, params)
        if result is not None:
            results.append(result)
            all_results.append(all_param_results)
        else:
            results.append(None)
            all_results.append([])
            cache_hits = False

    return results, all_results, cache_hits


# Update parameter cache with result, using a rotating system of 3 slots
def update_parameter_cache(smiles, parameters, result):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)
    current_timestamp = datetime.now().isoformat()

    # Convert result to JSON string
    result_json = json.dumps(result) if result is not None else None

    # Check if this parameter set exists
    cursor.execute(
        '''
    SELECT output_json_1, output_json_2, output_json_3, 
           timestamp_1, timestamp_2, timestamp_3 
    FROM parameter_cache
    WHERE smiles = ? AND parameters = ?
    ''', (smiles, params_hash))

    existing = cursor.fetchone()

    if existing:
        print(
            f"found existing cache entry for {smiles} with params {params_hash}"
        )
        # Determine which slot to update (the oldest one or the first empty one)
        outputs = existing[:3]  # The three output JSONs
        timestamps = existing[3:]  # The three timestamps

        # Find first empty slot
        empty_slot = None
        for i, output in enumerate(outputs):
            if output is None:
                empty_slot = i + 1  # Slots are 1-indexed in the database
                break

        # If no empty slot, find oldest
        if empty_slot is None:
            # Convert timestamps to datetime objects for comparison
            dt_timestamps = []
            for ts in timestamps:
                try:
                    dt_timestamps.append(
                        datetime.fromisoformat(ts) if ts else datetime.max)
                except ValueError:
                    dt_timestamps.append(datetime.max)

            # Find index of oldest timestamp
            oldest_idx = dt_timestamps.index(min(dt_timestamps))
            empty_slot = oldest_idx + 1  # Slots are 1-indexed

        # Update the selected slot
        cursor.execute(
            f'''
        UPDATE parameter_cache
        SET output_json_{empty_slot} = ?, timestamp_{empty_slot} = ?
        WHERE smiles = ? AND parameters = ?
        ''', (result_json, current_timestamp, smiles, params_hash))
    else:
        # Insert new record with the result in the first slot
        cursor.execute(
            '''
        INSERT INTO parameter_cache 
        (smiles, parameters, output_json_1, output_json_2, output_json_3, 
         timestamp_1, timestamp_2, timestamp_3)
        VALUES (?, ?, ?, NULL, NULL, ?, NULL, NULL)
        ''', (smiles, params_hash, result_json, current_timestamp))

    conn.commit()
    conn.close()


# Get available worker nodes
def get_available_worker_nodes():
    config = load_config()
    worker_nodes = config['worker_nodes']
    available_nodes = []

    for node in worker_nodes:
        try:
            health_url = f"http://{node['host']}:{node['port']}/api/health"
            headers = {'X-API-KEY': node['api_key']}
            response = requests.get(health_url, headers=headers, timeout=2)
            if response.status_code == 200:
                available_nodes.append(node)
        except requests.exceptions.RequestException:
            continue

    return available_nodes


# Forward request to selected worker node
def forward_to_worker(worker_node, endpoint, data):
    url = f"http://{worker_node['host']}:{worker_node['port']}/api/{endpoint}"
    headers = {
        'X-API-KEY': worker_node['api_key'],
        'Content-Type': 'application/json'
    }

    try:
        response = requests.post(url, headers=headers, json=data, timeout=60)
        return response.json(), response.status_code
    except requests.exceptions.RequestException as e:
        return {"error": f"Worker node communication error: {str(e)}"}, 500


# Process retrosynthesis request asynchronously for a specific worker with specific parameters
def process_retrosynthesis(request_id, smiles, parameters, worker_node,
                           worker_index):
    start_time = time.time()

    # Add the SMILES to the parameters if not already there
    if "smiles" not in parameters:
        parameters["smiles"] = smiles

    try:
        # Forward request to worker node
        result, status_code = forward_to_worker(worker_node, "retrosynthesis",
                                                parameters)

        # Calculate execution time
        execution_time = time.time() - start_time

        if status_code == 200:
            # Update cache with this worker's result
            update_parameter_cache(smiles, parameters, result)

            # Add parameter info to the result
            result["parameters_used"] = parameters

            log_event(request_id,
                      smiles,
                      parameters,
                      'completed',
                      worker_node['host'],
                      execution_time,
                      worker_index=worker_index)
        else:
            log_event(request_id,
                      smiles,
                      parameters,
                      'failed',
                      worker_node['host'],
                      execution_time,
                      error_message=result.get('error', 'Unknown error'),
                      worker_index=worker_index)
    except Exception as e:
        execution_time = time.time() - start_time
        log_event(request_id,
                  smiles,
                  parameters,
                  'failed',
                  worker_node['host'],
                  execution_time,
                  error_message=str(e),
                  worker_index=worker_index)
    finally:
        # Remove from running jobs
        complete_running_job(request_id, worker_index)


@app.route('/api/retrosynthesis', methods=['POST'])
@require_api_key
def retrosynthesis_api():
    """
    Headnode endpoint to handle retrosynthesis requests and distribute them to multiple worker nodes
    with different parameter sets.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({
            "error":
            "SMILES string is required. Please include a 'smiles' field"
        }), 400

    smiles = data['smiles']

    # Check for parameters_list, which should contain up to 8 parameter sets
    if 'parameters_list' not in data or not isinstance(data['parameters_list'],
                                                       list):
        return jsonify({
            "error":
            "parameters_list is required and must be a list of parameter sets"
        }), 400

    parameters_list = data['parameters_list']

    # Limit to maximum 8 parameter sets
    if len(parameters_list) > 8:
        return jsonify({
            "error":
            "Maximum of 8 parameter sets allowed in parameters_list"
        }), 400

    # Generate a unique request ID
    request_id = str(uuid.uuid4())

    # Check if all parameters are already in the cache
    cached_results, all_cached_results, all_cached = check_all_parameter_cache(
        smiles, parameters_list)

    # Get available worker nodes
    available_nodes = get_available_worker_nodes()
    if not available_nodes:
        log_event(request_id, smiles, parameters_list, 'no_workers', None, 0,
                  'No available worker nodes')
        return jsonify({"error": "No worker nodes available"}), 503

    # Track which workers are processing which parameter sets
    active_workers = []

    # Start a job for each parameter set, regardless of cache status
    for i, params in enumerate(parameters_list):
        # Check if this specific job is already running
        running_job = is_job_running(request_id, smiles, params, i)
        if running_job:
            job_id, worker, start_time = running_job
            active_workers.append({
                "worker_index": i,
                "worker_node": worker,
                "start_time": start_time,
                "request_id": job_id,
                "status": "already_processing",
                "parameter_set": params,
                "cached": cached_results[i] is not None
            })
            continue

        # Select the worker node for this parameter set
        worker_node = available_nodes[i % len(available_nodes)]

        # Register the job as running
        if register_running_job(request_id, smiles, params,
                                worker_node['host'], i):
            # Log the start of the job
            log_event(request_id,
                      smiles,
                      params,
                      'started',
                      worker_node['host'],
                      0,
                      worker_index=i)

            # Start processing in a separate thread
            threading.Thread(target=process_retrosynthesis,
                             args=(request_id, smiles, params, worker_node,
                                   i)).start()

            active_workers.append({
                "worker_index": i,
                "worker_node": worker_node['host'],
                "status": "processing",
                "parameter_set": params,
                "cached": cached_results[i] is not None
            })

    # Return cached results if available, while calculations continue in background
    if any(r is not None for r in cached_results):
        # Add parameter info to each result
        for i, result in enumerate(cached_results):
            if result:
                result["parameters_used"] = parameters_list[i]

        return jsonify({
            "status": "partial_results_available",
            "message":
            "Some results are available from cache, additional calculations in progress",
            "request_id": request_id,
            "results": cached_results,
            "all_cached_results":
            all_cached_results,  # Include all cached results
            "from_cache": True,
            "background_processing": True,
            "worker_status": active_workers
        }), 200
    else:
        return jsonify({
            "status": "processing",
            "message": "Retrosynthesis request distributed across workers",
            "request_id": request_id,
            "worker_status": active_workers,
        }), 202


@app.route('/api/result/<request_id>', methods=['GET'])
@require_api_key
def get_result(request_id):
    """
    Endpoint to check the status of a multi-worker request by its ID.
    """
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Get all jobs related to this request
    cursor.execute(
        '''
    SELECT smiles, parameters, status, worker_node FROM logs
    WHERE request_id = ? AND status IN ('started', 'completed', 'failed')
    ORDER BY timestamp DESC
    ''', (request_id, ))

    all_jobs = cursor.fetchall()

    if not all_jobs:
        conn.close()
        return jsonify({"error": "Request ID not found"}), 404

    # Extract smiles and all parameter sets used
    smiles_set = set()
    parameters_by_index = {}
    status_by_index = {}

    for job in all_jobs:
        smiles, params_str, status, worker_index = job
        smiles_set.add(smiles)

        if worker_index is not None and worker_index not in parameters_by_index:
            try:
                params = json.loads(params_str)
                if isinstance(params,
                              dict):  # This is an individual parameter set
                    parameters_by_index[worker_index] = params
                    status_by_index[worker_index] = status
            except:
                pass

    if not smiles_set:
        conn.close()
        return jsonify({"error":
                        "No valid data found for this request ID"}), 404

    smiles = list(smiles_set)[0]  # Use the first SMILES string found

    # Sort parameter sets by index
    parameters_list = []
    max_index = max(parameters_by_index.keys()) if parameters_by_index else -1
    max_index = int(max_index.split("[")[1].split(']')[0]) if isinstance(
        max_index, str) else max_index

    for i in range(max_index + 1):
        if i in parameters_by_index:
            parameters_list.append(parameters_by_index[i])

    # Check for running jobs
    cursor.execute(
        '''
    SELECT worker_index, worker_node, parameters, start_time FROM running_jobs
    WHERE request_id = ?
    ''', (request_id, ))

    running_jobs = cursor.fetchall()
    conn.close()

    # If we have running jobs, build status information
    if running_jobs:
        running_workers = []
        for job in running_jobs:
            index, worker, params_hash, start_time = job
            try:
                params = parameters_by_index.get(index, {"unknown": True})
            except:
                params = {"unknown": True}

            running_workers.append({
                "worker_index": index,
                "worker_node": worker,
                "start_time": start_time,
                "parameter_set": params
            })

        # Get partial results from cache
        results = []
        for i, params in enumerate(parameters_list):
            if params:
                cached_result = check_parameter_cache(smiles, params)
                if cached_result:
                    cached_result["parameters_used"] = params
                results.append(cached_result)
            else:
                results.append(None)

        # Determine if we have some results but not all
        result_count = sum(1 for r in results if r is not None)
        expected_count = len(parameters_list)

        return jsonify({
            "status": "processing",
            "message": f"Job is still running on {len(running_jobs)} workers",
            "request_id": request_id,
            "running_workers": running_workers,
            "partial_results": results,
            "progress": f"{result_count}/{expected_count} completed"
        }), 202

    # No running jobs - get all cached results
    results = []
    failed_jobs = []

    for i, params in enumerate(parameters_list):
        if params:
            cached_result = check_parameter_cache(smiles, params)

            if cached_result:
                cached_result["parameters_used"] = params
                results.append(cached_result)
            else:
                results.append(None)
                # This job might have failed
                if status_by_index.get(i) == 'failed':
                    failed_jobs.append({
                        "worker_index": i,
                        "parameter_set": params,
                        "status": "failed"
                    })

    # If we have at least one result, return partial success
    if any(r is not None for r in results):
        return jsonify({
            "status":
            "completed",
            "request_id":
            request_id,
            "results":
            results,
            "failed_jobs":
            failed_jobs if failed_jobs else None,
            "success_count":
            sum(1 for r in results if r is not None),
            "total_count":
            len(parameters_list)
        }), 200

    # If all jobs failed, return error status
    return jsonify({
        "status": "failed",
        "message": "All jobs failed or no results found",
        "request_id": request_id,
        "failed_jobs": failed_jobs
    }), 500


@app.route('/api/parameter/retrosynthesis', methods=['POST'])
@require_api_key
def single_parameter_retrosynthesis():
    """
    Endpoint to perform retrosynthesis with a single parameter set.
    This is useful for direct parameter testing or individual refinement.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({
            "error":
            "SMILES string is required. Please include a 'smiles' field"
        }), 400

    smiles = data['smiles']

    # Use provided parameters or default if none given
    parameters = data.get('parameters', {"smiles": smiles})
    if "smiles" not in parameters:
        parameters["smiles"] = smiles

    # Check if this exact parameter set is already in cache
    cached_result = check_parameter_cache(smiles, parameters)

    # Generate a unique request ID
    request_id = str(uuid.uuid4())

    # Select a worker node
    available_nodes = get_available_worker_nodes()
    if not available_nodes:
        log_event(request_id, smiles, parameters, 'no_workers', None, 0,
                  'No available worker nodes')
        return jsonify({"error": "No worker nodes available"}), 503

    worker_node = available_nodes[0]  # Use the first available node

    # Check if job is already running
    running_job = is_job_running(request_id, smiles, parameters, 0)
    if running_job:
        job_id, worker, start_time = running_job
        # Return cached result if available while job is running
        if cached_result:
            cached_result["parameters_used"] = parameters
            return jsonify({
                "status": "partial_results_available",
                "message":
                f"This request is already being processed by worker {worker}, partial results available",
                "request_id": job_id,
                "result": cached_result,
                "from_cache": True
            }), 200
        else:
            return jsonify({
                "status": "processing",
                "message":
                f"This request is already being processed by worker {worker}",
                "request_id": job_id
            }), 202

    # Register the job
    if register_running_job(request_id, smiles, parameters,
                            worker_node['host'], 0):
        # Log the start
        log_event(request_id, smiles, parameters, 'started',
                  worker_node['host'], 0)

        # Start processing
        threading.Thread(target=process_retrosynthesis,
                         args=(request_id, smiles, parameters, worker_node,
                               0)).start()

        # Return cached result if available while job starts
        if cached_result:
            cached_result["parameters_used"] = parameters
            return jsonify({
                "status": "partial_results_available",
                "message":
                "Additional calculation in progress, partial results available",
                "request_id": request_id,
                "worker_node": worker_node['host'],
                "result": cached_result,
                "from_cache": True
            }), 200
        else:
            return jsonify({
                "status": "processing",
                "message": "Retrosynthesis request has been submitted",
                "request_id": request_id,
                "worker_node": worker_node['host']
            }), 202
    else:
        return jsonify({"error": "Failed to register job"}), 500


@app.route('/api/check_cache', methods=['POST'])
@require_api_key
def check_cache_endpoint():
    """
    Endpoint to check if specific parameters are cached without running the calculation.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({"error": "SMILES string is required"}), 400

    smiles = data['smiles']

    if 'parameters_list' in data and isinstance(data['parameters_list'], list):
        # Check multiple parameter sets
        results, all_cached = check_all_parameter_cache(
            smiles, data['parameters_list'])

        # Add parameter info to results
        for i, result in enumerate(results):
            if result:
                result["parameters_used"] = data['parameters_list'][i]

        return jsonify({
            "all_cached":
            all_cached,
            "results":
            results,
            "parameters_cached":
            sum(1 for r in results if r is not None),
            "parameters_total":
            len(data['parameters_list'])
        }), 200

    elif 'parameters' in data:
        # Check single parameter set
        result = check_parameter_cache(smiles, data['parameters'])
        if result:
            result["parameters_used"] = data['parameters']
            return jsonify({"cached": True, "result": result}), 200
        else:
            return jsonify({"cached": False}), 200

    else:
        return jsonify(
            {"error":
             "Either parameters or parameters_list must be provided"}), 400


@app.route('/api/health', methods=['GET'])
@require_api_key
def health():
    """
    Endpoint to check the health of the headnode API.
    """
    # Check if worker nodes are available
    available_nodes = get_available_worker_nodes()
    config = load_config()
    worker_nodes = config['worker_nodes']

    return jsonify({
        "status": "healthy",
        "worker_nodes": {
            "total": len(worker_nodes),
            "available": len(available_nodes)
        }
    }), 200


@app.route('/api/clear_cache', methods=['POST'])
@require_api_key
def clear_specific_cache():
    """
    Endpoint to clear the cache for a specific molecule and parameter set.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({"error": "Molecule SMILES is required"}), 400

    smiles = data['smiles']
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    if 'parameters' in data:
        # Clear a specific parameter set
        params_hash = hash_parameters(data['parameters'])
        cursor.execute(
            'DELETE FROM parameter_cache WHERE smiles = ? AND parameters = ?',
            (smiles, params_hash))
    elif 'parameters_list' in data and isinstance(data['parameters_list'],
                                                  list):
        # Clear multiple parameter sets
        for params in data['parameters_list']:
            params_hash = hash_parameters(params)
            cursor.execute(
                'DELETE FROM parameter_cache WHERE smiles = ? AND parameters = ?',
                (smiles, params_hash))
    else:
        # Clear all cache entries for this molecule
        cursor.execute('DELETE FROM parameter_cache WHERE smiles = ?',
                       (smiles, ))

    conn.commit()
    conn.close()

    return jsonify({
        "status":
        "success",
        "message":
        "Cache cleared for specified molecule and parameters"
    }), 200


@app.route('/api/stats', methods=['GET'])
@require_api_key
def get_stats():
    """
    Endpoint to get statistics about the system.
    """
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Get parameter cache stats
    cursor.execute('SELECT COUNT(*) FROM parameter_cache')
    param_cache_count = cursor.fetchone()[0]

    # Get running jobs count
    cursor.execute('SELECT COUNT(*) FROM running_jobs')
    running_count = cursor.fetchone()[0]

    # Get recent jobs
    cursor.execute('''
    SELECT request_id, timestamp, smiles, status, worker_node, execution_time
    FROM logs
    ORDER BY timestamp DESC
    LIMIT 10
    ''')

    recent_jobs = []
    for row in cursor.fetchall():
        recent_jobs.append({
            'request_id': row[0],
            'timestamp': row[1],
            'smiles': row[2],
            'status': row[3],
            'worker_node': row[4],
            'execution_time': row[5]
        })

    # Get most frequently cached molecules
    cursor.execute('''
    SELECT smiles, COUNT(*) as count
    FROM parameter_cache
    GROUP BY smiles
    ORDER BY count DESC
    LIMIT 5
    ''')

    frequent_molecules = []
    for row in cursor.fetchall():
        frequent_molecules.append({'smiles': row[0], 'cache_entries': row[1]})

    conn.close()

    return jsonify({
        "parameter_cache_entries": param_cache_count,
        "running_jobs": running_count,
        "recent_jobs": recent_jobs,
        "frequent_molecules": frequent_molecules
    }), 200


if __name__ == '__main__':
    # Create the database if it doesn't exist
    if not os.path.exists('retrosynthesis.db'):
        init_db()

    # Load configuration
    config = load_config()

    # Start the Flask app
    app.run(host=config.get('host', '0.0.0.0'),
            port=config.get('port', 5001),
            debug=config.get('debug', False))
