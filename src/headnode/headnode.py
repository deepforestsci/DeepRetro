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

    # Create cache table
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
        request_id TEXT UNIQUE,
        smiles TEXT,
        parameters TEXT,
        start_time TEXT,
        worker_node TEXT,
        UNIQUE(smiles, parameters)
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
              error_message=None):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()
    timestamp = datetime.now().isoformat()

    # Convert parameters to JSON string if it's a dictionary
    if isinstance(parameters, dict):
        parameters = json.dumps(parameters)

    cursor.execute(
        '''
    INSERT INTO logs (request_id, timestamp, smiles, parameters, status, worker_node, execution_time, error_message)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ''', (request_id, timestamp, smiles, parameters, status, worker_node,
          execution_time, error_message))

    conn.commit()
    conn.close()


# Check if a job is already running
def is_job_running(smiles, parameters):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)

    cursor.execute(
        '''
    SELECT request_id, worker_node, start_time FROM running_jobs
    WHERE smiles = ? AND parameters = ?
    ''', (smiles, params_hash))

    result = cursor.fetchone()
    conn.close()

    return result


# Register a running job
def register_running_job(request_id, smiles, parameters, worker_node):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)
    start_time = datetime.now().isoformat()

    try:
        cursor.execute(
            '''
        INSERT INTO running_jobs (request_id, smiles, parameters, start_time, worker_node)
        VALUES (?, ?, ?, ?, ?)
        ''', (request_id, smiles, params_hash, start_time, worker_node))
        conn.commit()
        success = True
    except sqlite3.IntegrityError:
        # Job is already running
        success = False

    conn.close()
    return success


# Remove a job from running_jobs when complete
def complete_running_job(request_id):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    cursor.execute(
        '''
    DELETE FROM running_jobs WHERE request_id = ?
    ''', (request_id, ))

    conn.commit()
    conn.close()


# Check cache for existing results
def check_cache(smiles, parameters):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)

    cursor.execute(
        '''
    SELECT output_json_1, output_json_2, output_json_3, output_json_4, output_json_5, output_json_6, output_json_7, output_json_8
    FROM cache
    WHERE smiles = ? AND parameters = ?
    ''', (smiles, params_hash))

    result = cursor.fetchone()
    conn.close()

    if result:
        outputs = []
        for output_json in result:
            if output_json:
                outputs.append(json.loads(output_json))
            else:
                outputs.append(None)
        return outputs
    return None


# Update cache with new results
def update_cache(smiles, parameters, results):
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    params_hash = hash_parameters(parameters)
    timestamp = datetime.now().isoformat()

    # Prepare results for storage
    output_jsons = []
    for i in range(8):
        if i < len(results) and results[i] is not None:
            output_jsons.append(json.dumps(results[i]))
        else:
            output_jsons.append(None)

    cursor.execute(
        '''
    INSERT OR REPLACE INTO cache 
    (smiles, parameters, output_json_1, output_json_2, output_json_3, output_json_4, output_json_5, output_json_6, output_json_7, output_json_8, timestamp)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (smiles, params_hash, *output_jsons, timestamp))

    conn.commit()
    conn.close()


# Select a worker node using round-robin or load-based selection
def select_worker_node():
    config = load_config()
    worker_nodes = config['worker_nodes']

    # Simple round-robin for now, but could be enhanced with load balancing
    # Check each node for health and availability
    for node in worker_nodes:
        print(f"Checking node: {node['host']}:{node['port']}")
        health_url = f"http://{node['host']}:{node['port']}/api/health"
        try:
            health_url = f"http://{node['host']}:{node['port']}/api/health"
            headers = {'X-API-KEY': node['api_key']}
            response = requests.request("GET",
                                        health_url,
                                        headers=headers,
                                        data={},
                                        timeout=100)
            print(f"Response: {response}")
            if response.status_code == 200:
                return node
        except requests.exceptions.RequestException:
            continue

    # If no nodes are available
    return None


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


# Process retrosynthesis request asynchronously
def process_retrosynthesis(request_id, smiles, parameters, worker_node):
    start_time = time.time()

    try:
        # Forward request to worker node
        result, status_code = forward_to_worker(worker_node, "retrosynthesis",
                                                parameters)

        # Calculate execution time
        execution_time = time.time() - start_time

        if status_code == 200:
            # Update cache with results
            update_cache(smiles, parameters, [result])
            log_event(request_id, smiles, parameters, 'completed',
                      worker_node['host'], execution_time)
        else:
            log_event(request_id,
                      smiles,
                      parameters,
                      'failed',
                      worker_node['host'],
                      execution_time,
                      error_message=result.get('error', 'Unknown error'))
    except Exception as e:
        execution_time = time.time() - start_time
        log_event(request_id,
                  smiles,
                  parameters,
                  'failed',
                  worker_node['host'],
                  execution_time,
                  error_message=str(e))
    finally:
        # Remove from running jobs
        complete_running_job(request_id)


@app.route('/api/retrosynthesis', methods=['POST'])
@require_api_key
def retrosynthesis_api():
    """
    Headnode endpoint to handle retrosynthesis requests and distribute them to worker nodes.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({
            "error":
            "SMILES string is required. Please include a 'smiles' field"
        }), 400

    smiles = data['smiles']

    # Generate a unique request ID
    request_id = str(uuid.uuid4())

    # Check if this is already in the cache
    cached_results = check_cache(smiles, data)
    if cached_results and cached_results[0]:
        log_event(request_id, smiles, data, 'cache_hit', None, 0)
        return jsonify(cached_results[0]), 200

    # Check if job is already running
    running_job = is_job_running(smiles, data)
    if running_job:
        job_id, worker, start_time = running_job
        return jsonify({
            "status": "processing",
            "message":
            f"This molecule with the provided parameters is already being processed by worker {worker} since {start_time}",
            "request_id": job_id
        }), 202

    # Select a worker node
    worker_node = select_worker_node()
    if not worker_node:
        log_event(request_id, smiles, data, 'no_workers', None, 0,
                  'No available worker nodes')
        return jsonify({"error": "No worker nodes available"}), 503

    # Register the job as running
    if not register_running_job(request_id, smiles, data, worker_node['host']):
        # Another request just started processing this job
        running_job = is_job_running(smiles, data)
        if running_job:
            job_id, worker, start_time = running_job
            return jsonify({
                "status": "processing",
                "message":
                f"This molecule with the provided parameters is already being processed by worker {worker} since {start_time}",
                "request_id": job_id
            }), 202

    # Log the start of the job
    log_event(request_id, smiles, data, 'started', worker_node['host'], 0)

    # Start processing in a separate thread
    threading.Thread(target=process_retrosynthesis,
                     args=(request_id, smiles, data, worker_node)).start()

    return jsonify({
        "status": "processing",
        "message": "Retrosynthesis request has been queued",
        "request_id": request_id,
        "worker_node": worker_node['host']
    }), 202


@app.route('/api/result/<request_id>', methods=['GET'])
@require_api_key
def get_result(request_id):
    """
    Endpoint to check the status of a request by its ID.
    """
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Check if job is still running
    cursor.execute(
        '''
    SELECT smiles, parameters, worker_node, start_time FROM running_jobs
    WHERE request_id = ?
    ''', (request_id, ))

    running_job = cursor.fetchone()

    if running_job:
        smiles, params_hash, worker, start_time = running_job
        return jsonify({
            "status": "processing",
            "message":
            f"Job is still running on worker {worker} since {start_time}",
            "request_id": request_id
        }), 202

    # Check logs for the job
    cursor.execute(
        '''
    SELECT smiles, parameters, status, worker_node, timestamp, error_message
    FROM logs
    WHERE request_id = ? AND status IN ('completed', 'failed')
    ORDER BY timestamp DESC
    LIMIT 1
    ''', (request_id, ))

    job_log = cursor.fetchone()

    if not job_log:
        conn.close()
        return jsonify({"error": "Request ID not found"}), 404

    smiles, params_str, status, worker, timestamp, error = job_log

    if status == 'failed':
        conn.close()
        return jsonify({
            "status": "failed",
            "message": f"Job failed: {error}",
            "request_id": request_id
        }), 500

    # If completed, fetch from cache
    try:
        parameters = json.loads(params_str)
    except:
        parameters = {"smiles": smiles}

    cached_results = check_cache(smiles, parameters)

    conn.close()

    if cached_results and cached_results[0]:
        return jsonify(cached_results[0]), 200
    else:
        return jsonify({
            "status": "completed",
            "message": "Job completed but results not found in cache",
            "request_id": request_id
        }), 404


@app.route('/api/health', methods=['GET'])
@require_api_key
def health():
    """
    Endpoint to check the health of the headnode API.
    """
    # Check if worker nodes are available
    config = load_config()
    worker_nodes = config['worker_nodes']

    available_workers = 0
    for node in worker_nodes:
        try:
            health_url = f"http://{node['host']}:{node['port']}/api/health"
            headers = {'X-API-KEY': node['api_key']}
            response = requests.get(health_url, headers=headers, timeout=2)
            if response.status_code == 200:
                available_workers += 1
        except:
            continue

    return jsonify({
        "status": "healthy",
        "worker_nodes": {
            "total": len(worker_nodes),
            "available": available_workers
        }
    }), 200


@app.route('/api/clear_cache', methods=['POST'])
@require_api_key
def clear_specific_cache():
    """
    Endpoint to clear the cache for a specific molecule.
    """
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({"error": "Molecule SMILES is required"}), 400

    smiles = data['smiles']

    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    if 'parameters' in data:
        # Clear specific parameters for this molecule
        params_hash = hash_parameters(data['parameters'])
        cursor.execute('DELETE FROM cache WHERE smiles = ? AND parameters = ?',
                       (smiles, params_hash))
    else:
        # Clear all cache entries for this molecule
        cursor.execute('DELETE FROM cache WHERE smiles = ?', (smiles, ))

    conn.commit()
    conn.close()

    return jsonify({
        "status": "success",
        "message": "Cache cleared for specified molecule"
    }), 200


@app.route('/api/stats', methods=['GET'])
@require_api_key
def get_stats():
    """
    Endpoint to get statistics about the system.
    """
    conn = sqlite3.connect('retrosynthesis.db')
    cursor = conn.cursor()

    # Get cache stats
    cursor.execute('SELECT COUNT(*) FROM cache')
    cache_count = cursor.fetchone()[0]

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

    conn.close()

    return jsonify({
        "cache_entries": cache_count,
        "running_jobs": running_count,
        "recent_jobs": recent_jobs
    }), 200


# Initialize database on startup
# @app.before_first_request
# def before_first_request():
#     if not os.path.exists('retrosynthesis.db'):
#         init_db()

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
