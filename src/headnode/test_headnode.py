import requests
import json
import time
import argparse


def submit_retrosynthesis(url, api_key, smiles, parameters_list=None):
    """
    Submit a retrosynthesis job to the headnode server with multiple parameter sets.
    """
    headers = {'X-API-KEY': api_key, 'Content-Type': 'application/json'}

    data = {'smiles': smiles}

    # If parameters_list is provided, use it; otherwise, create a default one
    if parameters_list:
        data['parameters_list'] = parameters_list
    else:
        # Create a default parameter set
        default_params = {'smiles': smiles}
        data['parameters_list'] = [default_params]

    print(
        f"Submitting request for SMILES: {smiles} with {len(data['parameters_list'])} parameter sets"
    )
    response = requests.post(f"{url}/api/retrosynthesis",
                             headers=headers,
                             json=data)

    return response.json(), response.status_code


def submit_single_parameter_retrosynthesis(url,
                                           api_key,
                                           smiles,
                                           parameters=None):
    """
    Submit a single parameter retrosynthesis job to the headnode server.
    """
    headers = {'X-API-KEY': api_key, 'Content-Type': 'application/json'}

    data = {'smiles': smiles}

    # Add parameters if provided
    if parameters:
        data['parameters'] = parameters

    print(f"Submitting single parameter request for SMILES: {smiles}")
    response = requests.post(f"{url}/api/parameter/retrosynthesis",
                             headers=headers,
                             json=data)

    return response.json(), response.status_code


def check_job_status(url, api_key, request_id):
    """
    Check the status of a job.
    """
    headers = {'X-API-KEY': api_key}

    response = requests.get(f"{url}/api/result/{request_id}", headers=headers)

    return response.json(), response.status_code


def check_cache(url, api_key, smiles, parameters_list=None):
    """
    Check if results for specific parameters are already in cache.
    """
    headers = {'X-API-KEY': api_key, 'Content-Type': 'application/json'}

    data = {'smiles': smiles}
    if parameters_list:
        data['parameters_list'] = parameters_list

    response = requests.post(f"{url}/api/check_cache",
                             headers=headers,
                             json=data)

    return response.json(), response.status_code


def wait_for_completion(url, api_key, request_id, max_wait=300, interval=5):
    """
    Wait for a job to complete, polling at regular intervals.
    """
    elapsed = 0
    while elapsed < max_wait:
        result, status_code = check_job_status(url, api_key, request_id)

        if status_code == 200:
            print("Job completed successfully")
            print(
                f"Success count: {result.get('success_count', 0)}/{result.get('total_count', 0)}"
            )
            return result
        elif status_code == 500:
            print(f"Job failed: {result.get('message', 'Unknown error')}")
            return None
        elif status_code == 202:
            progress = result.get('progress', 'Unknown')
            print(
                f"Job still processing: {result.get('message', 'In progress')} - {progress}"
            )

            # Print partial results if available
            if 'partial_results' in result and result['partial_results']:
                print(
                    f"Partial results available for {sum(1 for r in result['partial_results'] if r is not None)} parameter sets"
                )
        else:
            print(f"Unexpected status code: {status_code}")
            return None

        time.sleep(interval)
        elapsed += interval

    print(f"Timed out after waiting {max_wait} seconds")
    return None


def get_server_stats(url, api_key):
    """
    Get statistics about the server.
    """
    headers = {'X-API-KEY': api_key}

    response = requests.get(f"{url}/api/stats", headers=headers)

    return response.json(), response.status_code


def check_server_health(url, api_key):
    """
    Check if the server is healthy.
    """
    headers = {'X-API-KEY': api_key}

    response = requests.get(f"{url}/api/health", headers=headers)

    return response.json(), response.status_code


def main():
    parser = argparse.ArgumentParser(
        description='Test client for retrosynthesis headnode')
    parser.add_argument('--url',
                        default='http://localhost:5001',
                        help='Headnode server URL')
    parser.add_argument('--api-key',
                        default='your-secure-api-key',
                        help='API key')
    parser.add_argument('--smiles', required=True, help='SMILES string')

    # Command options
    parser.add_argument(
        '--command',
        choices=['run', 'check-cache', 'stats', 'health', 'single'],
        default='run',
        help='Command to execute')

    # Parameter options
    parser.add_argument('--advanced-model',
                        action='store_true',
                        help='Use advanced model')
    parser.add_argument('--advanced-prompt',
                        action='store_true',
                        help='Use advanced prompt')
    parser.add_argument('--model-version',
                        default='USPTO',
                        help='AiZynthFinder model version')
    parser.add_argument('--stability-flag',
                        action='store_true',
                        help='Enable stability check')
    parser.add_argument('--hallucination-check',
                        action='store_true',
                        help='Enable hallucination check')

    # Multi-parameter run options
    parser.add_argument('--run-variants',
                        action='store_true',
                        help='Run multiple parameter variants')
    parser.add_argument('--max-wait',
                        type=int,
                        default=300,
                        help='Maximum wait time in seconds')
    parser.add_argument('--poll-interval',
                        type=int,
                        default=5,
                        help='Polling interval in seconds')

    args = parser.parse_args()

    # Basic parameter set
    base_params = {
        'advanced_model': str(args.advanced_model).lower(),
        'advanced_prompt': str(args.advanced_prompt).lower(),
        'model_version': args.model_version,
        'stability_flag': str(args.stability_flag).lower(),
        'hallucination_check': str(args.hallucination_check).lower(),
        'smiles': args.smiles
    }

    # Handle different commands
    if args.command == 'check-cache':
        if args.run_variants:
            # Generate parameter variants for cache check
            parameters_list = [
                base_params, {
                    **base_params, 'advanced_model': 'true'
                }, {
                    **base_params, 'advanced_prompt': 'true'
                }, {
                    **base_params, 'stability_flag': 'true'
                }
            ]
            result, status_code = check_cache(args.url, args.api_key,
                                              args.smiles, parameters_list)
        else:
            # Check single parameter set
            result, status_code = check_cache(args.url, args.api_key,
                                              args.smiles, [base_params])

        if status_code == 200:
            print("Cache check results:")
            print(f"All cached: {result.get('all_cached', False)}")
            print(
                f"Parameters cached: {result.get('parameters_cached', 0)}/{result.get('parameters_total', 0)}"
            )
            print(json.dumps(result, indent=2))
        else:
            print(
                f"Error checking cache: {result.get('error', 'Unknown error')}"
            )

    elif args.command == 'stats':
        result, status_code = get_server_stats(args.url, args.api_key)
        if status_code == 200:
            print("Server statistics:")
            print(json.dumps(result, indent=2))
        else:
            print(
                f"Error getting stats: {result.get('error', 'Unknown error')}")

    elif args.command == 'health':
        result, status_code = check_server_health(args.url, args.api_key)
        if status_code == 200:
            print("Server health:")
            print(json.dumps(result, indent=2))
        else:
            print(
                f"Error checking health: {result.get('error', 'Unknown error')}"
            )

    elif args.command == 'single':
        # Submit a single parameter job
        result, status_code = submit_single_parameter_retrosynthesis(
            args.url, args.api_key, args.smiles, base_params)

        if status_code == 200:
            # Job completed immediately (cache hit)
            print("Result retrieved from cache:")
            print(json.dumps(result, indent=2))
        elif status_code == 202:
            # Job is processing
            request_id = result.get('request_id')
            print(f"Job submitted with request ID: {request_id}")
            print(f"Status: {result.get('message')}")

            # Wait for completion
            final_result = wait_for_completion(args.url, args.api_key,
                                               request_id, args.max_wait,
                                               args.poll_interval)
            if final_result:
                print("Final result:")
                print(json.dumps(final_result, indent=2))
        else:
            # Error
            print(f"Error: {result.get('error', 'Unknown error')}")

    else:  # Default 'run' command
        if args.run_variants:
            # Generate parameter variants
            parameters_list = [
                base_params, {
                    **base_params, 'advanced_model': 'true'
                }, {
                    **base_params, 'advanced_prompt': 'true'
                }, {
                    **base_params, 'stability_flag': 'true'
                }, {
                    **base_params, 'hallucination_check': 'true'
                }, {
                    **base_params, 'advanced_model': 'true',
                    'advanced_prompt': 'true'
                }, {
                    **base_params, 'stability_flag': 'true',
                    'hallucination_check': 'true'
                }, {
                    **base_params, 'advanced_model': 'true',
                    'advanced_prompt': 'true',
                    'stability_flag': 'true'
                }
            ]
        else:
            # Use single parameter set
            parameters_list = [base_params]

        # Submit the job
        result, status_code = submit_retrosynthesis(args.url, args.api_key,
                                                    args.smiles,
                                                    parameters_list)

        print("--------------------")
        print(f"Response status code: {status_code}")
        print("--------------------")
        print(f"Response: {result}")
        print("--------------------")

        if status_code == 200:
            # Job completed immediately (cache hit)
            print("Results retrieved from cache:")
            print(json.dumps(result, indent=2))
        elif status_code == 202:
            # Job is processing
            request_id = result.get('request_id')
            print(f"Job submitted with request ID: {request_id}")
            print(f"Status: {result.get('message')}")

            # Print worker status
            if 'worker_status' in result:
                print("\nWorker Status:")
                for worker in result['worker_status']:
                    print(
                        f"  Worker {worker.get('worker_index')}: {worker.get('status')} on {worker.get('worker_node', 'N/A')}"
                    )

            # Print partial results if available
            if 'partial_results' in result and result[
                    'partial_results'] and any(result['partial_results']):
                print("\nPartial results available from cache")

            # Wait for completion
            final_result = wait_for_completion(args.url, args.api_key,
                                               request_id, args.max_wait,
                                               args.poll_interval)
            if final_result:
                print("Final result:")
                print(json.dumps(final_result, indent=2))
        else:
            # Error
            print(f"Error: {result.get('error', 'Unknown error')}")


if __name__ == "__main__":
    main()
