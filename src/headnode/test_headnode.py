import requests
import json
import time
import argparse


def submit_retrosynthesis(url, api_key, smiles, options=None):
    """
    Submit a retrosynthesis job to the headnode server.
    """
    headers = {'X-API-KEY': api_key, 'Content-Type': 'application/json'}

    data = {'smiles': smiles}

    # Add any additional options
    if options:
        data.update(options)

    print(f"Submitting request for SMILES: {smiles}")
    response = requests.post(f"{url}/api/retrosynthesis",
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


def wait_for_completion(url, api_key, request_id, max_wait=300, interval=5):
    """
    Wait for a job to complete, polling at regular intervals.
    """
    elapsed = 0
    while elapsed < max_wait:
        result, status_code = check_job_status(url, api_key, request_id)

        if status_code == 200:
            print("Job completed successfully")
            return result
        elif status_code == 500:
            print(f"Job failed: {result.get('message', 'Unknown error')}")
            return None
        elif status_code == 202:
            print(
                f"Job still processing: {result.get('message', 'In progress')}"
            )
        else:
            print(f"Unexpected status code: {status_code}")
            return None

        time.sleep(interval)
        elapsed += interval

    print(f"Timed out after waiting {max_wait} seconds")
    return None


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

    args = parser.parse_args()

    # Prepare additional options
    options = {
        'advanced_model': str(args.advanced_model).lower(),
        'advanced_prompt': str(args.advanced_prompt).lower(),
        'model_version': args.model_version,
        'stability_flag': str(args.stability_flag).lower(),
        'hallucination_check': str(args.hallucination_check).lower()
    }

    # Submit the job
    result, status_code = submit_retrosynthesis(args.url, args.api_key,
                                                args.smiles, options)

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
        final_result = wait_for_completion(args.url, args.api_key, request_id)
        if final_result:
            print("Final result:")
            print(json.dumps(final_result, indent=2))
    else:
        # Error
        print(f"Error: {result.get('error', 'Unknown error')}")


if __name__ == "__main__":
    main()
