import os
import boto3
from dotenv import load_dotenv
from botocore.exceptions import NoCredentialsError

load_dotenv()

# Configuration
S3_BUCKET = 'sicic-logs-dump'
FOLDER_PATH = '/home/ubuntu/recursiveLLM/logs/'

def upload_files_to_s3():
    s3_client = boto3.client(
    's3',
    aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
    aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
    # aws_session_token=SESSION_TOKEN
)
    
    for root, dirs, files in os.walk(FOLDER_PATH):
        for file in files:
            local_file_path = os.path.join(root, file)
            s3_key = os.path.relpath(local_file_path, FOLDER_PATH)
            
            if ".metaflow" not in local_file_path:
                try:
                    s3_client.upload_file(local_file_path, S3_BUCKET, s3_key)
                    print(f"Uploaded {local_file_path} to s3://{S3_BUCKET}/{s3_key}")
                except FileNotFoundError:
                    print(f"The file {local_file_path} was not found.")
                except NoCredentialsError:
                    print("Credentials not available.")

if __name__ == "__main__":
    upload_files_to_s3()