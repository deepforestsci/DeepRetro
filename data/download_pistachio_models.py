import os
import boto3
from dotenv import load_dotenv
from botocore.exceptions import NoCredentialsError

load_dotenv()

# Configuration
S3_BUCKET = 'pistachio-training-dump'
FOLDER_PATH = "/home/ubuntu/recursiveLLM/aizynthfinder/models/"

download_list = [
    "50_percent/uspto_expansion.onnx",
    "50_percent/uspto_ringbreaker_expansion.onnx",
    "50_percent/uspto_unique_templates.csv.gz",
    "50_percent/uspto_ringbreaker_unique_templates.csv.gz",
    "25_percent_backup/uspto_expansion.onnx",
    "25_percent_backup/uspto_ringbreaker_expansion.onnx",
    "25_percent_backup/uspto_unique_templates.csv.gz",
    "25_percent_backup/uspto_ringbreaker_unique_templates.csv.gz",
    "full/model_100e/uspto_expansion.onnx",
    "full/uspto_ringbreaker_expansion.onnx",
    "full/uspto_unique_templates.csv.gz",
    "full/uspto_ringbreaker_unique_templates.csv.gz",
    "full/model_150e/uspto_expansion.onnx",
    "full/uspto_ringbreaker_expansion.onnx",
    "full/uspto_unique_templates.csv.gz",
    "full/uspto_ringbreaker_unique_templates.csv.gz",
]

download_locations = [
    "Pistachio_50/uspto_expansion.onnx",
    "Pistachio_50/uspto_ringbreaker_expansion.onnx",
    "Pistachio_50/uspto_unique_templates.csv.gz",
    "Pistachio_50/uspto_ringbreaker_unique_templates.csv.gz",
    "Pistachio_25/uspto_expansion.onnx",
    "Pistachio_25/uspto_ringbreaker_expansion.onnx",
    "Pistachio_25/uspto_unique_templates.csv.gz",
    "Pistachio_25/uspto_ringbreaker_unique_templates.csv.gz",
    "Pistachio_100/uspto_expansion.onnx",
    "Pistachio_100/uspto_ringbreaker_expansion.onnx",
    "Pistachio_100/uspto_unique_templates.csv.gz",
    "Pistachio_100/uspto_ringbreaker_unique_templates.csv.gz",
    "Pistachio_100+/uspto_expansion.onnx",
    "Pistachio_100+/uspto_ringbreaker_expansion.onnx",
    "Pistachio_100+/uspto_unique_templates.csv.gz",
    "Pistachio_100+/uspto_ringbreaker_unique_templates.csv.gz",
]

folder_list = [
    "Pistachio_50", "Pistachio_25", "Pistachio_100", "Pistachio_100+"
]


def download_files_from_s3():
    s3_client = boto3.client(
        's3',
        aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        # aws_session_token=SESSION_TOKEN
    )

    for file, dest in zip(download_list, download_locations):

        s3_key = os.path.join(FOLDER_PATH, dest)
        print(f"Downloading {file} to {s3_key}")
        try:
            s3_client.download_file(S3_BUCKET, "models/" + file, s3_key)
            print(f"Downloaded s3://{S3_BUCKET}/models/{file} to {s3_key}")
        except FileNotFoundError:
            print(f"The file {file} was not found.")
        except NoCredentialsError:
            print("Credentials not available.")


def copy_other_files():
    uspto_path = FOLDER_PATH + "USPTO/"
    files = ["config.yml"]
    # files = ["uspto_filter_model.onnx", "zinc_stock.hd5", "config.yml"]
    # copy files from USPTO folder to folders in folder_list
    for folder in folder_list:
        for file in files:
            os.system(f"cp {uspto_path}{file} {FOLDER_PATH}{folder}/{file}")


def modify_config_files():
    for folder in folder_list:
        with open(FOLDER_PATH + folder + "/config.yml", "r") as f:
            lines = f.readlines()
        lines[2] = lines[2].replace("USPTO", folder)
        lines[2] = lines[2].replace("uspto_model.onnx", "uspto_expansion.onnx")
        lines[3] = lines[3].replace("USPTO", folder)
        lines[3] = lines[3].replace("uspto_templates.csv.gz",
                                    "uspto_unique_templates.csv.gz")
        lines[5] = lines[5].replace("USPTO", folder)
        lines[5] = lines[5].replace("uspto_ringbreaker_model.onnx",
                                    "uspto_ringbreaker_expansion.onnx")
        lines[6] = lines[6].replace("USPTO", folder)
        lines[6] = lines[6].replace(
            "uspto_ringbreaker_templates.csv.gz",
            "uspto_ringbreaker_unique_templates.csv.gz")

        with open(FOLDER_PATH + folder + "/config.yml", "w") as f:
            f.writelines(lines)


if __name__ == "__main__":
    download_files_from_s3()
    copy_other_files()
    modify_config_files()
