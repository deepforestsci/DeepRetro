import os
from dotenv import load_dotenv

# This will load variables from a .env file in the current directory (project root)
# or parent directories.
load_dotenv()

# Attempt to get the API_KEY from the environment
API_KEY_FROM_ENV = os.getenv('API_KEY')

# This is the debug print statement you uncommented in src/api.py
print(f"DEBUG: API_KEY loaded from environment: {API_KEY_FROM_ENV}")

if not API_KEY_FROM_ENV:
    print("--------------------------------------------------------------------")
    print("CRITICAL INFO: The 'API_KEY' environment variable was NOT found or is empty.")
    print("Please ensure:")
    print("  1. You have a .env file in your project root directory (recursiveLLM/).")
    print("     (The same directory where you are running this check_api_key.py script from)")
    print("  2. The .env file contains a line like: API_KEY='your_actual_secret_key'")
    print("  3. OR, you have set API_KEY as a system environment variable globally.")
    print("--------------------------------------------------------------------")
else:
    print("--------------------------------------------------------------------")
    print(f"SUCCESS: API_KEY found with value: '{API_KEY_FROM_ENV}'")
    print("This means python-dotenv and os.getenv are working as expected for API_KEY.")
    print("--------------------------------------------------------------------")

# You can also add the check from src/api.py to see its exact message
if not API_KEY_FROM_ENV:
    print("\nSimulating src/api.py exit condition:")
    print("CRITICAL ERROR: The 'API_KEY' environment variable is not set.")
    print("Please set this variable in your .env file or your system environment.")
    print("Example: API_KEY='your-chosen-secret-key'")
    # exit(1) # No need to exit in this test script
else:
    print("\nSimulating src/api.py successful load: API_KEY is present.") 