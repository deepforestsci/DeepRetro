# DeepRetro: A Recursive Retrosynthesis LLM Project

<!-- TODO: Add a concise 1-2 paragraph summary of the project, its goals, and key features. -->

**Full Documentation:** [Read the Docs](YOUR_READTHEDOCS_LINK_HERE) <!-- Replace with actual link once available -->

<!-- Badges will go here: PyPI, Build Status, Docs Status, License -->

## ToDos
- [x] Add Flask API for retrosynthesis predictions
- [ ] Change logging style to jobID based logging
- [ ] Add Dockerfile for easy deployment
- [ ] Add tests for the retrosynthesis predictions
- [ ] Add CI/CD pipeline for the project
- [ ] Add documentation for the project
- [ ] Add clean logging for the project


## Getting Started: Self-Hosting DeepRetro

For detailed installation, advanced usage, and API information, please see our [full documentation](YOUR_READTHEDOCS_LINK_HERE).

These instructions will guide you through setting up and running your own instance of DeepRetro locally.

### Prerequisites

*   Python (version 3.9 or as specified in `environment.yml`)
*   Conda (for managing Python environments)
*   Access to an internet connection for downloading dependencies and models.
*   Your own API keys for Large Language Model (LLM) services if you intend to use models from Anthropic, OpenAI, etc. (see Configuration section).

### 1. Clone the Repository

First, clone this repository to your local machine:

```bash
git clone https://github.com/your-username/recursiveLLM.git # Replace with your repo URL
cd recursiveLLM
```

### 2. Set Up Python Environment

We use Conda to manage dependencies. Create and activate the environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate dfs_si_challenge
```
This command will install all necessary Python packages, including Flask, AiZynthFinder, RDKit, and libraries for LLM interaction.

### 3. Configure API Keys and Application Settings

This application requires API keys to function correctly. Some are for external services (like LLMs), and one is for securing your local backend instance.

1.  **Create a `.env` file:**
    In the root directory of the project, copy the template file `.env.template` to a new file named `.env`:
    ```bash
    cp .env.template .env
    ```

2.  **Edit your `.env` file:**
    Open the newly created `.env` file with a text editor and fill in the required values:

    *   **`API_KEY`**: This is a secret key **you choose** to protect your local backend server. Any string of characters will work, but make it reasonably complex (e.g., `myRandomSecretKeyForDeepRetro123!`). Your frontend will use this key to authenticate with your backend.
    *   **`ANTHROPIC_API_KEY`**: If you plan to use Anthropic's Claude models, enter your Anthropic API key here.
    *   **`OPENAI_API_KEY`**: If you plan to use OpenAI models (like GPT), enter your OpenAI API key here.
    *   **`LANGFUSE_SECRET_KEY`**, **`LANGFUSE_PUBLIC_KEY`**, **`LANGFUSE_HOST`**: If you use Langfuse for tracing and observability, provide your Langfuse credentials.
    *   Review other variables like `AZ_MODEL_CONFIG_PATH` and `RXN_CLASSIFICATION_MODEL_PATH` and ensure they point to the correct locations if you have custom model paths. The defaults should generally work with the provided structure.

    **Example `.env` file content:**
    ```env
    API_KEY='your_chosen_secret_for_local_backend'
    ANTHROPIC_API_KEY='sk-ant-xxxxxxxxxxxxxxxxxxxx'
    OPENAI_API_KEY='sk-xxxxxxxxxxxxxxxxxxxx'
    LANGFUSE_SECRET_KEY='ls_sk_xxxxxxxxxxxx'
    LANGFUSE_PUBLIC_KEY='ls_pk_xxxxxxxxxxxx'
    LANGFUSE_HOST='https://us.cloud.langfuse.com'
    AZ_MODEL_CONFIG_PATH='aizynthfinder/models/config.yml'
    RXN_CLASSIFICATION_MODEL_PATH='reaction_prediction/rfc.pkl'
    ```
    **Important:** Do NOT commit your `.env` file to Git. It's already listed in `.gitignore`.

### 4. Start the Backend Server

With your Conda environment (`dfs_si_challenge`) activated and your `.env` file configured, start the Flask backend server:

```bash
python src/api.py
```
You should see output indicating the server is running (e.g., on `http://127.0.0.1:5000/`). The server will also print a debug message confirming the `API_KEY` it loaded from your `.env` file.

### 5. Access and Configure the Frontend

1.  **Open the application:** Open the `viewer/index.html` file in your web browser.
2.  **Enter API Key (First Time Only):**
    *   When you first open `index.html`, your browser will display a prompt asking: "Please enter the API Key for your self-hosted backend..."
    *   Enter the **exact same `API_KEY`** that you defined in your `.env` file in Step 3 (e.g., `'your_chosen_secret_for_local_backend'`).
    *   Click "OK".
    *   The frontend will store this key in your browser's local storage, so you won't be prompted again unless you clear your browser's storage for this page.

### 6. Use DeepRetro

You should now be able to use the DeepRetro application! Enter a SMILES string and click "Analyze" to start a retrosynthesis prediction.

**Troubleshooting:**
*   **`ModuleNotFoundError` for `aizynthfinder` (or others):** Ensure your Conda environment (`dfs_si_challenge`) is activated and that the `conda env create -f environment.yml` command completed successfully. You might need to force reinstall a problematic package using pip within the Conda environment (e.g., `pip install --upgrade --force-reinstall aizynthfinder`).
*   **Backend `API_KEY` not found error:** Double-check that your `.env` file is in the project root, is named correctly (`.env`, not `env` or `.env.txt`), and contains the `API_KEY='yourkey'` line.
*   **Frontend shows errors or can't connect:**
    *   Ensure your backend server (`python src/api.py`) is running.
    *   Verify the API key you entered in the frontend prompt exactly matches the `API_KEY` in your `.env` file.
    *   Check the browser's developer console (F12) for any error messages.