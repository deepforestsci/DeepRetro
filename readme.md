# DeepRetro - AI-Powered Retrosynthesis Tool

DeepRetro is an advanced retrosynthesis tool that combines Large Language Models (LLMs) with AiZynthFinder to predict chemical reaction pathways. It provides both a web-based interface for interactive analysis and a REST API for programmatic access.

## What is DeepRetro?

DeepRetro performs retrosynthesis analysis by:
- Taking a target molecule (in SMILES format) as input
- Using AI models (Claude 3, DeepSeek-R1) to predict possible precursor molecules
- Integrating with AiZynthFinder for reaction pathway validation
- Providing interactive visualization of reaction trees
- Supporting multiple model configurations and validation checks

## Quick Setup

### Prerequisites

- Python 3.9
- Conda or Miniconda
- Git

### Step 1: Clone and Setup

```bash
git clone <repository-url>
cd recursiveLLM

# Create conda environment
conda env create -f environment.yml
conda activate dfs_si_challenge
```

### Step 2: Install Models

```bash
# Create models directory
mkdir aizynthfinder/models

# Download USPTO models (free)
python -c "from aizynthfinder.utils.download_public_data import download_public_data; download_public_data('aizynthfinder/models/')"
```

### Step 3: Configure Environment

Create a `.env` file in the project root:

```bash
# Backend API Key (required)
API_KEY=your-secure-backend-api-key

# LLM API Keys (required for your chosen models)
ANTHROPIC_API_KEY=your-anthropic-key
OPENAI_API_KEY=your-openai-key
FIREWORKS_API_KEY=your-fireworks-key
```

### Step 4: Start the Application

```bash
# Start backend (in one terminal)
python src/api.py

# Start frontend (in another terminal)
cd viewer
python -m http.server 8000
```

### Step 5: Access the Application

1. Open your browser and go to `http://localhost:8000`
2. Enter your API key when prompted
3. Start analyzing molecules!

## Usage

### Web Interface

1. **Smart Retrosynthesis**: Enter a SMILES string and click "Analyze"
2. **View Pathway**: Upload existing JSON pathway files
3. **Advanced Settings**: Configure model types and validation flags
4. **Interactive Visualization**: Explore reaction pathways with molecular structures
5. **Partial Rerun**: Edit and rerun specific steps in the pathway
6. **File Upload**: Load and visualize previously generated pathways

### Key Features

- **Interactive Pathway Visualization**: D3.js-based tree visualization with molecular structures
- **Partial Rerun Analysis**: Edit any step and regenerate the pathway from that point
- **Multiple Model Support**: Switch between Claude 3, DeepSeek-R1, and GPT-4o models
- **Validation Checks**: Stability and hallucination detection for reliable results
- **File Management**: Upload, view, and edit JSON pathway files
- **Advanced Settings**: Configure model parameters and validation flags

### API Endpoints

- `POST /api/retrosynthesis` - Perform retrosynthesis analysis
- `POST /api/rerun_retrosynthesis` - Rerun complete analysis
- `POST /api/partial_rerun` - Rerun from specific step
- `GET /api/health` - Health check
- `POST /api/clear_molecule_cache` - Clear molecule cache

### Example API Request

```bash
curl -X POST http://localhost:5000/api/retrosynthesis \
  -H "Content-Type: application/json" \
  -H "X-API-KEY: your-api-key" \
  -d '{
    "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
    "model_type": "claude37",
    "advanced_prompt": true,
    "model_version": "USPTO",
    "stability_flag": true,
    "hallucination_check": true
  }'
```

## Model Configuration

### Supported LLM Models
- **Claude 3 Opus**: `claude-3-opus-20240229`
- **Claude 3.5 Sonnet**: `anthropic/claude-3-5-sonnet-20241022`
- **DeepSeek-R1**: `fireworks_ai/accounts/fireworks/models/deepseek-r1`
- **GPT-4o**: `gpt-4o`
- **GPT-4o Mini**: `gpt-4o-mini`

### Supported AiZynthFinder Models
- **USPTO**: Standard USPTO reaction database (Free, publicly available)
- **Pistachio_25**: Pistachio model with 25% coverage (Requires permissions)
- **Pistachio_50**: Pistachio model with 50% coverage (Requires permissions)
- **Pistachio_100**: Pistachio model with 100% coverage (Requires permissions)
- **Pistachio_100+**: Enhanced Pistachio model (Requires permissions)

> **Note**: Only the USPTO model is freely available. Pistachio models require special access permissions.

### Validation Features
- **Stability Check**: Validates molecular stability
- **Hallucination Check**: Detects AI-generated invalid structures
- **Advanced Prompt**: Uses enhanced prompting for better results

## Troubleshooting

### Common Issues

1. **API Key Errors**
   - Ensure backend `.env` file has `API_KEY` set
   - Verify frontend is using the same API key

2. **Model Download Issues**
   - Ensure you have internet connection for model downloads
   - Check `aizynthfinder/models/` directory exists

3. **Port Conflicts**
   - Backend default: 5000
   - Frontend default: 8000

4. **LLM API Errors**
   - Verify API keys are valid and have sufficient credits
   - Check rate limits for your chosen LLM provider

## Development

### Running Tests

```bash
pip install -r tests/requirements_tests.txt
python -m pytest tests/
```

## License

This project is licensed under the MIT License.

## Contributing

We welcome contributions! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.