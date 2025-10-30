# Agentic LLM Pipeline Implementation Summary

## Overview

Successfully implemented an agentic LLM pipeline with Anthropic tool calling API and comprehensive Langfuse logging. The system now uses iterative tool validation to reduce hallucinations and improve SMILES accuracy.

## What Was Implemented

### 1. Tools Module (`src/utils/tools.py`)
- **SMILES Validator Tool**: Validates SMILES strings using RDKit and returns canonical forms
- **Tool Registry**: Centralized management of tool definitions and executors
- **Tool Dispatcher**: Executes tools by name with error handling
- Follows Anthropic's tool definition format

### 2. Agentic LLM Functions (`src/utils/llm.py`)
- **`agentic_call_LLM()`**: New function implementing agent loop with tool calling
  - Supports up to `max_iterations` (default: 5) iterations
  - Handles tool use blocks from Claude
  - Executes tools and feeds results back to LLM
  - Tracks detailed tool usage statistics
- **`agentic_llm_pipeline()`**: Pipeline wrapper maintaining retry logic
  - Compatible with existing `llm_pipeline()` interface
  - Returns additional tool statistics
  - Supports all existing features (stability checks, hallucination checks, etc.)

### 3. Langfuse Integration
Direct Langfuse SDK integration for comprehensive logging:
- **Trace Level**: Overall retrosynthesis job with molecule and parameters
- **Iteration Spans**: Each agent iteration with timing and tool counts
- **Generation Logs**: LLM calls with token usage and responses
- **Tool Spans**: Individual tool executions with inputs and outputs
- **Event Logs**: Final statistics and completion events

All Langfuse errors are caught gracefully to prevent pipeline failures.

### 4. Prompt Engineering (`src/variables.py`)
Updated all system and user prompts to guide tool usage:
- Instructions to use `validate_smiles` tool before finalizing predictions
- Workflow guidance for iterative validation and refinement
- Applied to: `SYS_PROMPT`, `USER_PROMPT`, `SYS_PROMPT_V4`, `USER_PROMPT_V4`, `SYS_PROMPT_OPENAI`, `SYS_PROMPT_DEEPSEEK`

### 5. Full System Integration
- **`rec_prithvi.py`**: Updated recursive retrosynthesis to use agentic pipeline
  - New parameters: `use_agentic_pipeline` (default: True), `max_tool_iterations` (default: 5)
  - Logs tool statistics at each recursion level
- **`prithvi.py`**: Updated entry point to pass through agentic parameters
- **`main.py`**: Updated main function with new parameters
- **`api.py`**: Updated all API endpoints:
  - `/api/retrosynthesis`
  - `/api/rerun_retrosynthesis`
  - `/api/partial_rerun`

### 6. Configuration (`config/advanced_settings.json`)
Added new configuration sections:
```json
{
  "defaults": {
    "use_agentic_pipeline": true,
    "max_tool_iterations": 5
  },
  "agentic_settings": {
    "enable_tools": true,
    "available_tools": ["validate_smiles"]
  },
  "langfuse_settings": {
    "enable_tool_logging": true
  }
}
```

### 7. Testing
Created comprehensive test suites:
- **`tests/test_tools.py`**: Tests for SMILES validator and tool infrastructure
- **`tests/test_agentic_llm.py`**: Tests for agentic LLM functions with mocking

## How to Use

### API Usage

#### Standard Retrosynthesis with Agentic Pipeline (Default)
```json
POST /api/retrosynthesis
{
  "smiles": "CCCO",
  "model_type": "claude4sonnet",
  "use_agentic_pipeline": true,
  "max_tool_iterations": 5
}
```

#### Disable Agentic Pipeline (Use Original)
```json
POST /api/retrosynthesis
{
  "smiles": "CCCO",
  "use_agentic_pipeline": false
}
```

#### Advanced Configuration
```json
POST /api/retrosynthesis
{
  "smiles": "CCCO",
  "model_type": "claude4opus",
  "advanced_prompt": true,
  "stability_flag": true,
  "hallucination_check": true,
  "use_agentic_pipeline": true,
  "max_tool_iterations": 10
}
```

### Python Usage

```python
from src.main import main

# With agentic pipeline (default)
result = main(
    smiles="CCCO",
    llm="claude-opus-4-20250514",
    use_agentic_pipeline=True,
    max_tool_iterations=5
)

# Without agentic pipeline
result = main(
    smiles="CCCO",
    use_agentic_pipeline=False
)
```

### Direct Function Usage

```python
from src.utils.llm import agentic_llm_pipeline

pathways, explanations, confidence, tool_stats = agentic_llm_pipeline(
    molecule="CCCO",
    LLM="claude-opus-4-20250514",
    max_tool_iterations=5,
    trace_id="my_trace_id"  # Optional: for Langfuse hierarchical logging
)

print(f"Tool statistics: {tool_stats}")
# Output: {
#   "total_iterations": 2,
#   "total_tool_calls": 3,
#   "successful_validations": 3,
#   "failed_validations": 0,
#   "tool_calls_by_name": {"validate_smiles": 3}
# }
```

## Benefits

### 1. Reduced Hallucinations
- LLM validates each proposed SMILES before returning
- Invalid SMILES are caught and corrected iteratively
- Canonical forms ensure consistency

### 2. Self-Correction
- Agent can iterate to fix validation errors
- Max iterations prevent infinite loops
- Graceful degradation on repeated failures

### 3. Better Accuracy
- Validated SMILES ensure chemical validity
- RDKit canonicalization provides standard forms
- Fewer downstream errors from invalid structures

### 4. Observability
- Comprehensive Langfuse logging shows:
  - Which tools were called and when
  - What inputs and outputs tools received
  - How many iterations were needed
  - Token usage per iteration
- Helps debug and improve prompts

### 5. Extensibility
- Easy to add new tools:
  1. Define tool schema in `src/utils/tools.py`
  2. Implement execution function
  3. Register in `TOOL_DEFINITIONS` and `TOOL_EXECUTORS`
- Tools automatically available to all models

## Langfuse Dashboard

To view tool usage in Langfuse:
1. Go to your Langfuse dashboard
2. Find traces named "agentic_retrosynthesis"
3. Expand to see:
   - Iteration spans showing agent loops
   - Tool call spans with validate_smiles executions
   - Input SMILES and validation results
   - Token usage and timing

## Backwards Compatibility

The system maintains full backwards compatibility:
- Default behavior is agentic (can be disabled)
- All existing parameters still work
- API responses have the same format
- Caching automatically differentiates tool-enhanced vs standard calls

## Performance Considerations

- **Latency**: Agentic calls take longer due to tool iterations (typically 2-3 iterations)
- **Token Usage**: Higher token usage from conversation turns with tools
- **Caching**: Helps reduce costs for repeated molecules
- **Max Iterations**: Set lower (e.g., 3) for faster responses

## Future Enhancements

Potential tools to add:
1. **Atom Balance Checker**: Verify atom conservation
2. **Reaction Feasibility Checker**: Validate proposed reactions
3. **Functional Group Analyzer**: Identify reactive groups
4. **Stereochemistry Validator**: Check stereochemical consistency
5. **Literature Reaction Lookup**: Find similar known reactions
6. **Commercial Availability Checker**: Verify precursor availability

## Testing

Run the test suites:
```bash
# Activate virtual environment
source .venv/bin/activate

# Run all tests
pytest tests/test_tools.py tests/test_agentic_llm.py -v

# Run specific test
pytest tests/test_tools.py::TestSMILESValidatorTool::test_valid_smiles -v
```

## Environment Variables

Add to `.env` if needed:
```bash
# Enable/disable tool use (default: true)
ENABLE_TOOL_USE=true

# Max tool iterations (default: 5)
MAX_TOOL_ITERATIONS=5

# Enable Langfuse tool logging (default: true)
LANGFUSE_ENABLE_TOOL_LOGGING=true
```

## Files Modified

### New Files
- `src/utils/tools.py` - Tool definitions and executors
- `tests/test_tools.py` - Tool tests
- `tests/test_agentic_llm.py` - Agentic pipeline tests

### Modified Files
- `src/utils/llm.py` - Added agentic functions
- `src/variables.py` - Updated prompts for tool usage
- `src/rec_prithvi.py` - Integration with recursive retrosynthesis
- `src/prithvi.py` - Pass through agentic parameters
- `src/main.py` - Updated main function
- `src/api.py` - Updated API endpoints
- `config/advanced_settings.json` - Added agentic configuration

### Unchanged Files (Backwards Compatible)
- `src/cache.py` - Automatic cache differentiation
- `src/utils/utils_molecule.py` - Reused for validation
- `src/utils/az.py` - No changes needed
- All other utility modules

## Notes

1. **Langfuse SDK Required**: Ensure `langfuse` Python package is installed
   ```bash
   pip install langfuse
   ```

2. **Claude Models**: Works best with Claude 3.5+ models that support tool calling

3. **Error Handling**: All Langfuse logging failures are caught to prevent pipeline disruption

4. **Cost**: Agentic calls cost more due to additional iterations and token usage

## Support

For issues or questions:
1. Check Langfuse dashboard for tool execution logs
2. Review test files for usage examples
3. Set `use_agentic_pipeline=False` to test without tools
4. Check log files in `logs/` directory for detailed execution traces

