# DeepRetro Testing Suite Documentation

## Overview

This document provides a detailed explanation of the testing suite in the DeepRetro project, focusing on `test_api.py` and `test_cache.py`. The testing suite is designed to ensure the proper functioning of the API endpoints and caching mechanisms that are critical for the retrosynthesis application.

## Table of Contents

1. [Testing Environment Setup](#testing-environment-setup)
2. [Test API (test_api.py)](#test-api)
    - [API Functions Testing](#api-functions-testing)
    - [API Endpoints Testing](#api-endpoints-testing)
    - [Parameter Handling Tests](#parameter-handling-tests)
    - [Error Handling Tests](#error-handling-tests)
3. [Test Cache (test_cache.py)](#test-cache)
    - [Cache Key Generation Testing](#cache-key-generation-testing)
    - [Cache Results Decorator Testing](#cache-results-decorator-testing)
    - [Cache Clearing Functionality Testing](#cache-clearing-functionality-testing)
4. [Integration Between API and Cache](#integration-between-api-and-cache)
5. [The Chemistry Behind the Tests](#the-chemistry-behind-the-tests)
    - [The Retrosynthesis Workflow](#the-retrosynthesis-workflow)
    - [Role of Caching in Retrosynthesis](#role-of-caching-in-retrosynthesis)
    - [API's Role in the Synthesis Planning Process](#apis-role-in-the-synthesis-planning-process)
6. [Testing Methodology](#testing-methodology)

## Testing Environment Setup

The testing suite uses pytest and unittest frameworks. Before running tests, an environment can be created from `tests/test_env.yml` which contains all the necessary dependencies.

To run tests:
```
pytest tests/test_api.py  # Run without verbosity
pytest tests/test_api.py -v  # Run with verbosity
```

The tests leverage mocking extensively to isolate components and test them independently without requiring actual external resources like LLMs or database connections.

## Test API

The `test_api.py` file tests the Flask API endpoints and functions defined in `src/api.py`. It uses unittest and mock objects to simulate API requests and isolate the API functions from external dependencies.

### API Functions Testing

#### File Operations

- **`test_save_result_creates_file_with_correct_content`**: Verifies that the `save_result` function creates a file with the correct content structure.
- **`test_load_result_file_not_exists`**: Tests the behavior of `load_result` when the file does not exist.
- **`test_load_result_file_exists_and_valid`**: Ensures `load_result` correctly loads valid JSON data.
- **`test_load_result_file_malformed_json`**: Tests error handling when the file contains malformed JSON.

### API Endpoints Testing

The test suite thoroughly tests all API endpoints:

#### Health Endpoint

- **`test_health_endpoint_no_api_key`**: Verifies 401 status when no API key is provided.
- **`test_health_endpoint_wrong_api_key`**: Verifies 401 status when incorrect API key is provided.
- **`test_health_endpoint_correct_api_key`**: Verifies 200 status and proper response with correct API key.
- **`test_health_endpoint_options_request_bypasses_auth`**: Tests that OPTIONS requests bypass authentication as required for CORS.

#### Clear Molecule Cache Endpoint

- **`test_clear_molecule_cache_no_api_key`**: Verifies authorization is required.
- **`test_clear_molecule_cache_no_molecule_param`**: Tests handling of missing parameters.
- **`test_clear_molecule_cache_success`**: Verifies proper behavior when all requirements are met.

#### Retrosynthesis Endpoint

- **`test_retrosynthesis_no_api_key`**: Verifies authorization is required.
- **`test_retrosynthesis_no_smiles_param`**: Tests handling when SMILES is missing.
- **`test_retrosynthesis_invalid_smiles_string`**: Verifies validation of SMILES input.
- **`test_retrosynthesis_default_parameters`**: Tests behavior with default parameters.
- **`test_retrosynthesis_custom_parameters_valid`**: Tests behavior with custom valid parameters.
- **`test_retrosynthesis_invalid_model_type_falls_back_to_default`**: Verifies fallback behavior for invalid model type.
- **`test_retrosynthesis_advanced_prompt_false_string`**: Tests handling of advanced prompt flag as string.
- **`test_retrosynthesis_main_function_exception_handling`**: Tests error handling for exceptions in the main function.

#### Rerun Retrosynthesis Endpoint

- Tests for authentication, parameter validation, and SMILES validation
- Tests for successful execution with default and custom parameters
- Tests for exception handling

#### Partial Rerun Endpoint

- Tests for authentication, parameter validation
- Tests for various scenarios like missing stored results, SMILES mismatch, step not found
- Tests for target step validation (missing products, invalid product format)
- Tests for different types of partial reruns (middle, root, leaf steps)
- Tests for handling empty sub-synthesis results and exceptions

### Parameter Handling Tests

The test suite verifies that the API correctly handles various parameters:

- **Model Type**: Handles "claude3", "claude37", and "deepseek", with appropriate fallback
- **Advanced Prompt**: Properly processes boolean strings "True"/"False"
- **Model Version**: Validates against `AZ_MODEL_LIST`
- **Stability Flag**: Processes "True"/"False" strings
- **Hallucination Check**: Processes "True"/"False" strings

### Error Handling Tests

Tests ensure that the API provides appropriate error responses for:
- Authentication failures
- Invalid or missing parameters
- Invalid SMILES strings
- Exceptions in the core processing logic

## Test Cache

The `test_cache.py` file tests the caching mechanism defined in `src/cache.py`. It uses pytest fixtures and mocking to create an isolated environment for testing the cache functionality.

### Cache Key Generation Testing

- **`test_generate_cache_key_format_and_consistency`**: Verifies that cache keys have the expected format and are consistently generated for the same inputs.
- **`test_generate_cache_key_arg_sensitivity`**: Tests that different inputs produce different cache keys, including changes to positional args, kwargs values, and function names.

### Cache Results Decorator Testing

- **`test_cache_results_decorator_cache_miss_and_hit`**: Verifies that functions are only executed once for the same inputs and subsequent calls use cached results.
- **`test_cache_results_decorator_different_args`**: Tests that different arguments result in different cache entries and function calls.

### Cache Clearing Functionality Testing

- **`test_clear_entire_cache`**: Tests that the entire cache can be cleared.
- **`test_clear_cache_for_molecule`**: Tests selective clearing of cache entries related to a specific molecule.
- **`test_clear_cache_for_molecule_no_matches`**: Tests behavior when no matching entries are found.
- **`test_clear_cache_for_molecule_empty_cache`**: Tests behavior with an empty cache.

## Integration Between API and Cache

The tests demonstrate the integration between the API and cache:

1. The API uses the cache decorator to store results of expensive computations
2. The API provides endpoints to selectively clear the cache
3. When rerunning, the API can clear cached results to ensure fresh computations

## The Chemistry Behind the Tests

### The Retrosynthesis Workflow

The DeepRetro project implements a computer-aided retrosynthesis workflow that combines two powerful technologies:

1. **AiZynthFinder**: A rule-based retrosynthesis tool that applies known chemical transformations to break down complex molecules into simpler precursors
2. **Large Language Models (LLMs)**: Models like Claude and DeepSeek that bring chemical reasoning and knowledge to suggest synthesis pathways

The testing suite ensures that this hybrid approach functions correctly by verifying:

- **SMILES Validation**: Tests confirm that chemical structures are properly validated using RDKit's `MolFromSmiles` function before processing
- **Model Selection**: Tests verify that the appropriate AI models (both AiZynthFinder model versions and LLMs) can be selected and configured
- **Parameter Handling**: Tests ensure stability and hallucination checks can be enabled to improve synthesis route quality

The retrosynthesis process works as follows:

1. A target molecule (represented as a SMILES string) is submitted to the API
2. The API validates the structure and selects appropriate models
3. AiZynthFinder suggests possible reaction pathways based on transformation rules
4. The LLM evaluates these pathways, suggests improvements, and provides chemical reasoning
5. Results are returned, showing a hierarchical synthesis route from commercially available starting materials to the target

### Role of Caching in Retrosynthesis

The caching system plays a crucial role in making the retrosynthesis process efficient:

1. **Computational Efficiency**: Retrosynthesis calculations, especially for complex molecules, are computationally expensive. The caching system prevents redundant calculations for molecules that have already been analyzed.

2. **LLM Call Optimization**: API calls to large language models are both expensive and time-consuming. The `cache_results` decorator is applied to functions like:
   - `call_LLM`: Caches responses from language models for specific molecule inputs
   - `run_az_with_img`: Caches AiZynthFinder results for each molecule
   - `validity_check`: Caches the validity verification of proposed reaction pathways

3. **Molecule-Specific Cache Clearing**: The `clear_cache_for_molecule` function allows selective clearing of cached results for a specific molecule. This is important when:
   - New models or parameters become available
   - Previous synthesis attempts had errors
   - Users want to explore alternative routes

The tests verify this caching behavior by ensuring:
- Cache hits return stored results without recomputation
- Cache misses trigger actual function execution
- Cache keys properly differentiate between different inputs
- Cache entries can be selectively cleared

### API's Role in the Synthesis Planning Process

The API serves as the interface between users and the retrosynthesis engine, with several key functions tested:

1. **Full Retrosynthesis Planning**: The `/api/retrosynthesis` endpoint handles full synthesis planning for a target molecule, with tests verifying:
   - Parameter customization (model selection, stability checks, etc.)
   - Proper execution of the main retrosynthesis function
   - Storage of results for later partial reruns

2. **Selective Reprocessing**: The `/api/rerun_retrosynthesis` endpoint allows clearing the cache and rerunning the entire synthesis, important when new models or parameters are available.

3. **Partial Reruns**: The `/api/partial_rerun` endpoint enables targeted reprocessing of specific steps in a synthesis route. Tests verify:
   - Different rerun scenarios (middle steps, leaf nodes, root nodes)
   - Proper integration of new sub-synthesis results into the existing synthesis tree
   - Preservation of unmodified branches of the synthesis tree

This sophisticated API structure allows chemists to iteratively refine synthesis plans by:
- Starting with a complete retrosynthesis
- Identifying problematic steps or suboptimal routes
- Selectively rerunning those specific parts with different parameters
- Preserving satisfactory portions of the synthesis plan

## Testing Methodology

The test suite demonstrates several testing best practices:

1. **Isolated Testing**: Using mocks to isolate components
2. **Comprehensive Coverage**: Testing both happy paths and error conditions
3. **Parameter Validation**: Ensuring all parameters are properly validated
4. **Error Handling**: Verifying appropriate error responses
5. **Integration Testing**: Testing how components work together

The testing architecture reflects the application's design:

- **API Layer**: Handles HTTP requests, parameter validation, and response formatting
- **Caching Layer**: Optimizes performance by storing results of expensive operations
- **Core Business Logic**: The retrosynthesis functionality (tested indirectly through API tests)

This modular approach allows for targeted testing and makes it easier to identify issues when tests fail. 