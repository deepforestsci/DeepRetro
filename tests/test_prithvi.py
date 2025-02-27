import pytest
from unittest.mock import patch, MagicMock
from src.prithvi import run_prithvi, add_metadata

@pytest.fixture
def mock_context_logger():
    with patch('src.prithvi.context_logger') as mock_logger:
        yield mock_logger

@pytest.fixture
def mock_structlog():
    with patch('src.prithvi.structlog') as mock_log:
        yield mock_log

@pytest.fixture
def mock_add_job_specific_handler():
    with patch('src.prithvi.add_job_specific_handler') as mock_handler:
        yield mock_handler

@pytest.fixture
def mock_rec_run_prithvi():
    with patch('src.prithvi.rec_run_prithvi') as mock_run:
        yield mock_run

@pytest.fixture
def mock_format_output():
    with patch('src.prithvi.format_output') as mock_format:
        yield mock_format

@pytest.fixture
def mock_reagent_agent():
    with patch('src.prithvi.reagent_agent') as mock_agent:
        yield mock_agent

@pytest.fixture
def mock_conditions_agent():
    with patch('src.prithvi.conditions_agent') as mock_agent:
        yield mock_agent

@pytest.fixture
def mock_literature_agent():
    with patch('src.prithvi.literature_agent') as mock_agent:
        yield mock_agent

def test_run_prithvi_success(mock_context_logger, mock_structlog, mock_add_job_specific_handler, mock_rec_run_prithvi, mock_format_output):
    mock_rec_run_prithvi.return_value = ({'steps': []}, None)
    mock_format_output.return_value = {'steps': []}
    mock_add_job_specific_handler.return_value = MagicMock()

    result = run_prithvi(molecule="CCO")

    assert result == {'steps': []}
    mock_rec_run_prithvi.assert_called_once()
    mock_format_output.assert_called_once()

def test_run_prithvi_custom_llm(mock_context_logger, mock_structlog, mock_add_job_specific_handler, mock_rec_run_prithvi, mock_format_output):
    mock_rec_run_prithvi.return_value = ({'steps': []}, None)
    mock_format_output.return_value = {'steps': []}
    mock_add_job_specific_handler.return_value = MagicMock()

    result = run_prithvi(molecule="CCO", llm="custom-llm")

    assert result == {'steps': []}
    mock_rec_run_prithvi.assert_called_once_with(molecule="CCO", job_id=pytest.ANY, llm="custom-llm", az_model="USPTO")

def test_run_prithvi_custom_az_model(mock_context_logger, mock_structlog, mock_add_job_specific_handler, mock_rec_run_prithvi, mock_format_output):
    mock_rec_run_prithvi.return_value = ({'steps': []}, None)
    mock_format_output.return_value = {'steps': []}
    mock_add_job_specific_handler.return_value = MagicMock()

    result = run_prithvi(molecule="CCO", az_model="custom-az-model")

    assert result == {'steps': []}
    mock_rec_run_prithvi.assert_called_once_with(molecule="CCO", job_id=pytest.ANY, llm="claude-3-opus-20240229", az_model="custom-az-model")

def test_run_prithvi_exception_handling(mock_context_logger, mock_structlog, mock_add_job_specific_handler, mock_rec_run_prithvi, mock_format_output):
    mock_rec_run_prithvi.side_effect = Exception("Test Exception")
    mock_add_job_specific_handler.return_value = MagicMock()

    with pytest.raises(Exception, match="Test Exception"):
        run_prithvi(molecule="CCO")

    mock_context_logger.reset.assert_called_once()

def test_add_metadata(mock_reagent_agent, mock_conditions_agent, mock_literature_agent):
    mock_reagent_agent.return_value = (True, ['reagent1'])
    mock_conditions_agent.return_value = (True, 'condition1')
    mock_literature_agent.return_value = (True, 'literature1')

    output_data = {
        'steps': [
            {
                'reactants': ['reactant1'],
                'products': ['product1'],
                'reagents': [],
                'reactionmetrics': [{}]
            }
        ]
    }

    result = add_metadata(output_data)

    assert result['steps'][0]['reagents'] == ['reagent1']
    assert result['steps'][0]['conditions'] == 'condition1'
    assert result['steps'][0]['reactionmetrics'][0]['closestliterature'] == 'literature1'