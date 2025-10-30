"""Unit tests for agentic LLM pipeline"""
import pytest
from unittest.mock import Mock, patch, MagicMock
from src.utils.llm import agentic_call_LLM, agentic_llm_pipeline


class TestAgenticCallLLM:
    """Test cases for agentic_call_LLM function"""
    
    @patch('src.utils.llm.completion')
    @patch('src.utils.llm.execute_tool')
    def test_simple_completion_no_tools(self, mock_execute_tool, mock_completion):
        """Test agentic call that completes without using tools"""
        # Mock LLM response without tool use
        mock_response = Mock()
        mock_response.choices = [Mock()]
        mock_response.choices[0].message = Mock()
        mock_response.choices[0].message.content = """<cot>
<thinking>Analysis here</thinking>
</cot>

<json>
{
  "data": [["CCO", "CC=O"]],
  "explanation": ["Oxidation"],
  "confidence_scores": [0.9]
}
</json>"""
        mock_response.choices[0].finish_reason = "end_turn"
        mock_response.choices[0].message.tool_calls = None
        mock_response.usage = Mock(input_tokens=100, output_tokens=50, total_tokens=150)
        mock_completion.return_value = mock_response
        
        # Call function
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_iterations=3
        )
        
        # Assertions
        assert status_code == 200
        assert "json" in res_text.lower()
        assert tool_stats["total_iterations"] == 1
        assert tool_stats["total_tool_calls"] == 0
    
    @patch('src.utils.llm.completion')
    @patch('src.utils.llm.execute_tool')
    def test_with_tool_use(self, mock_execute_tool, mock_completion):
        """Test agentic call that uses tools"""
        # Mock tool execution result
        mock_execute_tool.return_value = {
            "success": True,
            "valid": True,
            "canonical_smiles": "CCO",
            "error": None,
            "context": "unknown"
        }
        
        # First response: requests tool use
        mock_response1 = Mock()
        mock_response1.choices = [Mock()]
        mock_response1.choices[0].message = Mock()
        mock_response1.choices[0].message.content = "Let me validate this SMILES"
        mock_response1.choices[0].finish_reason = "tool_use"
        
        # Create mock tool call
        mock_tool_call = Mock()
        mock_tool_call.function = Mock()
        mock_tool_call.function.name = "validate_smiles"
        mock_tool_call.function.arguments = '{"smiles": "CCO"}'
        mock_tool_call.id = "call_123"
        mock_response1.choices[0].message.tool_calls = [mock_tool_call]
        mock_response1.usage = Mock(input_tokens=100, output_tokens=50, total_tokens=150)
        
        # Second response: final answer
        mock_response2 = Mock()
        mock_response2.choices = [Mock()]
        mock_response2.choices[0].message = Mock()
        mock_response2.choices[0].message.content = """<cot>
<thinking>Validated successfully</thinking>
</cot>

<json>
{
  "data": [["CCO", "CC=O"]],
  "explanation": ["Oxidation"],
  "confidence_scores": [0.9]
}
</json>"""
        mock_response2.choices[0].finish_reason = "end_turn"
        mock_response2.choices[0].message.tool_calls = None
        mock_response2.usage = Mock(input_tokens=120, output_tokens=60, total_tokens=180)
        
        mock_completion.side_effect = [mock_response1, mock_response2]
        
        # Call function
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_iterations=3
        )
        
        # Assertions
        assert status_code == 200
        assert tool_stats["total_iterations"] == 2
        assert tool_stats["total_tool_calls"] == 1
        assert tool_stats["successful_validations"] == 1
        assert tool_stats["failed_validations"] == 0
        assert mock_execute_tool.called
    
    @patch('src.utils.llm.completion')
    def test_max_iterations_reached(self, mock_completion):
        """Test that max iterations limit is respected"""
        # Mock responses that always request tools
        mock_response = Mock()
        mock_response.choices = [Mock()]
        mock_response.choices[0].message = Mock()
        mock_response.choices[0].message.content = "Using tools"
        mock_response.choices[0].finish_reason = "tool_use"
        
        mock_tool_call = Mock()
        mock_tool_call.function = Mock()
        mock_tool_call.function.name = "validate_smiles"
        mock_tool_call.function.arguments = '{"smiles": "CCO"}'
        mock_tool_call.id = "call_123"
        mock_response.choices[0].message.tool_calls = [mock_tool_call]
        mock_response.usage = Mock(input_tokens=100, output_tokens=50, total_tokens=150)
        
        mock_completion.return_value = mock_response
        
        # Call function with low max_iterations
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_iterations=2
        )
        
        # Should reach max iterations
        assert tool_stats["total_iterations"] == 2
    
    @patch('src.utils.llm.completion')
    def test_error_handling(self, mock_completion):
        """Test error handling in agentic call"""
        # Mock LLM to raise exception
        mock_completion.side_effect = Exception("API Error")
        
        # Call function
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_iterations=3
        )
        
        # Should return error status
        assert status_code == 400


class TestAgenticLLMPipeline:
    """Test cases for agentic_llm_pipeline function"""
    
    @patch('src.utils.llm.agentic_call_LLM')
    @patch('src.utils.llm.validity_check')
    def test_successful_pipeline(self, mock_validity_check, mock_agentic_call):
        """Test successful execution of agentic pipeline"""
        # Mock agentic call response
        mock_agentic_call.return_value = (
            200,
            """<cot>
<thinking>Analysis</thinking>
</cot>

<json>
{
  "data": [["CCO", "CC=O"]],
  "explanation": ["Oxidation of ethanol"],
  "confidence_scores": [0.9]
}
</json>""",
            {
                "total_iterations": 1,
                "total_tool_calls": 1,
                "successful_validations": 1,
                "failed_validations": 0,
                "tool_calls_by_name": {"validate_smiles": 1}
            }
        )
        
        # Mock validity check
        mock_validity_check.return_value = (
            [["CCO", "CC=O"]],
            ["Oxidation of ethanol"],
            [0.9]
        )
        
        # Call pipeline
        pathways, explanations, confidence, tool_stats = agentic_llm_pipeline(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_tool_iterations=5
        )
        
        # Assertions
        assert len(pathways) == 1
        assert len(explanations) == 1
        assert len(confidence) == 1
        assert tool_stats["total_tool_calls"] == 1
        assert tool_stats["successful_validations"] == 1
    
    @patch('src.utils.llm.agentic_call_LLM')
    @patch('src.utils.llm.validity_check')
    def test_retry_on_failure(self, mock_validity_check, mock_agentic_call):
        """Test pipeline retries on failure"""
        # First call fails, second succeeds
        mock_agentic_call.side_effect = [
            (400, "", {"total_iterations": 0, "total_tool_calls": 0, "successful_validations": 0, "failed_validations": 0}),
            (200, """<cot>
<thinking>Analysis</thinking>
</cot>

<json>
{
  "data": [["CCO"]],
  "explanation": ["Direct precursor"],
  "confidence_scores": [0.8]
}
</json>""", {"total_iterations": 1, "total_tool_calls": 0, "successful_validations": 0, "failed_validations": 0})
        ]
        
        mock_validity_check.return_value = ([["CCO"]], ["Direct precursor"], [0.8])
        
        # Call pipeline
        pathways, explanations, confidence, tool_stats = agentic_llm_pipeline(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            max_tool_iterations=5
        )
        
        # Should have retried
        assert len(pathways) == 1
        assert mock_agentic_call.call_count == 2
    
    @patch('src.utils.llm.agentic_call_LLM')
    def test_tool_statistics_aggregation(self, mock_agentic_call):
        """Test that tool statistics are properly aggregated"""
        # Multiple runs with different tool stats
        mock_agentic_call.side_effect = [
            (400, "", {"total_iterations": 1, "total_tool_calls": 2, "successful_validations": 1, "failed_validations": 1}),
            (200, """<cot>
<thinking>Analysis</thinking>
</cot>

<json>
{
  "data": [["CCO"]],
  "explanation": ["Test"],
  "confidence_scores": [0.8]
}
</json>""", {"total_iterations": 2, "total_tool_calls": 3, "successful_validations": 3, "failed_validations": 0})
        ]
        
        with patch('src.utils.llm.validity_check') as mock_validity:
            mock_validity.return_value = ([["CCO"]], ["Test"], [0.8])
            
            pathways, explanations, confidence, tool_stats = agentic_llm_pipeline(
                molecule="CCCO",
                LLM="claude-opus-4-20250514",
                max_tool_iterations=5
            )
            
            # Check aggregation
            assert tool_stats["total_tool_calls"] == 5  # 2 + 3
            assert tool_stats["total_iterations"] == 3  # 1 + 2
            assert tool_stats["successful_validations"] == 4  # 1 + 3
            assert tool_stats["failed_validations"] == 1  # 1 + 0


class TestLangfuseIntegration:
    """Test cases for Langfuse logging integration"""
    
    @patch('src.utils.llm.langfuse_client')
    @patch('src.utils.llm.LANGFUSE_AVAILABLE', True)
    @patch('src.utils.llm.completion')
    def test_langfuse_trace_creation(self, mock_completion, mock_langfuse):
        """Test that Langfuse trace is created"""
        # Mock Langfuse client
        mock_trace = Mock()
        mock_langfuse.trace.return_value = mock_trace
        
        # Mock LLM response
        mock_response = Mock()
        mock_response.choices = [Mock()]
        mock_response.choices[0].message = Mock()
        mock_response.choices[0].message.content = "<json>...</json>"
        mock_response.choices[0].finish_reason = "end_turn"
        mock_response.choices[0].message.tool_calls = None
        mock_response.usage = Mock(input_tokens=100, output_tokens=50, total_tokens=150)
        mock_completion.return_value = mock_response
        
        # Call function
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule="CCCO",
            LLM="claude-opus-4-20250514",
            trace_id="test_trace_123"
        )
        
        # Check that trace was created
        mock_langfuse.trace.assert_called_once()
        call_args = mock_langfuse.trace.call_args
        assert call_args[1]["id"] == "test_trace_123"
        assert call_args[1]["name"] == "agentic_retrosynthesis"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

