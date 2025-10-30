"""Unit tests for tools module"""
import pytest
from src.utils.tools import (
    validate_smiles_tool,
    execute_tool,
    get_available_tools,
    get_tool_names,
    SMILES_VALIDATOR_TOOL
)


class TestSMILESValidatorTool:
    """Test cases for SMILES validator tool"""
    
    def test_valid_smiles(self):
        """Test validation of valid SMILES strings"""
        result = validate_smiles_tool("CCO")
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
        assert result["error"] is None
    
    def test_invalid_smiles(self):
        """Test validation of invalid SMILES strings"""
        result = validate_smiles_tool("InvalidSMILES@#$")
        assert result["valid"] is False
        assert result["canonical_smiles"] is None
        assert result["error"] is not None
    
    def test_benzene_canonicalization(self):
        """Test canonical form of benzene"""
        result = validate_smiles_tool("c1ccccc1")
        assert result["valid"] is True
        # Check that canonical form is provided
        assert result["canonical_smiles"] in ["c1ccccc1", "C1=CC=CC=C1"]
    
    def test_with_context(self):
        """Test validation with context parameter"""
        result = validate_smiles_tool("CCO", context="reactant 1")
        assert result["valid"] is True
        assert result["context"] == "reactant 1"
    
    def test_complex_molecule(self):
        """Test validation of complex molecule"""
        # Aspirin
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        result = validate_smiles_tool(smiles)
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
    
    def test_stereochemistry(self):
        """Test validation preserves stereochemistry"""
        # L-alanine
        smiles = "C[C@H](N)C(=O)O"
        result = validate_smiles_tool(smiles)
        assert result["valid"] is True
        assert "@" in result["canonical_smiles"] or "[C@@H]" in result["canonical_smiles"]


class TestToolExecutor:
    """Test cases for tool execution dispatcher"""
    
    def test_execute_validate_smiles(self):
        """Test executing validate_smiles tool"""
        result = execute_tool("validate_smiles", {"smiles": "CCO"})
        assert result["success"] is True
        assert result["valid"] is True
    
    def test_unknown_tool(self):
        """Test executing unknown tool"""
        result = execute_tool("unknown_tool", {})
        assert result["success"] is False
        assert "Unknown tool" in result["error"]
        assert "available_tools" in result
    
    def test_tool_error_handling(self):
        """Test tool execution error handling"""
        # Missing required parameter
        result = execute_tool("validate_smiles", {})
        assert result["success"] is False
        assert "error" in result
    
    def test_invalid_smiles_execution(self):
        """Test executing tool with invalid SMILES"""
        result = execute_tool("validate_smiles", {"smiles": "InvalidSMILES@#$"})
        # Should return result but indicate invalid
        assert "valid" in result
        assert result["valid"] is False


class TestToolRegistry:
    """Test cases for tool registry functions"""
    
    def test_get_available_tools(self):
        """Test getting available tools"""
        tools = get_available_tools()
        assert isinstance(tools, list)
        assert len(tools) > 0
        # Check format
        for tool in tools:
            assert "name" in tool
            assert "description" in tool
            assert "input_schema" in tool
    
    def test_get_tool_names(self):
        """Test getting tool names"""
        names = get_tool_names()
        assert isinstance(names, list)
        assert "validate_smiles" in names
    
    def test_tool_definition_format(self):
        """Test that tool definitions follow Anthropic format"""
        assert "name" in SMILES_VALIDATOR_TOOL
        assert "description" in SMILES_VALIDATOR_TOOL
        assert "input_schema" in SMILES_VALIDATOR_TOOL
        
        schema = SMILES_VALIDATOR_TOOL["input_schema"]
        assert schema["type"] == "object"
        assert "properties" in schema
        assert "smiles" in schema["properties"]
        assert "required" in schema


class TestToolIntegration:
    """Integration tests for tools"""
    
    def test_validate_multiple_smiles(self):
        """Test validating multiple SMILES in sequence"""
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "InvalidSMILES"]
        results = [validate_smiles_tool(s) for s in smiles_list]
        
        # First three should be valid
        assert results[0]["valid"] is True
        assert results[1]["valid"] is True
        assert results[2]["valid"] is True
        
        # Last one should be invalid
        assert results[3]["valid"] is False
    
    def test_canonical_forms_consistency(self):
        """Test that canonical forms are consistent"""
        # Different representations of same molecule
        smiles1 = "c1ccccc1"  # benzene
        smiles2 = "C1=CC=CC=C1"  # benzene (different representation)
        
        result1 = validate_smiles_tool(smiles1)
        result2 = validate_smiles_tool(smiles2)
        
        assert result1["valid"] is True
        assert result2["valid"] is True
        # Canonical forms should be the same
        assert result1["canonical_smiles"] == result2["canonical_smiles"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

