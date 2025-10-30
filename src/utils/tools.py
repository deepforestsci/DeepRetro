"""Tool definitions and execution handlers for agentic LLM pipeline"""
import os
from typing import Dict, Any, Callable, Optional
from rdkit import Chem
from dotenv import load_dotenv

from src.utils.utils_molecule import is_valid_smiles
from src.utils.job_context import logger as context_logger

load_dotenv()

ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING", "true").lower() == "false" else True


def log_message(message: str, logger=None):
    """Log the message
    
    Parameters
    ----------
    message : str
        The message to be logged
    logger : _type_, optional
        The logger object, by default None
    """
    if logger is not None:
        logger.info(message)
    else:
        print(message)


# ============================================================================
# Tool Definitions (Anthropic Format)
# ============================================================================

SMILES_VALIDATOR_TOOL = {
    "name": "validate_smiles",
    "description": (
        "Validates a SMILES string to check if it represents a valid chemical structure. "
        "Returns the canonical SMILES form if valid, or an error message if invalid. "
        "Use this tool before proposing SMILES strings in your retrosynthesis analysis "
        "to ensure all molecules are chemically valid."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "smiles": {
                "type": "string",
                "description": "The SMILES string to validate"
            },
            "context": {
                "type": "string",
                "description": "Optional context about this SMILES (e.g., 'reactant 1', 'product')"
            }
        },
        "required": ["smiles"]
    }
}


# ============================================================================
# Tool Execution Functions
# ============================================================================

def validate_smiles_tool(smiles: str, context: Optional[str] = None) -> Dict[str, Any]:
    """Execute SMILES validation tool
    
    Parameters
    ----------
    smiles : str
        SMILES string to validate
    context : Optional[str], optional
        Context information about the SMILES, by default None
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - valid: bool - Whether the SMILES is valid
        - canonical_smiles: str - Canonical form (if valid)
        - error: str - Error message (if invalid)
        - context: str - Echo back the context
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    
    result = {
        "valid": False,
        "canonical_smiles": None,
        "error": None,
        "context": context or "unknown"
    }
    
    try:
        # Check if valid using existing function
        if not is_valid_smiles(smiles):
            result["error"] = "Invalid SMILES string - cannot be parsed by RDKit"
            log_message(f"SMILES validation failed for '{smiles}': Invalid format", logger)
            return result
        
        # Get canonical form
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Failed to convert SMILES to molecule object"
            log_message(f"SMILES validation failed for '{smiles}': Mol conversion failed", logger)
            return result
        
        canonical = Chem.MolToSmiles(mol, canonical=True)
        result["valid"] = True
        result["canonical_smiles"] = canonical
        
        log_message(f"SMILES validation successful: '{smiles}' -> '{canonical}'", logger)
        
    except Exception as e:
        result["error"] = f"Exception during validation: {str(e)}"
        log_message(f"SMILES validation exception for '{smiles}': {e}", logger)
    
    return result


# ============================================================================
# Tool Registry
# ============================================================================

# Map tool names to their definitions
TOOL_DEFINITIONS = {
    "validate_smiles": SMILES_VALIDATOR_TOOL,
}

# Map tool names to their execution functions
TOOL_EXECUTORS: Dict[str, Callable] = {
    "validate_smiles": validate_smiles_tool,
}


# ============================================================================
# Tool Execution Dispatcher
# ============================================================================

def execute_tool(tool_name: str, tool_input: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a tool by name with given input
    
    Parameters
    ----------
    tool_name : str
        Name of the tool to execute
    tool_input : Dict[str, Any]
        Input parameters for the tool
    
    Returns
    -------
    Dict[str, Any]
        Result from tool execution, always includes 'success' key
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    
    # Check if tool exists
    if tool_name not in TOOL_EXECUTORS:
        error_msg = f"Unknown tool: {tool_name}"
        log_message(error_msg, logger)
        return {
            "success": False,
            "error": error_msg,
            "available_tools": list(TOOL_EXECUTORS.keys())
        }
    
    try:
        # Execute the tool
        log_message(f"Executing tool '{tool_name}' with input: {tool_input}", logger)
        executor = TOOL_EXECUTORS[tool_name]
        result = executor(**tool_input)
        
        # Add success flag if not present
        if "success" not in result:
            result["success"] = result.get("valid", True)
        
        log_message(f"Tool '{tool_name}' execution completed: {result}", logger)
        return result
        
    except Exception as e:
        error_msg = f"Tool execution failed: {str(e)}"
        log_message(f"Tool '{tool_name}' execution error: {e}", logger)
        return {
            "success": False,
            "error": error_msg,
            "exception_type": type(e).__name__
        }


def get_available_tools() -> list:
    """Get list of all available tool definitions for LLM
    
    Returns
    -------
    list
        List of tool definition dictionaries in Anthropic format
    """
    return list(TOOL_DEFINITIONS.values())


def get_tool_names() -> list:
    """Get list of available tool names
    
    Returns
    -------
    list
        List of tool name strings
    """
    return list(TOOL_DEFINITIONS.keys())

