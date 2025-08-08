"""Simple Langfuse configuration helper."""

import json
import os
from pathlib import Path
from typing import Dict, Any


def get_langfuse_metadata(session_type: str = "default") -> Dict[str, Any]:
    """Get Langfuse metadata from config file.
    
    Args:
        session_type: Type of session ("default", "metadata", "retrosynthesis")
        
    Returns:
        Metadata dictionary for Langfuse
        
    Raises:
        FileNotFoundError: If config file doesn't exist
        KeyError: If session_type doesn't exist in config
        json.JSONDecodeError: If config file is invalid JSON
    """
    # Try to load config file
    config_path = Path(__file__).parent.parent.parent / "config" / "langfuse_config.json"
    
    if not config_path.exists():
        raise FileNotFoundError(f"Langfuse config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Get the session config
    sessions = config.get('sessions', {})
    
    if session_type in sessions:
        # Return the full config for this session
        return sessions[session_type].copy()
    else:
        # Session type not found
        raise KeyError(f"Session type '{session_type}' not found in config. Available sessions: {list(sessions.keys())}") 