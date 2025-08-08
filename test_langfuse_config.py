#!/usr/bin/env python3
"""Test script for Langfuse configuration helper."""

import sys
import os
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from utils.langfuse_config import get_langfuse_metadata


def test_langfuse_config():
    """Test the Langfuse configuration helper."""
    print("Testing Langfuse configuration helper...")
    
    # Test default session
    print("\n1. Testing default session:")
    metadata = get_langfuse_metadata()
    print(f"   Default metadata: {metadata}")
    assert metadata["session_id"] == "default"
    assert metadata["generation_name"] == "retrosynthesis"
    
    # Test metadata session
    print("\n2. Testing metadata session:")
    metadata = get_langfuse_metadata("metadata")
    print(f"   Metadata session: {metadata}")
    assert metadata["session_id"] == "metadata_agent"
    assert metadata["generation_name"] == "metadata_prediction"
    assert metadata["trace_name"] == "metadata_pipeline"
    
    # Test retrosynthesis session
    print("\n3. Testing retrosynthesis session:")
    metadata = get_langfuse_metadata("retrosynthesis")
    print(f"   Retrosynthesis session: {metadata}")
    assert metadata["session_id"] == "retrosynthesis_pipeline"
    assert metadata["generation_name"] == "retrosynthesis"
    assert metadata["trace_name"] == "retrosynthesis_pipeline"
    
    # Test unknown session (should raise KeyError)
    print("\n4. Testing unknown session:")
    try:
        metadata = get_langfuse_metadata("unknown")
        print("Expected KeyError but got success")
        assert False, "Should have raised KeyError"
    except KeyError as e:
        print(f"Correctly raised KeyError: {e}")
    
    print("\n Langfuse configuration helper test completed!")


if __name__ == "__main__":
    test_langfuse_config() 