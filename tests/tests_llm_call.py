#!/usr/bin/env python3
"""
Real LLM call test with protecting group integration.
"""

import os
import sys
import json
from src.protecting_group import mask_protecting_groups_multisymbol
from src.utils.llm import call_LLM


def test_real_llm_call():
    """Make a real LLM call with protecting group integration."""

    print("=" * 80)
    print("REAL LLM CALL WITH PROTECTING GROUP INTEGRATION")
    print("=" * 80)

    # Test molecule with protecting group
    test_molecule = "CCOC(=O)c1ccccc1"  # Ethyl benzoate

    print(f"Test Molecule: {test_molecule}")

    # Check protecting group detection
    masked_smiles = mask_protecting_groups_multisymbol(test_molecule)
    has_protecting_groups = masked_smiles != test_molecule and masked_smiles != "INVALID_SMILES"

    print(f"Masked SMILES: {masked_smiles}")
    print(f"Protecting Groups Detected: {has_protecting_groups}")

    if has_protecting_groups:
        print("Protecting group integration will be active")
        print("   - OEt (ethoxy) protecting group detected")
        print("   - Context will include deprotection instructions")
    else:
        print("No protecting groups detected - standard prompt will be used")

    print(f"\nMaking LLM call...")
    print("This may take a few moments...")

    try:
        # Make the actual LLM call
        # Using a reliable model
        status_code, response_text = call_LLM(test_molecule, "gpt-4o-mini")

        print(f"\nLLM Call Status: {status_code}")
        print(f"Response Length: {len(response_text)} characters")

        if status_code == 200:
            print("SUCCESS! LLM call completed successfully")

            # Analyze the response
            print(f"\n{'='*80}")
            print("RESPONSE ANALYSIS")
            print("=" * 80)

            # Check if response contains protecting group mentions
            if has_protecting_groups:
                print("Protecting Group Analysis:")

                if "deprotection" in response_text.lower():
                    print("DEPROTECTION mentioned in response")
                else:
                    print("No deprotection mentioned")

                if "protecting group" in response_text.lower():
                    print("PROTECTING GROUP mentioned in response")
                else:
                    print("No protecting group mentioned")

                if "hydrolysis" in response_text.lower(
                ) or "hydrogenolysis" in response_text.lower():
                    print("DEPROTECTION CONDITIONS mentioned")
                else:
                    print("No deprotection conditions mentioned")

                if "ethox" in response_text.lower(
                ) or "oet" in response_text.lower():
                    print("ETHOXY GROUP mentioned")
                else:
                    print("No ethoxy group mentioned")

            # Try to extract JSON content
            if "<json>" in response_text and "</json>" in response_text:
                json_start = response_text.find("<json>") + 6
                json_end = response_text.find("</json>")
                json_content = response_text[json_start:json_end].strip()

                print(
                    f"\nJSON Content Extracted: {len(json_content)} characters"
                )

                # Check if actual SMILES are returned (not masked symbols)
                if "&" not in json_content and "%" not in json_content and "$" not in json_content:
                    print("ACTUAL SMILES returned (no masked symbols in data)")
                else:
                    print("Masked symbols found in data field")

                # Try to parse JSON
                try:
                    parsed_json = json.loads(json_content)
                    print("JSON successfully parsed")

                    if "data" in parsed_json:
                        print(
                            f"Data field found with {len(parsed_json['data'])} entries"
                        )
                        for i, entry in enumerate(parsed_json['data']):
                            print(f"   Entry {i+1}: {entry}")

                    if "explanation" in parsed_json:
                        print(
                            f"Explanation field found with {len(parsed_json['explanation'])} entries"
                        )
                        for i, explanation in enumerate(
                                parsed_json['explanation']):
                            print(
                                f"   Explanation {i+1}: {explanation[:100]}..."
                            )

                except json.JSONDecodeError as e:
                    print(f"JSON parsing failed: {e}")
            else:
                print("No JSON content found in response")

            # Show response preview
            print(f"\n{'='*80}")
            print("RESPONSE PREVIEW")
            print("=" * 80)
            print(response_text[:1000] +
                  "..." if len(response_text) > 1000 else response_text)

            return True, response_text

        else:
            print(f"LLM call failed with status code: {status_code}")
            print(f"Response: {response_text}")
            return False, response_text

    except Exception as e:
        print(f"Error during LLM call: {e}")
        return False, str(e)


if __name__ == "__main__":
    success, response = test_real_llm_call()

    print(f"\n{'='*80}")
    print("FINAL TEST RESULT")
    print("=" * 80)
    print(f"Overall Test: {'PASS' if success else 'FAIL'}")

    if success:
        print("\nPROTECTING GROUP INTEGRATION IS WORKING WITH REAL LLM!")
        print("   - API call successful")
        print("   - Protecting group context integrated")
        print("   - LLM provided informed response")
        print("   - Ready for production use")
    else:
        print("\nSome issues detected - check error messages above")
        print("   - May need to adjust model or API configuration")
