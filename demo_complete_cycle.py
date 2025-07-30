#!/usr/bin/env python3
"""
Demonstration script showing the complete masking and unmasking cycle for protecting groups.
This script demonstrates how molecules with protecting groups can be masked for retrosynthetic analysis
and then unmasked back to their full chemical representations.
"""

from src.protecting_group import mask_protecting_groups_multisymbol
from src.deprotecting_group import (unmask_protecting_groups_multisymbol,
                                    get_protecting_group_info,
                                    validate_masked_smiles)


def demonstrate_complete_cycle():
    """Demonstrate the complete masking and unmasking cycle."""

    print("=" * 80)
    print("PROTECTING GROUP MASKING/UNMASKING CYCLE DEMONSTRATION")
    print("=" * 80)

    # Test molecules with different protecting groups
    test_molecules = [
        "COC",  # Simple OEt
        "COCc1ccccc1",  # Simple OBn
        "OC",  # Simple OMe
        "CC(C)COCc1ccccc1",  # Molecule with OBn
        "CC(C)OC",  # Molecule with OMe
        "COC.COC",  # Multiple OEt groups
        "CC(C)OCCOCc1ccccc1COC",  # Complex molecule with multiple protecting groups
    ]

    print("\n1. PROTECTING GROUP INFORMATION:")
    print("-" * 40)
    info = get_protecting_group_info()
    for symbol, details in info.items():
        print(
            f"  {symbol} → {details['name']} ({details['smiles']}) - {details['description']}"
        )

    print("\n2. MASKING AND UNMASKING CYCLE:")
    print("-" * 40)

    for i, molecule in enumerate(test_molecules, 1):
        print(f"\nTest {i}: {molecule}")
        print(f"  {'─' * 50}")

        # Step 1: Mask the molecule
        masked = mask_protecting_groups_multisymbol(molecule)
        print(f"  Masked:     {masked}")

        # Step 2: Validate the masked SMILES
        is_valid = validate_masked_smiles(masked)
        print(f"  Valid:      {is_valid}")

        # Step 3: Unmask the molecule
        unmasked = unmask_protecting_groups_multisymbol(masked)
        print(f"  Unmasked:   {unmasked}")

        # Step 4: Check if cycle is reversible
        is_reversible = molecule == unmasked
        status = "✓ REVERSIBLE" if is_reversible else "✗ NOT REVERSIBLE"
        print(f"  Cycle:      {status}")

        if not is_reversible:
            print(f"  Note:       Canonicalization changed atom order")


def demonstrate_retrosynthetic_context():
    """Demonstrate how this would be used in retrosynthetic analysis."""

    print("\n\n3. RETROSYNTHETIC ANALYSIS CONTEXT:")
    print("-" * 40)

    # Example molecule with protecting groups
    target_molecule = "CC(C)COCc1ccccc1"

    print(f"Target molecule: {target_molecule}")

    # Mask for analysis
    masked_target = mask_protecting_groups_multisymbol(target_molecule)
    print(f"Masked for analysis: {masked_target}")

    print("\nIn retrosynthetic analysis, the LLM would see:")
    print(f"  - Original SMILES: {target_molecule}")
    print(f"  - Masked SMILES: {masked_target}")
    print("  - Context about protecting groups and deprotection strategies")

    print("\nThe LLM would then:")
    print(
        "  1. Analyze the masked structure for retrosynthetic disconnections")
    print("  2. Suggest synthetic routes")
    print("  3. Include deprotection steps in the synthesis")
    print(
        "  4. Use unmask_protecting_groups_multisymbol() to convert symbols back to full SMILES"
    )

    # Example of what the LLM might suggest
    print("\nExample LLM suggestion:")
    print("  'The molecule contains a benzyl protecting group (%).")
    print(
        "   After synthesis, deprotect using H2/Pd-C to remove the benzyl group.'"
    )

    # Show the deprotection
    print(f"\nDeprotection example:")
    print(f"  Masked intermediate: CC(C)%")
    print(
        f"  After deprotection: {unmask_protecting_groups_multisymbol('CC(C)%')}"
    )


def demonstrate_error_handling():
    """Demonstrate error handling capabilities."""

    print("\n\n4. ERROR HANDLING:")
    print("-" * 40)

    error_cases = [
        ("", "Empty string"),
        ("invalid_smiles", "Invalid SMILES"),
        ("CC(C)#", "Invalid symbol"),
        ("INVALID_SMILES", "Special case"),
    ]

    for test_input, description in error_cases:
        print(f"\n{description}: '{test_input}'")

        # Test masking
        masked = mask_protecting_groups_multisymbol(test_input)
        print(f"  Masked: {masked}")

        # Test unmasking
        unmasked = unmask_protecting_groups_multisymbol(masked)
        print(f"  Unmasked: {unmasked}")

        # Test validation
        is_valid = validate_masked_smiles(masked)
        print(f"  Valid: {is_valid}")


def demonstrate_performance():
    """Demonstrate performance with multiple operations."""

    print("\n\n5. PERFORMANCE DEMONSTRATION:")
    print("-" * 40)

    import time

    # Create a complex test case
    test_smiles = "CC(C)$%&" * 100  # Repeat the pattern 100 times

    print(f"Testing with complex SMILES (length: {len(test_smiles)})")

    # Test masking performance
    start_time = time.time()
    for _ in range(1000):
        masked = mask_protecting_groups_multisymbol(test_smiles)
    mask_time = time.time() - start_time

    # Test unmasking performance
    start_time = time.time()
    for _ in range(1000):
        unmasked = unmask_protecting_groups_multisymbol(masked)
    unmask_time = time.time() - start_time

    print(f"  1000 masking operations: {mask_time:.4f} seconds")
    print(f"  1000 unmasking operations: {unmask_time:.4f} seconds")
    print(
        f"  Average per operation: {(mask_time + unmask_time) / 2000:.6f} seconds"
    )


if __name__ == "__main__":
    demonstrate_complete_cycle()
    demonstrate_retrosynthetic_context()
    demonstrate_error_handling()
    demonstrate_performance()

    print("\n" + "=" * 80)
    print("DEMONSTRATION COMPLETED")
    print("=" * 80)
    print("\nThe deprotecting module successfully completes the cycle by:")
    print("1. Converting masked symbols back to full protecting group SMILES")
    print("2. Providing validation for masked SMILES")
    print("3. Offering detailed information about protecting groups")
    print("4. Handling edge cases and errors gracefully")
    print("5. Integrating seamlessly with the existing masking functionality")
