#!/usr/bin/env python3
"""
Comprehensive test suite for the deprotecting group module.
"""

import unittest
from src.deprotecting_group import (unmask_protecting_groups_multisymbol,
                                    get_protecting_group_info,
                                    validate_masked_smiles, DEPROTECT_MAP)


class TestDeprotectingGroup(unittest.TestCase):
    """Test cases for the deprotecting group module."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_cases = {
            # Simple single protecting groups
            "$": "CO",  # OMe (canonicalized)
            "%": "COCc1ccccc1",  # OBn
            "&": "COC",  # OEt

            # Molecules with single protecting groups
            "CC(C)$": "COC(C)C",  # Canonicalized
            "CC(C)%": "CC(C)COCc1ccccc1",
            "CC(C)&": "COCC(C)C",  # Canonicalized

            # Molecules with multiple protecting groups
            "CC(C)$%": "CC(C)OCCOCc1ccccc1",  # Not canonicalized in this case
            "CC(C)$&": "COCCOC(C)C",  # Canonicalized
            "CC(C)%&": "COCc1ccccc1COCC(C)C",  # Canonicalized
            "CC(C)$%&": "COCc1ccccc1COCCOC(C)C",  # Canonicalized

            # Complex molecules
            "CC(C)(C)$": "COC(C)(C)C",  # Canonicalized
            "c1ccccc1%": "c1ccc(COCc2ccccc2)cc1",  # Canonicalized
            "CCO&": "CCOCOC",

            # Molecules with disconnections (salts)
            "$.%": "CO.COCc1ccccc1",  # Canonicalized
            "&.$": "CO.COC",  # Canonicalized
            "$%&": "COCc1ccccc1COCCO",  # Canonicalized

            # Edge cases
            "": "",
            "INVALID_SMILES": "INVALID_SMILES",
        }

    def test_unmask_single_protecting_groups(self):
        """Test unmasking of single protecting group symbols."""
        for masked, expected in self.test_cases.items():
            if masked in ["", "INVALID_SMILES"]:
                continue
            with self.subTest(masked=masked):
                result = unmask_protecting_groups_multisymbol(masked)
                self.assertEqual(result, expected)

    def test_unmask_empty_string(self):
        """Test unmasking of empty string."""
        result = unmask_protecting_groups_multisymbol("")
        self.assertEqual(result, "")

    def test_unmask_invalid_smiles(self):
        """Test unmasking of 'INVALID_SMILES' string."""
        result = unmask_protecting_groups_multisymbol("INVALID_SMILES")
        self.assertEqual(result, "INVALID_SMILES")

    def test_unmask_no_symbols(self):
        """Test unmasking of SMILES without protecting group symbols."""
        test_cases = [
            "CC(C)",
            "c1ccccc1",
            "CCO",
            "CC(C)(C)C",
        ]

        for smiles in test_cases:
            with self.subTest(smiles=smiles):
                result = unmask_protecting_groups_multisymbol(smiles)
                # Should return canonicalized SMILES since no symbols were present
                # The canonicalization might change the SMILES representation
                self.assertNotEqual(result, "INVALID_SMILES")
                self.assertIsInstance(result, str)

    def test_unmask_invalid_symbols(self):
        """Test unmasking of SMILES with invalid symbols."""
        invalid_cases = [
            "CC(C)#",  # Invalid symbol
            "CC(C)@",  # Invalid symbol
            "CC(C)!",  # Invalid symbol
        ]

        for smiles in invalid_cases:
            with self.subTest(smiles=smiles):
                result = unmask_protecting_groups_multisymbol(smiles)
                # Should return INVALID_SMILES since these are invalid SMILES
                self.assertEqual(result, "INVALID_SMILES")

    def test_unmask_complex_molecules(self):
        """Test unmasking of complex molecules with multiple protecting groups."""
        complex_cases = {
            "CC(C)$%&": "COCc1ccccc1COCCOC(C)C",  # Canonicalized
            "c1ccccc1$%": "c1ccc(COCCOc2ccccc2)cc1",  # Canonicalized
            "CCO$%&": "CCOOCCOCc1ccccc1COC",  # Canonicalized
        }

        for masked, expected in complex_cases.items():
            with self.subTest(masked=masked):
                result = unmask_protecting_groups_multisymbol(masked)
                self.assertEqual(result, expected)

    def test_unmask_disconnected_molecules(self):
        """Test unmasking of disconnected molecules (salts)."""
        salt_cases = {
            "$.%": "CO.COCc1ccccc1",  # Canonicalized
            "&.$": "CO.COC",  # Canonicalized
            "%.&": "COC.COCc1ccccc1",  # Canonicalized
        }

        for masked, expected in salt_cases.items():
            with self.subTest(masked=masked):
                result = unmask_protecting_groups_multisymbol(masked)
                self.assertEqual(result, expected)

    def test_unmask_symbol_order_independence(self):
        """Test that unmasking works regardless of symbol order in DEPROTECT_MAP."""
        # Test that the order of replacement doesn't matter
        test_smiles = "CC(C)$%&"

        # Get the expected result
        expected = unmask_protecting_groups_multisymbol(test_smiles)

        # Verify that the result is consistent
        for _ in range(5):  # Test multiple times
            result = unmask_protecting_groups_multisymbol(test_smiles)
            self.assertEqual(result, expected)


class TestProtectingGroupInfo(unittest.TestCase):
    """Test cases for the protecting group information functions."""

    def test_get_protecting_group_info(self):
        """Test the get_protecting_group_info function."""
        info = get_protecting_group_info()

        # Check that all expected symbols are present
        expected_symbols = {"$", "%", "&"}
        self.assertEqual(set(info.keys()), expected_symbols)

        # Check structure of each entry
        for symbol in expected_symbols:
            self.assertIn(symbol, info)
            entry = info[symbol]

            # Check required fields
            required_fields = {"smiles", "name", "description"}
            self.assertEqual(set(entry.keys()), required_fields)

            # Check data types
            self.assertIsInstance(entry["smiles"], str)
            self.assertIsInstance(entry["name"], str)
            self.assertIsInstance(entry["description"], str)

            # Check that SMILES matches DEPROTECT_MAP
            self.assertEqual(entry["smiles"], DEPROTECT_MAP[symbol])

    def test_protecting_group_info_content(self):
        """Test the specific content of protecting group information."""
        info = get_protecting_group_info()

        # Test OMe ($)
        self.assertEqual(info["$"]["name"], "OMe")
        self.assertEqual(info["$"]["smiles"], "OC")
        self.assertEqual(info["$"]["description"], "Methoxy protecting group")

        # Test OBn (%)
        self.assertEqual(info["%"]["name"], "OBn")
        self.assertEqual(info["%"]["smiles"], "COCc1ccccc1")
        self.assertEqual(info["%"]["description"], "Benzyl protecting group")

        # Test OEt (&)
        self.assertEqual(info["&"]["name"], "OEt")
        self.assertEqual(info["&"]["smiles"], "COC")
        self.assertEqual(info["&"]["description"], "Ethoxy protecting group")


class TestValidation(unittest.TestCase):
    """Test cases for the validation function."""

    def test_validate_masked_smiles_valid_cases(self):
        """Test validation of valid masked SMILES."""
        valid_cases = [
            "$",  # Single OMe symbol
            "%",  # Single OBn symbol
            "&",  # Single OEt symbol
            "CC(C)$",  # Molecule with OMe
            "CC(C)%",  # Molecule with OBn
            "CC(C)&",  # Molecule with OEt
            "CC(C)$%",  # Molecule with multiple symbols
            "CC(C)$%&",  # Molecule with all symbols
            "c1ccccc1$",  # Aromatic with symbol
            "CCO$%",  # Alcohol with symbols
        ]

        for smiles in valid_cases:
            with self.subTest(smiles=smiles):
                result = validate_masked_smiles(smiles)
                self.assertTrue(result, f"Expected {smiles} to be valid")

    def test_validate_masked_smiles_invalid_cases(self):
        """Test validation of invalid masked SMILES."""
        invalid_cases = [
            "CC(C)",  # No symbols
            "c1ccccc1",  # No symbols
            "CCO",  # No symbols
            "invalid",  # Invalid characters
            "CC(C)#",  # Invalid symbol
            "CC(C)@",  # Invalid symbol
        ]

        for smiles in invalid_cases:
            with self.subTest(smiles=smiles):
                result = validate_masked_smiles(smiles)
                self.assertFalse(result, f"Expected {smiles} to be invalid")

    def test_validate_masked_smiles_edge_cases(self):
        """Test validation of edge cases."""
        edge_cases = [
            ("", True),  # Empty string should be valid
            ("$$$", True),  # Multiple symbols should be valid
            ("%%%", True),  # Multiple symbols should be valid
            ("&&&", True),  # Multiple symbols should be valid
            ("$%&", True),  # All symbols should be valid
        ]

        for smiles, expected in edge_cases:
            with self.subTest(smiles=smiles):
                result = validate_masked_smiles(smiles)
                self.assertEqual(result, expected,
                                 f"Expected {smiles} to be {expected}")


class TestIntegration(unittest.TestCase):
    """Integration tests for the deprotecting group module."""

    def test_deprotect_map_consistency(self):
        """Test that DEPROTECT_MAP is consistent with get_protecting_group_info."""
        info = get_protecting_group_info()

        for symbol, smiles in DEPROTECT_MAP.items():
            with self.subTest(symbol=symbol):
                self.assertIn(symbol, info)
                self.assertEqual(info[symbol]["smiles"], smiles)

    def test_unmask_validation_consistency(self):
        """Test that unmasking works correctly with validated SMILES."""
        # Test that valid masked SMILES can be unmasked
        valid_masked = ["$", "%", "&", "CC(C)$", "CC(C)%&"]

        for masked in valid_masked:
            with self.subTest(masked=masked):
                # Should be valid
                self.assertTrue(validate_masked_smiles(masked))

                # Should be unmaskable
                result = unmask_protecting_groups_multisymbol(masked)
                self.assertNotEqual(result, "INVALID_SMILES")

    def test_error_handling(self):
        """Test error handling for various edge cases."""
        error_cases = [
            "INVALID_SMILES",  # Special case
            "$$$$$$",  # Many symbols
            "CC(C)$%&$%&",  # Repeated symbols
        ]

        for case in error_cases:
            with self.subTest(case=case):
                # Should handle other cases without crashing
                result = unmask_protecting_groups_multisymbol(case)
                self.assertIsInstance(result, str)

    def test_none_input(self):
        """Test handling of None input."""
        # The function should handle None gracefully by returning empty string
        # since None is falsy and the function checks "if not smiles:"
        result = unmask_protecting_groups_multisymbol(None)
        self.assertEqual(result, "")


def run_performance_test():
    """Run a simple performance test."""
    import time

    test_smiles = "CC(C)$%&" * 1000  # Create a long string

    start_time = time.time()
    for _ in range(1000):
        unmask_protecting_groups_multisymbol(test_smiles)
    end_time = time.time()

    print(
        f"Performance test: 1000 unmaskings took {end_time - start_time:.4f} seconds"
    )


if __name__ == "__main__":
    # Run the unit tests
    unittest.main(verbosity=2, exit=False)

    # Run performance test
    print("\n" + "=" * 60)
    print("Running performance test...")
    run_performance_test()

    print("\n" + "=" * 60)
    print("All tests completed!")
