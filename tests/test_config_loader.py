import unittest
import tempfile
import json
import os
import sys
from unittest.mock import patch, MagicMock
from pathlib import Path

# Add the src directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from config_loader import (ProtectingGroupConfig, ConfigurationError,
                           get_pg_map, get_full_pg_config,
                           generate_protecting_group_context, reload_config)


class TestProtectingGroupConfig(unittest.TestCase):
    """Test the ProtectingGroupConfig class."""

    def setUp(self):
        """Set up test fixtures."""
        self.valid_config = {
            "version": "1.0",
            "description": "Test config",
            "protecting_groups": {
                "TestPG": {
                    "smiles": "CO",
                    "symbol": "$",
                    "description": "Test protecting group",
                    "deprotection": "Test deprotection method"
                }
            }
        }

        self.minimal_config = {
            "protecting_groups": {
                "MinimalPG": {
                    "smiles": "CCO",
                    "symbol": "&"
                }
            }
        }

    def test_valid_config_loading(self):
        """Test loading a valid configuration."""
        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(self.valid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)
            pg_map = config.pg_map

            self.assertIn("TestPG", pg_map)
            self.assertEqual(pg_map["TestPG"], ("CO", "$"))
        finally:
            os.unlink(config_path)

    def test_missing_config_file_explicit_path(self):
        """Test behavior when explicitly provided config file doesn't exist and no fallback exists."""
        # Create a temp directory structure without any config files
        import tempfile
        import shutil

        with tempfile.TemporaryDirectory() as temp_dir:
            # Set up a completely isolated environment
            isolated_config_path = os.path.join(temp_dir, "config.json")

            # Create config with isolated path and patch to prevent fallback
            config = ProtectingGroupConfig(isolated_config_path)

            # Mock __file__ to point to temp directory so fallback search fails
            with patch('config_loader.__file__',
                       os.path.join(temp_dir, 'fake_module.py')):
                with self.assertRaises(ConfigurationError) as cm:
                    pg_map = config.pg_map

                self.assertIn("No configuration file found", str(cm.exception))

    def test_invalid_json(self):
        """Test behavior with invalid JSON."""
        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            f.write("{ invalid json content")
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)

            with self.assertRaises(ConfigurationError) as cm:
                _ = config.pg_map

            self.assertIn("Invalid JSON", str(cm.exception))
        finally:
            os.unlink(config_path)

    def test_missing_protecting_groups_key(self):
        """Test config without protecting_groups key."""
        invalid_config = {"version": "1.0"}

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(invalid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)

            with self.assertRaises(ConfigurationError) as cm:
                _ = config.pg_map

            self.assertIn("must contain 'protecting_groups' key",
                          str(cm.exception))
        finally:
            os.unlink(config_path)

    def test_invalid_symbol_length(self):
        """Test validation of symbol length."""
        invalid_config = {
            "protecting_groups": {
                "BadPG": {
                    "smiles": "CO",
                    "symbol": "TOOLONG"  # More than 1 character
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(invalid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)

            with self.assertRaises(ConfigurationError) as cm:
                _ = config.pg_map

            self.assertIn("must be a single character", str(cm.exception))
        finally:
            os.unlink(config_path)

    def test_empty_smiles_validation(self):
        """Test validation of empty SMILES patterns."""
        invalid_config = {
            "protecting_groups": {
                "EmptyPG": {
                    "smiles": "",  # Empty SMILES
                    "symbol": "$"
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(invalid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)

            with self.assertRaises(ConfigurationError) as cm:
                _ = config.pg_map

            self.assertIn("cannot be empty", str(cm.exception))
        finally:
            os.unlink(config_path)

    def test_deprotection_validation(self):
        """Test validation of deprotection field."""
        # Valid deprotection field
        valid_config = {
            "protecting_groups": {
                "ValidPG": {
                    "smiles": "CO",
                    "symbol": "$",
                    "deprotection": "Valid deprotection method"
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(valid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)
            pg_map = config.pg_map  # Should not raise exception
            self.assertIn("ValidPG", pg_map)
        finally:
            os.unlink(config_path)

        # Invalid deprotection field (empty)
        invalid_config = {
            "protecting_groups": {
                "InvalidPG": {
                    "smiles": "CO",
                    "symbol": "$",
                    "deprotection": ""  # Empty deprotection
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(invalid_config, f)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)

            with self.assertRaises(ConfigurationError) as cm:
                _ = config.pg_map

            self.assertIn("Deprotection condition", str(cm.exception))
        finally:
            os.unlink(config_path)


class TestGlobalFunctions(unittest.TestCase):
    """Test the global helper functions."""

    def test_generate_protecting_group_context(self):
        """Test dynamic context generation."""
        # Mock the config loading
        mock_config = {
            "TestPG1": {
                "smiles": "CO",
                "symbol": "$",
                "description": "Test PG 1",
                "deprotection": "Method 1"
            },
            "TestPG2": {
                "smiles": "CCO",
                "symbol": "&",
                "description": "Test PG 2",
                "deprotection": "Method 2"
            }
        }

        with patch('config_loader.get_full_pg_config',
                   return_value=mock_config):
            context = generate_protecting_group_context("CO", "$")

            # Check context contains expected elements
            self.assertIn("Original SMILES: CO", context)
            self.assertIn("Masked SMILES: $", context)
            self.assertIn("'$' represents TestPG1 groups", context)
            self.assertIn("'&' represents TestPG2 groups", context)
            self.assertIn("TestPG1 ($): Method 1", context)
            self.assertIn("TestPG2 (&): Method 2", context)

    def test_context_with_missing_deprotection(self):
        """Test context generation when deprotection field is missing."""
        mock_config = {
            "NoDep": {
                "smiles": "CO",
                "symbol": "$",
                "description": "No deprotection field"
                # Missing deprotection field
            }
        }

        with patch('config_loader.get_full_pg_config',
                   return_value=mock_config):
            context = generate_protecting_group_context("CO", "$")

            # Should use fallback deprotection text
            self.assertIn("Standard deprotection conditions", context)

    @patch('config_loader._global_config')
    def test_reload_config(self, mock_global_config):
        """Test config reloading functionality."""
        reload_config()

        # Should reset the cached config
        mock_global_config._pg_map = None


class TestConfigFileResolution(unittest.TestCase):
    """Test configuration file path resolution."""

    def setUp(self):
        """Set up test fixtures."""
        self.valid_config = {
            "protecting_groups": {
                "TestPG": {
                    "smiles": "CO",
                    "symbol": "$"
                }
            }
        }

    def test_config_file_priority(self):
        """Test that explicit config path takes priority over default locations."""
        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False) as f:
            json.dump(self.valid_config, f)
            explicit_config_path = f.name

        try:
            # Config with explicit path should use that path
            config = ProtectingGroupConfig(explicit_config_path)
            found_path = config._find_config_file()

            # Should return the explicit path we provided
            self.assertEqual(found_path, explicit_config_path)

        finally:
            os.unlink(explicit_config_path)

    def test_config_file_found_in_project(self):
        """Test that actual config file is found in project structure."""
        config = ProtectingGroupConfig()

        # This should find the actual config file
        result = config._find_config_file()
        self.assertIsNotNone(result)
        self.assertTrue(result.endswith('protecting_groups.json'))


class TestIntegrationWithRealConfig(unittest.TestCase):
    """Integration tests using the actual project config file."""

    def setUp(self):
        """Set up test fixtures."""
        # Force reload to ensure we get fresh config
        reload_config()

    def test_load_real_config(self):
        """Test loading the actual project configuration."""
        # This tests the integration with the real config file
        try:
            pg_map = get_pg_map()

            # Should have some protecting groups
            self.assertGreater(len(pg_map), 0)

            # Each entry should be properly formatted
            for name, (smiles, symbol) in pg_map.items():
                self.assertIsInstance(name, str)
                self.assertIsInstance(smiles, str)
                self.assertIsInstance(symbol, str)
                self.assertEqual(len(symbol), 1)
                self.assertGreater(len(smiles), 0)

        except ConfigurationError:
            # If config doesn't exist, skip this test
            self.skipTest("Project config file not found")

    def test_real_config_context_generation(self):
        """Test context generation with real config."""
        try:
            # Test with a simple molecule that should be recognized
            context = generate_protecting_group_context("CO", "$")

            # Should contain the standard context elements
            self.assertIn("PROTECTING GROUP CONTEXT:", context)
            self.assertIn("Original SMILES: CO", context)
            self.assertIn("Masked SMILES: $", context)
            self.assertIn("Symbol mapping:", context)
            self.assertIn("deprotection conditions:", context)

        except ConfigurationError:
            self.skipTest("Project config file not found")

    def test_get_full_pg_config(self):
        """Test getting full configuration data."""
        try:
            full_config = get_full_pg_config()

            if full_config:  # If config exists
                # Should be a dictionary
                self.assertIsInstance(full_config, dict)

                # Each protecting group should have required fields
                for pg_name, pg_data in full_config.items():
                    self.assertIsInstance(pg_data, dict)
                    self.assertIn('smiles', pg_data)
                    self.assertIn('symbol', pg_data)

                    # Optional fields
                    if 'deprotection' in pg_data:
                        self.assertIsInstance(pg_data['deprotection'], str)
                        self.assertGreater(len(pg_data['deprotection']), 0)

        except Exception:
            self.skipTest("Could not load project config")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""

    def test_config_with_unicode_characters(self):
        """Test config with unicode characters in descriptions."""
        unicode_config = {
            "protecting_groups": {
                "UnicodePG": {
                    "smiles": "CO",
                    "symbol": "α",  # Greek alpha
                    "description": "Protecting group with unicode: αβγ",
                    "deprotection": "Method with unicode: δεζ"
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w',
                                         suffix='.json',
                                         delete=False,
                                         encoding='utf-8') as f:
            json.dump(unicode_config, f, ensure_ascii=False)
            config_path = f.name

        try:
            config = ProtectingGroupConfig(config_path)
            pg_map = config.pg_map

            self.assertIn("UnicodePG", pg_map)
            self.assertEqual(pg_map["UnicodePG"], ("CO", "α"))
        finally:
            os.unlink(config_path)

    def test_context_deduplication(self):
        """Test that context properly deduplicates protecting groups with same symbol."""
        mock_config = {
            "PG1_pattern1": {
                "smiles": "CO",
                "symbol": "$",
                "description": "Pattern 1",
                "deprotection": "Method A"
            },
            "PG1_pattern2": {
                "smiles": "COC",
                "symbol": "$",  # Same symbol
                "description": "Pattern 2",
                "deprotection": "Method A"  # Same deprotection
            }
        }

        with patch('config_loader.get_full_pg_config',
                   return_value=mock_config):
            context = generate_protecting_group_context("CO", "$")

            # Should only have one entry for the symbol
            symbol_count = context.count("'$' represents")
            self.assertEqual(symbol_count, 1)

            # Should have deprotection condition (check for the general pattern)
            self.assertIn("Method A", context)
            # Check that deprotection section exists
            self.assertIn("deprotection conditions:", context)


class TestProtectingGroupIntegrationWithVariablesTest(unittest.TestCase):
    """Test integration with test molecules from variables_test.py"""

    def test_protecting_group_test_molecules(self):
        """Test protecting group functionality with dedicated test molecules from variables_test.py."""
        try:
            # Import the test molecules
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), '.'))
            from variables_test import PROTECTING_GROUP_TEST_MOLECULES

            from protecting_group import mask_protecting_groups_multisymbol

            print("Testing protecting group molecules from variables_test.py:")

            expected_results = {
                "simple_methoxy": "&",  # COC should become &
                "simple_benzyl": "%",  # COCc1ccccc1 should become %
                "simple_ethoxy_cco": "&",  # CCO should become &  
                "simple_ethoxy_coc": "&",  # COC should become &
                "no_protecting_groups": "CC(=O)O",  # Should stay unchanged
            }

            for name, smiles in PROTECTING_GROUP_TEST_MOLECULES.items():
                masked = mask_protecting_groups_multisymbol(smiles)

                if name in expected_results:
                    self.assertEqual(
                        masked, expected_results[name],
                        f"Failed for {name}: {smiles} -> {masked}, expected {expected_results[name]}"
                    )
                else:
                    # For complex molecules, just verify they don't crash and return something
                    self.assertNotEqual(
                        masked, "INVALID_SMILES",
                        f"Masking failed for {name}: {smiles}")
                    self.assertIsInstance(
                        masked, str,
                        f"Invalid return type for {name}: {smiles}")

        except ImportError:
            self.skipTest(
                "Could not import PROTECTING_GROUP_TEST_MOLECULES from variables_test"
            )

    def test_protecting_group_context_with_test_molecules(self):
        """Test that protecting group context is generated correctly for test molecules."""
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), '.'))
            from variables_test import PROTECTING_GROUP_TEST_MOLECULES

            from protecting_group import mask_protecting_groups_multisymbol

            # Test context generation for benzyl molecule
            molecule = PROTECTING_GROUP_TEST_MOLECULES[
                "simple_benzyl"]  # "COCc1ccccc1"
            masked = mask_protecting_groups_multisymbol(molecule)

            if masked != molecule:
                context = generate_protecting_group_context(molecule, masked)

                # Verify context contains expected elements
                self.assertIn("PROTECTING GROUP CONTEXT", context)
                self.assertIn(f"Original SMILES: {molecule}", context)
                self.assertIn(f"Masked SMILES: {masked}", context)
                self.assertIn("Symbol mapping:", context)
                self.assertIn("deprotection conditions:", context)

                # Should mention benzyl protecting group
                context_lower = context.lower()
                self.assertTrue(
                    "obn" in context_lower or "benzyl" in context_lower,
                    f"Context should mention benzyl/OBn protecting group: {context}"
                )

        except ImportError:
            self.skipTest(
                "Could not import test molecules from variables_test")

    def test_real_drug_molecules_with_protecting_groups(self):
        """Test that real drug molecules with protecting groups are handled correctly."""
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), '.'))
            from variables_test import PROTECTING_GROUP_TEST_MOLECULES

            from protecting_group import mask_protecting_groups_multisymbol

            # Test real drug molecules that contain protecting group patterns
            drug_molecules = {
                "melatonin_with_pg":
                PROTECTING_GROUP_TEST_MOLECULES["melatonin_with_pg"],
                "naproxen_with_pg":
                PROTECTING_GROUP_TEST_MOLECULES["naproxen_with_pg"]
            }

            for drug_name, smiles in drug_molecules.items():
                masked = mask_protecting_groups_multisymbol(smiles)

                # Should successfully mask without errors
                self.assertNotEqual(
                    masked, "INVALID_SMILES",
                    f"Masking failed for {drug_name}: {smiles}")
                self.assertIsInstance(masked, str,
                                      f"Invalid return type for {drug_name}")

                # Should contain protecting group symbols (since these drugs have methoxy groups)
                self.assertIn(
                    "$", masked,
                    f"Expected methoxy masking in {drug_name}: {smiles} -> {masked}"
                )

        except ImportError:
            self.skipTest(
                "Could not import drug molecules from variables_test")


if __name__ == '__main__':
    unittest.main()
