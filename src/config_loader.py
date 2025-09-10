"""
Configuration loader for DeepRetro protecting groups.

This module loads protecting group configurations from external JSON files.
Raises errors if configuration cannot be loaded properly.
"""

import json
from pathlib import Path
from typing import Dict, Tuple, Optional


class ConfigurationError(Exception):
    """Raised when configuration cannot be loaded or is invalid."""
    pass


class ProtectingGroupConfig:
    """Simple configuration manager for protecting groups."""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the config manager.
        
        Parameters
        ----------
        config_path : str, optional
            Path to the configuration file. If None, searches for default locations.
        """
        self.config_path = config_path
        self._pg_map: Optional[Dict[str, Tuple[str, str]]] = None

    @property
    def pg_map(self) -> Dict[str, Tuple[str, str]]:
        """Get the protecting groups mapping."""
        if self._pg_map is None:
            self._load_config()
        if self._pg_map is None:
            raise ConfigurationError(
                "Failed to load protecting groups configuration")
        return self._pg_map

    def _find_config_file(self) -> Optional[str]:
        """Find the configuration file in default locations."""
        if self.config_path and Path(self.config_path).exists():
            return self.config_path

        # Search in default locations
        project_root = Path(__file__).parent.parent
        possible_paths = [
            project_root / "config" / "protecting_groups.json",
            project_root / "protecting_groups.json",
        ]

        for path in possible_paths:
            if path.exists():
                return str(path)

        return None

    def _validate_config_structure(self, config_data: dict) -> None:
        """Validate the structure of the loaded configuration."""
        if not isinstance(config_data, dict):
            raise ConfigurationError("Configuration must be a JSON object")

        if "protecting_groups" not in config_data:
            raise ConfigurationError(
                "Configuration must contain 'protecting_groups' key")

        protecting_groups = config_data["protecting_groups"]
        if not isinstance(protecting_groups, dict):
            raise ConfigurationError("'protecting_groups' must be an object")

        if not protecting_groups:
            raise ConfigurationError("'protecting_groups' cannot be empty")

        for pg_name, pg_data in protecting_groups.items():
            if not isinstance(pg_data, dict):
                raise ConfigurationError(
                    f"Protecting group '{pg_name}' must be an object")

            if "smiles" not in pg_data or "symbol" not in pg_data:
                raise ConfigurationError(
                    f"Protecting group '{pg_name}' must have 'smiles' and 'symbol' keys"
                )

            if not isinstance(pg_data["smiles"], str) or not isinstance(
                    pg_data["symbol"], str):
                raise ConfigurationError(
                    f"'smiles' and 'symbol' for '{pg_name}' must be strings")

            if len(pg_data["symbol"]) != 1:
                raise ConfigurationError(
                    f"Symbol for '{pg_name}' must be a single character, got: '{pg_data['symbol']}'"
                )

            if not pg_data["smiles"].strip():
                raise ConfigurationError(
                    f"SMILES pattern for '{pg_name}' cannot be empty")

            # Validate deprotection field if present (optional)
            if "deprotection" in pg_data:
                if not isinstance(pg_data["deprotection"], str):
                    raise ConfigurationError(
                        f"'deprotection' for '{pg_name}' must be a string")
                if not pg_data["deprotection"].strip():
                    raise ConfigurationError(
                        f"Deprotection condition for '{pg_name}' cannot be empty"
                    )

    def _load_config(self) -> None:
        """Load configuration from file. Raises ConfigurationError if loading fails."""
        config_file = self._find_config_file()

        if not config_file:
            raise ConfigurationError(
                "No configuration file found. Expected 'config/protecting_groups.json' "
                "in project root or 'protecting_groups.json' in current directory."
            )

        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
        except json.JSONDecodeError as e:
            raise ConfigurationError(
                f"Invalid JSON in config file '{config_file}': {e}")
        except IOError as e:
            raise ConfigurationError(
                f"Cannot read config file '{config_file}': {e}")

        # Validate configuration structure
        self._validate_config_structure(config_data)

        # Convert from JSON format to the expected tuple format
        protecting_groups = config_data["protecting_groups"]
        self._pg_map = {
            name: (pg_data["smiles"], pg_data["symbol"])
            for name, pg_data in protecting_groups.items()
        }


# Global instance for easy access
_global_config = ProtectingGroupConfig()


def get_pg_map() -> Dict[str, Tuple[str, str]]:
    """
    Get the protecting groups mapping.
    
    This is the main function that should be used throughout the application
    to access protecting group configurations.
    
    Returns
    -------
    Dict[str, Tuple[str, str]]
        Dictionary mapping protecting group names to (smiles, symbol) tuples.
    """
    return _global_config.pg_map


def reload_config() -> None:
    """Force reload of the configuration."""
    _global_config._pg_map = None


def get_full_pg_config() -> Dict[str, Dict]:
    """
    Get the full protecting groups configuration including all fields.
    
    Returns
    -------
    Dict[str, Dict]
        Dictionary mapping protecting group names to full configuration data.
    """
    config_file = _global_config._find_config_file()
    if not config_file:
        # Return empty dict if no config file
        return {}

    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        return config_data.get("protecting_groups", {})
    except (json.JSONDecodeError, IOError):
        return {}


def generate_protecting_group_context(molecule: str,
                                      masked_smiles: str) -> str:
    """
    Generate dynamic protecting group context based on current configuration.
    
    Parameters
    ----------
    molecule : str
        Original SMILES string
    masked_smiles : str  
        SMILES string with protecting groups masked
    
    Returns
    -------
    str
        Formatted context string for LLM prompts
    """
    # Get the full config data (not just the pg_map tuples)
    full_config = get_full_pg_config()

    # Build symbol mapping and deprotection conditions from config
    symbol_mapping = {}  # symbol -> display_name
    deprotection_conditions = {}  # (display_name, symbol) -> condition

    for pg_name, pg_data in full_config.items():
        # Clean up protecting group name for display
        display_name = pg_name.replace('_', ' ').replace('OEt CCO',
                                                         'OEt').replace(
                                                             'OEt COC', 'OEt')
        if display_name.endswith(' CCO') or display_name.endswith(' COC'):
            display_name = display_name.rsplit(' ', 1)[0]

        symbol = pg_data.get('symbol', '?')

        # Store unique symbol mappings
        symbol_mapping[symbol] = display_name

        # Add deprotection conditions from config (or use fallback)
        condition_key = (display_name, symbol)
        if condition_key not in deprotection_conditions:
            deprotection = pg_data.get('deprotection',
                                       'Standard deprotection conditions')
            deprotection_conditions[
                condition_key] = f"   - {display_name} ({symbol}): {deprotection}"

    # Convert to formatted strings
    symbol_mapping_lines = [
        f"- '{symbol}' represents {display_name} groups"
        for symbol, display_name in sorted(symbol_mapping.items())
    ]
    symbol_mapping_str = '\n'.join(symbol_mapping_lines)
    deprotection_info = '\n'.join(sorted(deprotection_conditions.values()))

    return f"""
PROTECTING GROUP CONTEXT:
The molecule contains protecting groups that have been masked:
Original SMILES: {molecule}
Masked SMILES: {masked_smiles}

Symbol mapping:
{symbol_mapping_str}

IMPORTANT: In your retrosynthetic analysis:
1. Return ACTUAL SMILES strings (not the masked symbols) in the 'data' field
2. Include deprotection steps in your explanations when appropriate
3. Consider the protecting groups as synthetic handles that may need removal
4. Suggest specific deprotection conditions:
{deprotection_info}
5. Consider protecting group compatibility and orthogonal deprotection strategies"""
