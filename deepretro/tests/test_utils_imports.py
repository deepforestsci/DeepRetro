"""Check lazy utility exports in a fresh interpreter."""

import subprocess
import sys

import pytest


def test_utils_import_does_not_load_rdkit_or_numpy() -> None:
    """Importing the namespace does not import RDKit or NumPy."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import deepretro.utils; "
            "assert 'rdkit' not in sys.modules; assert 'numpy' not in sys.modules",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_exports_resolve_to_public_implementations() -> None:
    """Each declared export remains available through the utility namespace."""
    import deepretro.utils as utils
    from deepretro.utils import domain_features, metrics

    assert utils.NUM_DOMAIN_FEATURES == domain_features.NUM_DOMAIN_FEATURES
    assert (
        utils.extract_domain_features_single
        is domain_features.extract_domain_features_single
    )
    assert utils.find_optimal_threshold is metrics.find_optimal_threshold


def test_unknown_export_raises_attribute_error() -> None:
    """Unknown names retain normal module attribute behavior."""
    import deepretro.utils as utils

    with pytest.raises(AttributeError, match="does_not_exist"):
        getattr(utils, "does_not_exist")
