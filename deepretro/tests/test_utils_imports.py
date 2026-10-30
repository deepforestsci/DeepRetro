"""Check lazy utility exports in a fresh interpreter."""

import subprocess
import sys


def test_utils_import_does_not_load_optional_science_stack() -> None:
    """Namespace import stays usable before science dependencies are loaded."""
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import deepretro.utils; "
            "assert 'rdkit' not in sys.modules; assert 'numpy' not in sys.modules",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


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
    import pytest
    import deepretro.utils as utils

    with pytest.raises(AttributeError, match="does_not_exist"):
        getattr(utils, "does_not_exist")
