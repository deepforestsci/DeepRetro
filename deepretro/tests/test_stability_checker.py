"""Tests for deepretro.algorithms.stability_checker."""

import json
import os
import subprocess
import sys
import unittest

import pytest

from deepretro.algorithms.stability_checker import (
    check_molecule_stability,
    is_valid_smiles,
)


class TestIsValidSmiles(unittest.TestCase):
    def test_valid(self):
        self.assertTrue(is_valid_smiles("CCO"))

    def test_invalid(self):
        self.assertFalse(is_valid_smiles("not_a_molecule"))


class TestStabilityChecker(unittest.TestCase):
    def test_invalid_smiles(self):
        res = check_molecule_stability("bad")
        self.assertFalse(res["valid_structure"])
        self.assertEqual(res["stability_score"], 0)

    def test_ethanol_no_issues(self):
        res = check_molecule_stability("CCO")
        self.assertTrue(res["valid_structure"])
        self.assertEqual(res["issues"], [])

    def test_benzene_likely_stable(self):
        res = check_molecule_stability("c1ccccc1")
        self.assertEqual(res["assessment"], "Likely stable")
        self.assertGreaterEqual(res["stability_score"], 80)

    def test_aziridine_flagged(self):
        issues = check_molecule_stability("C1CN1")["issues"]
        self.assertTrue(any("Three-membered heterocycle" in i for i in issues))

    def test_azetidine_flagged(self):
        issues = check_molecule_stability("C1CNC1")["issues"]
        self.assertTrue(any("Four-membered heterocycle" in i for i in issues))

    def test_fused_small_rings(self):
        issues = check_molecule_stability("C1C2CC12")["issues"]
        self.assertTrue(any("Fused" in i for i in issues))

    def test_large_heterocycle(self):
        issues = check_molecule_stability("C1CCNCCCC1")["issues"]
        self.assertTrue(any("Large" in i and "heterocycle" in i for i in issues))

    def test_ethyl_cation_unstable(self):
        res = check_molecule_stability("[CH2+]C")
        self.assertEqual(res["assessment"], "Potentially unstable")
        self.assertTrue(any("carbocation" in i for i in res["issues"]))

    def test_allyl_cation_stabilised(self):
        issues = check_molecule_stability("[CH2+]C=C")["issues"]
        self.assertFalse(any("non-stabilized sp2 carbocation" in i for i in issues))

    def test_cyclobutadiene_anti_aromatic(self):
        issues = check_molecule_stability("C1=CC=C1")["issues"]
        self.assertTrue(any("anti-aromatic" in i for i in issues))


@pytest.mark.parametrize(
    ("smiles", "message"),
    [
        (
            "C1CNCCOCC1",
            "Large (8-membered) heterocycle with N, O (potentially unstable)",
        ),
        (
            "C1CCC2NOC12",
            "5-membered carbon ring fused with 4-membered ring containing N, O "
            "(strained system)",
        ),
    ],
    ids=["large-ring", "fused-cyclopentane"],
)
def test_ring_warnings_list_heteroatoms_alphabetically(
    smiles: str, message: str
) -> None:
    code = (
        "import json, sys;"
        "from deepretro.algorithms.stability_checker import check_molecule_stability;"
        "print(json.dumps(check_molecule_stability(sys.argv[1])['issues']))"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code, smiles],
        # Seed 10 makes a plain {"N", "O"} set iterate as O, N on Python
        # 3.10 to 3.13, whatever order the atoms are added in.
        env=dict(os.environ, PYTHONHASHSEED="10"),
        capture_output=True,
        text=True,
        check=True,
    )
    assert message in json.loads(proc.stdout)


if __name__ == "__main__":
    unittest.main()
