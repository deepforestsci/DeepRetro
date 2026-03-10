"""Tests for deepretro.algorithms.stability_checker."""

import unittest
from deepretro.algorithms.stability_checker import check_molecule_stability, is_valid_smiles


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


if __name__ == "__main__":
    unittest.main()
