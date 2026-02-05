import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add the src directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from protecting_group import (
    FUNCTIONAL_GROUPS,
    PROTECTING_GROUPS,
    ProtectionSite,
    PGSuggestion,
    ProtectionRecommendation,
    ProtectionSelection,
    ProtectionPlan,
    ActiveProtection,
    ProtectionState,
    canonicalize_and_map_atoms,
    _generate_site_id,
    identify_protection_sites,
    suggest_protecting_groups,
    get_protection_recommendations,
    format_recommendations_for_prompt,
    # HITL API functions
    get_sites_for_selection,
    get_pg_options_for_site,
    create_protection_selection,
    create_protection_plan,
    hitl_workflow_example,
    format_sites_for_display,
    format_pg_options_for_display,
)
from rdkit import Chem


class TestAtomMapping(unittest.TestCase):
    """Tests for the canonicalize_and_map_atoms function."""

    def test_basic_mapping(self):
        """Test that atom mapping returns valid data."""
        mapped_smiles, atom_map, mol = canonicalize_and_map_atoms("CCO")
        self.assertIsInstance(mapped_smiles, str)
        self.assertGreater(len(mapped_smiles), 0)
        self.assertIsInstance(atom_map, dict)
        self.assertGreater(len(atom_map), 0)
        self.assertIsNotNone(mol)

    def test_invalid_smiles(self):
        """Test that invalid SMILES returns empty results."""
        mapped_smiles, atom_map, mol = canonicalize_and_map_atoms(
            "invalid_smiles")
        self.assertEqual(mapped_smiles, "")
        self.assertEqual(atom_map, {})
        self.assertIsNone(mol)

    def test_empty_smiles(self):
        """Test that empty SMILES is handled gracefully."""
        # RDKit may parse empty string as a valid (empty) molecule,
        # so we just verify no crash and atom_map is empty or minimal
        mapped_smiles, atom_map, mol = canonicalize_and_map_atoms("")
        # Either it fails gracefully or returns an empty/trivial result
        self.assertIsInstance(mapped_smiles, str)
        self.assertIsInstance(atom_map, dict)

    def test_atom_map_numbers_are_positive(self):
        """Test that all atom map numbers are positive integers."""
        _, atom_map, _ = canonicalize_and_map_atoms("CCO")
        for map_num in atom_map:
            self.assertIsInstance(map_num, int)
            self.assertGreater(map_num, 0)

    def test_atom_map_has_symbol_info(self):
        """Test that atom map contains symbol information."""
        _, atom_map, _ = canonicalize_and_map_atoms("CCO")
        for map_num, info in atom_map.items():
            self.assertIn("symbol", info)
            self.assertIn("idx", info)

    def test_stability_across_representations(self):
        """Test that different SMILES of the same molecule produce same atom count."""
        # These are all ethanol
        smiles_variants = ["CCO", "OCC", "C(O)C"]
        atom_counts = set()
        for smiles in smiles_variants:
            _, atom_map, _ = canonicalize_and_map_atoms(smiles)
            atom_counts.add(len(atom_map))
        # All variants should have the same number of heavy atoms
        self.assertEqual(len(atom_counts), 1)

    def test_canonical_determinism(self):
        """Test that the same SMILES always produces the same mapped SMILES."""
        s1, _, _ = canonicalize_and_map_atoms("CCO")
        s2, _, _ = canonicalize_and_map_atoms("CCO")
        self.assertEqual(s1, s2)

    def test_mapped_smiles_contains_map_numbers(self):
        """Test that mapped SMILES contains atom map notation."""
        mapped_smiles, _, _ = canonicalize_and_map_atoms("CCO")
        # Should contain at least one :N] pattern
        self.assertRegex(mapped_smiles, r":\d+\]")

    def test_stereochemistry_molecule(self):
        """Test atom mapping with stereochemistry."""
        smiles = "C[C@H](O)C(=O)OC"
        mapped_smiles, atom_map, mol = canonicalize_and_map_atoms(smiles)
        self.assertGreater(len(mapped_smiles), 0)
        self.assertIsNotNone(mol)
        self.assertGreater(len(atom_map), 0)


class TestSiteIdGeneration(unittest.TestCase):
    """Tests for the _generate_site_id function."""

    def test_deterministic(self):
        """Test that same inputs produce same site_id."""
        id1 = _generate_site_id("primary_alcohol", (1, 2))
        id2 = _generate_site_id("primary_alcohol", (1, 2))
        self.assertEqual(id1, id2)

    def test_order_independent(self):
        """Test that atom map order doesn't matter (sorted internally)."""
        id1 = _generate_site_id("primary_alcohol", (2, 1))
        id2 = _generate_site_id("primary_alcohol", (1, 2))
        self.assertEqual(id1, id2)

    def test_different_fg_different_id(self):
        """Test that different functional groups produce different IDs."""
        id1 = _generate_site_id("primary_alcohol", (1, 2))
        id2 = _generate_site_id("secondary_alcohol", (1, 2))
        self.assertNotEqual(id1, id2)

    def test_different_atoms_different_id(self):
        """Test that different atoms produce different IDs."""
        id1 = _generate_site_id("primary_alcohol", (1, 2))
        id2 = _generate_site_id("primary_alcohol", (3, 4))
        self.assertNotEqual(id1, id2)

    def test_returns_string(self):
        """Test that site_id is a string."""
        site_id = _generate_site_id("primary_alcohol", (1, 2))
        self.assertIsInstance(site_id, str)
        self.assertGreater(len(site_id), 0)


class TestFunctionalGroupsDatabase(unittest.TestCase):
    """Tests for the FUNCTIONAL_GROUPS database."""

    def test_database_not_empty(self):
        """Test that the database has entries."""
        self.assertGreater(len(FUNCTIONAL_GROUPS), 0)

    def test_all_entries_have_required_fields(self):
        """Test that all entries have the required fields."""
        required_fields = [
            "smarts", "name", "category", "reactivity", "compatible_pgs"
        ]
        for fg_id, fg_data in FUNCTIONAL_GROUPS.items():
            for field in required_fields:
                self.assertIn(
                    field, fg_data,
                    f"Missing field '{field}' in functional group '{fg_id}'")

    def test_all_smarts_are_valid(self):
        """Test that all SMARTS patterns are valid RDKit patterns."""
        for fg_id, fg_data in FUNCTIONAL_GROUPS.items():
            pattern = Chem.MolFromSmarts(fg_data["smarts"])
            self.assertIsNotNone(
                pattern,
                f"Invalid SMARTS pattern for '{fg_id}': {fg_data['smarts']}")

    def test_reactivity_values_are_valid(self):
        """Test that reactivity values are from expected set."""
        valid_reactivity = {"low", "medium", "high", "very_high"}
        for fg_id, fg_data in FUNCTIONAL_GROUPS.items():
            self.assertIn(
                fg_data["reactivity"], valid_reactivity,
                f"Invalid reactivity '{fg_data['reactivity']}' for '{fg_id}'")

    def test_categories_exist(self):
        """Test that expected categories are present."""
        categories = {fg["category"] for fg in FUNCTIONAL_GROUPS.values()}
        expected = {
            "alcohol", "phenol", "amine", "carbonyl", "carboxylic_acid",
            "thiol", "alkyne", "diol"
        }
        for cat in expected:
            self.assertIn(cat, categories,
                          f"Missing category: {cat}. Found: {categories}")

    def test_compatible_pgs_reference_valid_pgs(self):
        """Test that compatible_pgs reference valid protecting groups."""
        for fg_id, fg_data in FUNCTIONAL_GROUPS.items():
            for pg_id in fg_data["compatible_pgs"]:
                self.assertIn(
                    pg_id, PROTECTING_GROUPS,
                    f"Unknown PG '{pg_id}' in compatible_pgs of '{fg_id}'")

    def test_new_patterns_exist(self):
        """Test that the newly added patterns are present."""
        new_patterns = [
            "ester", "lactone", "tertiary_amine", "secondary_amide",
            "tertiary_amide", "epoxide"
        ]
        for pattern_id in new_patterns:
            self.assertIn(pattern_id, FUNCTIONAL_GROUPS,
                          f"Missing new pattern: {pattern_id}")


class TestProtectingGroupsDatabase(unittest.TestCase):
    """Tests for the PROTECTING_GROUPS database."""

    def test_database_not_empty(self):
        """Test that the database has entries."""
        self.assertGreater(len(PROTECTING_GROUPS), 0)

    def test_all_entries_have_required_fields(self):
        """Test that all entries have the required fields."""
        required_fields = [
            "name", "abbreviation", "protection_reagent", "orthogonality_score",
            "cost_score"
        ]
        for pg_id, pg_data in PROTECTING_GROUPS.items():
            for field in required_fields:
                self.assertIn(
                    field, pg_data,
                    f"Missing field '{field}' in protecting group '{pg_id}'")

    def test_stability_has_expected_keys(self):
        """Test that stability dict has expected keys."""
        expected_keys = {"acid", "base", "hydrogenation", "oxidation"}
        for pg_id, pg_data in PROTECTING_GROUPS.items():
            stability = pg_data.get("stability", {})
            if stability:
                for key in expected_keys:
                    self.assertIn(
                        key, stability,
                        f"Missing stability key '{key}' in PG '{pg_id}'")

    def test_orthogonality_score_in_valid_range(self):
        """Test that orthogonality scores are between 0 and 1."""
        for pg_id, pg_data in PROTECTING_GROUPS.items():
            score = pg_data.get("orthogonality_score", 0.5)
            self.assertGreaterEqual(score, 0.0,
                                    f"Score below 0 for '{pg_id}'")
            self.assertLessEqual(score, 1.0, f"Score above 1 for '{pg_id}'")

    def test_cost_score_in_valid_range(self):
        """Test that cost scores are between 0 and 1."""
        for pg_id, pg_data in PROTECTING_GROUPS.items():
            score = pg_data.get("cost_score", 0.5)
            self.assertGreaterEqual(score, 0.0,
                                    f"Cost score below 0 for '{pg_id}'")
            self.assertLessEqual(score, 1.0,
                                 f"Cost score above 1 for '{pg_id}'")


class TestIdentifyProtectionSites(unittest.TestCase):
    """Tests for the identify_protection_sites function."""

    def test_invalid_smiles_returns_empty_list(self):
        """Test that invalid SMILES returns empty list."""
        result = identify_protection_sites("invalid_smiles")
        self.assertEqual(result, [])

    def test_empty_smiles_returns_empty_list(self):
        """Test that empty SMILES returns empty list."""
        result = identify_protection_sites("")
        self.assertEqual(result, [])

    def test_identifies_primary_alcohol(self):
        """Test identification of primary alcohol."""
        sites = identify_protection_sites("CCO")
        self.assertGreater(len(sites), 0)
        alcohol_sites = [s for s in sites if s.category == "alcohol"]
        self.assertGreater(len(alcohol_sites), 0)

    def test_identifies_secondary_alcohol(self):
        """Test identification of secondary alcohol."""
        sites = identify_protection_sites("CC(O)C")
        self.assertGreater(len(sites), 0)
        alcohol_sites = [s for s in sites if s.category == "alcohol"]
        self.assertGreater(len(alcohol_sites), 0)

    def test_identifies_phenol(self):
        """Test identification of phenol."""
        sites = identify_protection_sites("c1ccc(O)cc1")
        self.assertGreater(len(sites), 0)
        phenol_sites = [s for s in sites if s.category == "phenol"]
        self.assertGreater(len(phenol_sites), 0)

    def test_identifies_primary_amine(self):
        """Test identification of primary amine."""
        sites = identify_protection_sites("CCN")
        self.assertGreater(len(sites), 0)
        amine_sites = [s for s in sites if s.category == "amine"]
        self.assertGreater(len(amine_sites), 0)

    def test_identifies_aldehyde(self):
        """Test identification of aldehyde."""
        sites = identify_protection_sites("CC=O")
        self.assertGreater(len(sites), 0)
        carbonyl_sites = [s for s in sites if s.category == "carbonyl"]
        self.assertGreater(len(carbonyl_sites), 0)

    def test_identifies_carboxylic_acid(self):
        """Test identification of carboxylic acid."""
        sites = identify_protection_sites("CC(=O)O")
        self.assertGreater(len(sites), 0)
        acid_sites = [s for s in sites if s.category == "carboxylic_acid"]
        self.assertGreater(len(acid_sites), 0)

    def test_identifies_thiol(self):
        """Test identification of thiol."""
        sites = identify_protection_sites("CCS")
        self.assertGreater(len(sites), 0)
        thiol_sites = [s for s in sites if s.category == "thiol"]
        self.assertGreater(len(thiol_sites), 0)

    def test_identifies_terminal_alkyne(self):
        """Test identification of terminal alkyne."""
        sites = identify_protection_sites("CC#C")
        self.assertGreater(len(sites), 0)
        alkyne_sites = [s for s in sites if s.category == "alkyne"]
        self.assertGreater(len(alkyne_sites), 0)

    def test_identifies_vicinal_diol(self):
        """Test identification of 1,2-diol."""
        sites = identify_protection_sites("CC(O)CO")
        self.assertGreater(len(sites), 0)
        alcohol_sites = [s for s in sites if s.category in ["alcohol", "diol"]]
        self.assertGreater(len(alcohol_sites), 0)

    def test_returns_protection_site_objects(self):
        """Test that function returns ProtectionSite objects with new fields."""
        sites = identify_protection_sites("CCO")

        for site in sites:
            self.assertIsInstance(site, ProtectionSite)
            self.assertIsInstance(site.site_id, str)
            self.assertGreater(len(site.site_id), 0)
            self.assertIsInstance(site.atom_map_numbers, tuple)
            self.assertIsInstance(site.functional_group_id, str)
            self.assertIsInstance(site.functional_group_name, str)
            self.assertIsInstance(site.category, str)
            self.assertIsInstance(site.reactivity, str)
            self.assertIsInstance(site.compatible_pgs, list)

    def test_multiple_functional_groups(self):
        """Test identification of multiple functional groups in same molecule."""
        sites = identify_protection_sites("OCC=O")  # glycolaldehyde
        self.assertGreater(len(sites), 0)

        sites2 = identify_protection_sites("NCCO")  # ethanolamine
        categories2 = {s.category for s in sites2}
        self.assertGreater(len(sites2), 0)
        self.assertTrue(
            "amine" in categories2 or "alcohol" in categories2,
            f"Expected amine or alcohol category, got: {categories2}")

    def test_accepts_premapped_molecule(self):
        """Test that function works with pre-mapped molecule."""
        _, atom_map, mapped_mol = canonicalize_and_map_atoms("CCO")
        sites = identify_protection_sites("CCO",
                                          mapped_mol=mapped_mol,
                                          atom_map=atom_map)
        self.assertGreater(len(sites), 0)


class TestCategoryDedup(unittest.TestCase):
    """Tests for category-aware deduplication in identify_protection_sites."""

    def test_cross_category_sites_both_found(self):
        """Test that sites from different categories on the same molecule are all found."""
        # Molecule with alcohol + carbonyl (ketone): 4-hydroxy-2-butanone
        sites = identify_protection_sites("CC(=O)CCO")
        categories = {s.category for s in sites}

        # Should find both alcohol and carbonyl categories
        self.assertIn("alcohol", categories,
                      f"Missing alcohol. Found categories: {categories}")
        self.assertIn("carbonyl", categories,
                      f"Missing carbonyl. Found categories: {categories}")

    def test_alcohol_and_amine_both_found(self):
        """Test that alcohol and amine on the same molecule are both found."""
        sites = identify_protection_sites("NCCO")  # ethanolamine
        categories = {s.category for s in sites}
        self.assertIn("amine", categories)
        self.assertIn("alcohol", categories)

    def test_same_category_no_duplicates(self):
        """Test that within the same category, overlapping patterns don't duplicate."""
        # Simple primary alcohol - should get exactly one alcohol site
        sites = identify_protection_sites("CCO")
        alcohol_sites = [s for s in sites if s.category == "alcohol"]
        # Should have at most a few matches, not an explosion
        self.assertLessEqual(len(alcohol_sites), 3)

    def test_identifies_ester(self):
        """Test identification of ester group."""
        sites = identify_protection_sites("CC(=O)OCC")  # ethyl acetate
        categories = {s.category for s in sites}
        self.assertIn("ester", categories,
                      f"Missing ester. Found categories: {categories}")

    def test_identifies_tertiary_amine(self):
        """Test identification of tertiary amine."""
        sites = identify_protection_sites("CN(C)C")  # trimethylamine
        amine_sites = [s for s in sites if s.category == "amine"]
        self.assertGreater(len(amine_sites), 0)
        # Should find tertiary amine
        fg_ids = {s.functional_group_id for s in amine_sites}
        self.assertTrue(
            "tertiary_amine" in fg_ids
            or "tertiary_amine_aromatic" in fg_ids,
            f"Expected tertiary_amine, got: {fg_ids}")


class TestSuggestProtectingGroups(unittest.TestCase):
    """Tests for the suggest_protecting_groups function."""

    def test_returns_suggestions_for_valid_site(self):
        """Test that function returns suggestions for valid site."""
        sites = identify_protection_sites("CCO")
        self.assertGreater(len(sites), 0)
        suggestions = suggest_protecting_groups(sites[0])
        self.assertGreater(len(suggestions), 0)

    def test_returns_pg_suggestion_objects(self):
        """Test that function returns PGSuggestion objects."""
        sites = identify_protection_sites("CCO")
        suggestions = suggest_protecting_groups(sites[0])
        for sug in suggestions:
            self.assertIsInstance(sug, PGSuggestion)
            self.assertIsInstance(sug.protecting_group_id, str)
            self.assertIsInstance(sug.name, str)
            self.assertIsInstance(sug.abbreviation, str)
            self.assertIsInstance(sug.score, float)
            self.assertIsInstance(sug.protection_reagent, str)
            self.assertIsInstance(sug.deprotection_conditions, list)

    def test_suggestions_are_sorted_by_score(self):
        """Test that suggestions are sorted by score (highest first)."""
        sites = identify_protection_sites("CCO")
        suggestions = suggest_protecting_groups(sites[0])
        if len(suggestions) > 1:
            scores = [s.score for s in suggestions]
            self.assertEqual(scores, sorted(scores, reverse=True))

    def test_max_suggestions_parameter(self):
        """Test that max_suggestions parameter limits results."""
        sites = identify_protection_sites("CCO")
        suggestions = suggest_protecting_groups(sites[0], max_suggestions=2)
        self.assertLessEqual(len(suggestions), 2)

    def test_scores_are_in_valid_range(self):
        """Test that all scores are between 0 and 1."""
        sites = identify_protection_sites("CCO")
        suggestions = suggest_protecting_groups(sites[0])
        for sug in suggestions:
            self.assertGreaterEqual(sug.score, 0.0)
            self.assertLessEqual(sug.score, 1.0)


class TestGetProtectionRecommendations(unittest.TestCase):
    """Tests for the get_protection_recommendations function."""

    def test_auto_mode_returns_single_suggestion_per_site(self):
        """Test that auto mode returns single suggestion per site."""
        recs = get_protection_recommendations("CCO", mode="auto")
        for rec in recs:
            self.assertLessEqual(len(rec.suggestions), 1)

    def test_hitl_mode_returns_multiple_suggestions(self):
        """Test that HITL mode returns multiple suggestions."""
        recs = get_protection_recommendations("CCO", mode="hitl")
        self.assertGreater(len(recs), 0)

    def test_returns_protection_recommendation_objects(self):
        """Test that function returns ProtectionRecommendation objects."""
        recs = get_protection_recommendations("CCO", mode="auto")
        for rec in recs:
            self.assertIsInstance(rec, ProtectionRecommendation)
            self.assertIsInstance(rec.site, ProtectionSite)
            self.assertIsInstance(rec.suggestions, list)

    def test_invalid_smiles_returns_empty_list(self):
        """Test that invalid SMILES returns empty list."""
        recs = get_protection_recommendations("invalid", mode="auto")
        self.assertEqual(recs, [])

    def test_molecule_with_no_protectable_sites(self):
        """Test molecule with no protectable functional groups."""
        recs = get_protection_recommendations("C", mode="auto")
        self.assertEqual(len(recs), 0)

    def test_accepts_premapped_molecule(self):
        """Test that function works with pre-mapped molecule."""
        _, atom_map, mapped_mol = canonicalize_and_map_atoms("CCO")
        recs = get_protection_recommendations("CCO",
                                              mode="hitl",
                                              mapped_mol=mapped_mol,
                                              atom_map=atom_map)
        self.assertGreater(len(recs), 0)

    def test_existing_pgs_affects_scoring(self):
        """Test that passing existing_pgs changes scores (orthogonality)."""
        recs_without = get_protection_recommendations("CCO", mode="hitl")
        recs_with = get_protection_recommendations("CCO",
                                                   mode="hitl",
                                                   existing_pgs=["TBS"])
        # The scores should differ when existing PGs are present
        # (orthogonality penalty applied)
        if recs_without and recs_with:
            scores_without = [
                s.score for rec in recs_without for s in rec.suggestions
            ]
            scores_with = [
                s.score for rec in recs_with for s in rec.suggestions
            ]
            # At least some scores should differ
            self.assertTrue(
                any(a != b
                    for a, b in zip(scores_without, scores_with))
                or len(scores_without) != len(scores_with),
                "Expected different scores with existing PGs")


class TestFormatRecommendationsForPrompt(unittest.TestCase):
    """Tests for the format_recommendations_for_prompt function."""

    def test_empty_recommendations(self):
        """Test formatting with no recommendations."""
        result = format_recommendations_for_prompt([])
        self.assertIn("No protection sites", result)

    def test_auto_mode_formatting(self):
        """Test formatting in auto mode."""
        recs = get_protection_recommendations("CCO", mode="auto")
        result = format_recommendations_for_prompt(recs, mode="auto")
        self.assertIsInstance(result, str)
        self.assertIn("IDENTIFIED PROTECTION SITES", result)

    def test_hitl_mode_formatting(self):
        """Test formatting in HITL mode."""
        recs = get_protection_recommendations("CCO", mode="hitl")
        result = format_recommendations_for_prompt(recs, mode="hitl")
        self.assertIn("Available protecting groups", result)

    def test_includes_atom_map_numbers(self):
        """Test that formatted output includes atom map numbers."""
        recs = get_protection_recommendations("CCO", mode="auto")
        result = format_recommendations_for_prompt(recs, mode="auto")
        self.assertIn("Atom map numbers", result)

    def test_includes_site_id(self):
        """Test that formatted output includes site ID."""
        recs = get_protection_recommendations("CCO", mode="auto")
        result = format_recommendations_for_prompt(recs, mode="auto")
        self.assertIn("id:", result)

    def test_returns_string(self):
        """Test that function returns a string."""
        recs = get_protection_recommendations("CCO", mode="hitl")
        result = format_recommendations_for_prompt(recs, mode="hitl")
        self.assertIsInstance(result, str)


class TestDataClasses(unittest.TestCase):
    """Tests for the data classes."""

    def test_protection_site_creation(self):
        """Test ProtectionSite dataclass creation."""
        site = ProtectionSite(site_id="abc123",
                              atom_map_numbers=(1, 2),
                              functional_group_id="primary_alcohol",
                              functional_group_name="Primary Alcohol",
                              category="alcohol",
                              reactivity="high",
                              compatible_pgs=["TBS", "Bn"])

        self.assertEqual(site.site_id, "abc123")
        self.assertEqual(site.atom_map_numbers, (1, 2))
        self.assertEqual(site.functional_group_id, "primary_alcohol")
        self.assertEqual(site.category, "alcohol")
        self.assertEqual(site.reactivity, "high")
        self.assertEqual(len(site.compatible_pgs), 2)

    def test_pg_suggestion_creation(self):
        """Test PGSuggestion dataclass creation."""
        sug = PGSuggestion(protecting_group_id="TBS",
                           name="tert-Butyldimethylsilyl",
                           abbreviation="TBS",
                           score=0.85,
                           protection_reagent="TBSCl, imidazole",
                           deprotection_conditions=["TBAF"],
                           stability={
                               "acid": "low",
                               "base": "high"
                           },
                           orthogonality_score=0.9)

        self.assertEqual(sug.protecting_group_id, "TBS")
        self.assertEqual(sug.score, 0.85)
        self.assertEqual(len(sug.deprotection_conditions), 1)

    def test_protection_recommendation_creation(self):
        """Test ProtectionRecommendation dataclass creation."""
        site = ProtectionSite(site_id="abc123",
                              atom_map_numbers=(1, ),
                              functional_group_id="primary_alcohol",
                              functional_group_name="Primary Alcohol",
                              category="alcohol",
                              reactivity="high")

        sug = PGSuggestion(protecting_group_id="TBS",
                           name="tert-Butyldimethylsilyl",
                           abbreviation="TBS",
                           score=0.85,
                           protection_reagent="TBSCl, imidazole")

        rec = ProtectionRecommendation(site=site, suggestions=[sug])

        self.assertIsInstance(rec.site, ProtectionSite)
        self.assertEqual(len(rec.suggestions), 1)


class TestProtectionState(unittest.TestCase):
    """Tests for the ProtectionState dataclass."""

    def test_empty_state_creation(self):
        """Test creating an empty protection state."""
        state = ProtectionState()
        self.assertEqual(len(state.active_protections), 0)
        self.assertEqual(state.get_active_pg_ids(), [])
        self.assertEqual(state.get_currently_active(), [])

    def test_add_protection(self):
        """Test adding a protection."""
        state = ProtectionState()
        state.add_protection(site_id="abc123",
                             pg_id="TBS",
                             atom_map_numbers=(1, 2),
                             step="1")
        self.assertEqual(len(state.active_protections), 1)
        self.assertEqual(state.get_active_pg_ids(), ["TBS"])

    def test_remove_protection(self):
        """Test removing a protection."""
        state = ProtectionState()
        state.add_protection("abc123", "TBS", (1, 2), "1")
        result = state.remove_protection("abc123", "3")
        self.assertTrue(result)
        self.assertEqual(state.get_active_pg_ids(), [])
        self.assertEqual(state.active_protections[0].step_removed, "3")

    def test_remove_nonexistent_protection(self):
        """Test removing a non-existent protection returns False."""
        state = ProtectionState()
        result = state.remove_protection("nonexistent", "1")
        self.assertFalse(result)

    def test_get_currently_active(self):
        """Test getting currently active protections."""
        state = ProtectionState()
        state.add_protection("site1", "TBS", (1, 2), "1")
        state.add_protection("site2", "Bn", (3, 4), "1")
        state.remove_protection("site1", "3")

        active = state.get_currently_active()
        self.assertEqual(len(active), 1)
        self.assertEqual(active[0].protecting_group_id, "Bn")

    def test_get_active_pg_ids(self):
        """Test getting active PG IDs."""
        state = ProtectionState()
        state.add_protection("site1", "TBS", (1, 2), "1")
        state.add_protection("site2", "Bn", (3, 4), "2")

        ids = state.get_active_pg_ids()
        self.assertEqual(set(ids), {"TBS", "Bn"})

    def test_serialization_round_trip(self):
        """Test to_dict/from_dict round trip."""
        state = ProtectionState()
        state.add_protection("abc123", "TBS", (1, 2), "1")
        state.add_protection("def456", "Bn", (3, 4), "2")
        state.remove_protection("abc123", "3")

        data = state.to_dict()
        restored = ProtectionState.from_dict(data)

        self.assertEqual(len(restored.active_protections), 2)
        self.assertEqual(restored.active_protections[0].site_id, "abc123")
        self.assertEqual(restored.active_protections[0].step_removed, "3")
        self.assertEqual(restored.active_protections[1].protecting_group_id,
                         "Bn")
        self.assertIsNone(restored.active_protections[1].step_removed)

    def test_from_dict_empty(self):
        """Test from_dict with empty/None input."""
        state = ProtectionState.from_dict(None)
        self.assertEqual(len(state.active_protections), 0)

        state2 = ProtectionState.from_dict({})
        self.assertEqual(len(state2.active_protections), 0)

    def test_serialization_preserves_atom_map_numbers(self):
        """Test that serialization preserves atom_map_numbers as tuples."""
        state = ProtectionState()
        state.add_protection("abc", "TBS", (5, 7), "1")

        data = state.to_dict()
        self.assertEqual(data["active_protections"][0]["atom_map_numbers"],
                         [5, 7])

        restored = ProtectionState.from_dict(data)
        self.assertEqual(restored.active_protections[0].atom_map_numbers,
                         (5, 7))


class TestSMARTSPatternMatching(unittest.TestCase):
    """Tests for specific SMARTS pattern matching."""

    def test_primary_alcohol_pattern(self):
        """Test primary alcohol SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["primary_alcohol"]["smarts"])
        mol1 = Chem.MolFromSmiles("CCO")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_phenol_pattern(self):
        """Test phenol SMARTS pattern."""
        pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUPS["phenol"]["smarts"])
        mol1 = Chem.MolFromSmiles("c1ccc(O)cc1")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)
        mol2 = Chem.MolFromSmiles("CCO")
        self.assertEqual(len(mol2.GetSubstructMatches(pattern)), 0)

    def test_aldehyde_pattern(self):
        """Test aldehyde SMARTS pattern."""
        pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUPS["aldehyde"]["smarts"])
        mol1 = Chem.MolFromSmiles("CC=O")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_carboxylic_acid_pattern(self):
        """Test carboxylic acid SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["carboxylic_acid"]["smarts"])
        mol1 = Chem.MolFromSmiles("CC(=O)O")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_primary_amine_pattern(self):
        """Test primary amine SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["primary_amine"]["smarts"])
        mol1 = Chem.MolFromSmiles("CCN")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_thiol_pattern(self):
        """Test thiol SMARTS pattern."""
        pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUPS["thiol"]["smarts"])
        mol1 = Chem.MolFromSmiles("CCS")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_terminal_alkyne_pattern(self):
        """Test terminal alkyne SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["terminal_alkyne"]["smarts"])
        mol1 = Chem.MolFromSmiles("CC#C")
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_ester_pattern(self):
        """Test ester SMARTS pattern."""
        pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUPS["ester"]["smarts"])
        mol1 = Chem.MolFromSmiles("CC(=O)OCC")  # ethyl acetate
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_tertiary_amine_pattern(self):
        """Test tertiary amine SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["tertiary_amine"]["smarts"])
        mol1 = Chem.MolFromSmiles("CN(C)C")  # trimethylamine
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)

    def test_secondary_amide_pattern(self):
        """Test secondary amide SMARTS pattern."""
        pattern = Chem.MolFromSmarts(
            FUNCTIONAL_GROUPS["secondary_amide"]["smarts"])
        mol1 = Chem.MolFromSmiles("CC(=O)NC")  # N-methylacetamide
        self.assertGreater(len(mol1.GetSubstructMatches(pattern)), 0)


class TestHITLAPI(unittest.TestCase):
    """Tests for the Human-in-the-Loop API functions."""

    def test_get_sites_for_selection_returns_list_of_dicts(self):
        """Test that get_sites_for_selection returns list of dicts."""
        sites = get_sites_for_selection("CCO")

        self.assertIsInstance(sites, list)
        self.assertGreater(len(sites), 0)

        for site in sites:
            self.assertIsInstance(site, dict)
            self.assertIn('index', site)
            self.assertIn('site', site)
            self.assertIn('site_id', site)
            self.assertIn('functional_group', site)
            self.assertIn('category', site)
            self.assertIn('reactivity', site)
            self.assertIn('atom_map_numbers', site)

    def test_get_sites_for_selection_indices_are_sequential(self):
        """Test that site indices are sequential starting from 0."""
        sites = get_sites_for_selection("CC(O)C=O")
        indices = [s['index'] for s in sites]
        expected = list(range(len(sites)))
        self.assertEqual(indices, expected)

    def test_get_sites_for_selection_empty_for_invalid_smiles(self):
        """Test that invalid SMILES returns empty list."""
        sites = get_sites_for_selection("invalid")
        self.assertEqual(sites, [])

    def test_site_ids_are_unique(self):
        """Test that site IDs are unique across sites."""
        sites = get_sites_for_selection("NCCO")
        site_ids = [s['site_id'] for s in sites]
        self.assertEqual(len(site_ids), len(set(site_ids)))

    def test_get_pg_options_for_site_returns_list_of_dicts(self):
        """Test that get_pg_options_for_site returns list of dicts."""
        sites = get_sites_for_selection("CCO")
        self.assertGreater(len(sites), 0)

        options = get_pg_options_for_site(sites[0]['site'])

        self.assertIsInstance(options, list)
        self.assertGreater(len(options), 0)

        for opt in options:
            self.assertIsInstance(opt, dict)
            self.assertIn('index', opt)
            self.assertIn('suggestion', opt)
            self.assertIn('name', opt)
            self.assertIn('abbreviation', opt)
            self.assertIn('score', opt)
            self.assertIn('protection_reagent', opt)
            self.assertIn('deprotection', opt)

    def test_get_pg_options_respects_max_options(self):
        """Test that max_options parameter limits results."""
        sites = get_sites_for_selection("CCO")
        options = get_pg_options_for_site(sites[0]['site'], max_options=3)
        self.assertLessEqual(len(options), 3)

    def test_create_protection_selection(self):
        """Test creating a ProtectionSelection."""
        sites = get_sites_for_selection("CCO")
        options = get_pg_options_for_site(sites[0]['site'])

        selection = create_protection_selection(sites[0], options[0])

        self.assertIsInstance(selection, ProtectionSelection)
        self.assertEqual(selection.site_index, 0)
        self.assertIsInstance(selection.site, ProtectionSite)
        self.assertIsInstance(selection.selected_pg, PGSuggestion)

    def test_create_protection_plan(self):
        """Test creating a ProtectionPlan."""
        smiles = "CCO"
        sites = get_sites_for_selection(smiles)
        options = get_pg_options_for_site(sites[0]['site'])
        selection = create_protection_selection(sites[0], options[0])

        plan = create_protection_plan(smiles, [selection])

        self.assertIsInstance(plan, ProtectionPlan)
        self.assertEqual(plan.smiles, smiles)
        self.assertEqual(len(plan.selections), 1)

    def test_protection_plan_get_all_pgs(self):
        """Test ProtectionPlan.get_all_pgs method."""
        smiles = "CCO"
        sites = get_sites_for_selection(smiles)
        options = get_pg_options_for_site(sites[0]['site'])
        selection = create_protection_selection(sites[0], options[0])

        plan = create_protection_plan(smiles, [selection])
        pgs = plan.get_all_pgs()

        self.assertIsInstance(pgs, list)
        self.assertEqual(len(pgs), 1)
        self.assertEqual(pgs[0], options[0]['suggestion'].protecting_group_id)

    def test_protection_plan_format_summary(self):
        """Test ProtectionPlan.format_summary method."""
        smiles = "CCO"
        sites = get_sites_for_selection(smiles)
        options = get_pg_options_for_site(sites[0]['site'])
        selection = create_protection_selection(sites[0], options[0])

        plan = create_protection_plan(smiles, [selection])
        summary = plan.format_summary()

        self.assertIsInstance(summary, str)
        self.assertIn(smiles, summary)
        self.assertIn("Protection Plan", summary)

    def test_hitl_workflow_example(self):
        """Test the hitl_workflow_example function."""
        result = hitl_workflow_example("CCO")

        self.assertIsInstance(result, str)
        self.assertIn("HITL Workflow", result)
        self.assertIn("STEP 1", result)
        self.assertIn("STEP 2", result)

    def test_format_sites_for_display(self):
        """Test format_sites_for_display function."""
        sites = get_sites_for_selection("CCO")
        result = format_sites_for_display(sites)

        self.assertIsInstance(result, str)
        self.assertIn("Available Protection Sites", result)
        self.assertIn("[0]", result)
        self.assertIn("Atom map numbers", result)

    def test_format_sites_for_display_empty(self):
        """Test format_sites_for_display with empty list."""
        result = format_sites_for_display([])
        self.assertIn("No protection sites", result)

    def test_format_pg_options_for_display(self):
        """Test format_pg_options_for_display function."""
        sites = get_sites_for_selection("CCO")
        options = get_pg_options_for_site(sites[0]['site'])
        result = format_pg_options_for_display(options)

        self.assertIsInstance(result, str)
        self.assertIn("Available Protecting Groups", result)
        self.assertIn("[0]", result)

    def test_format_pg_options_for_display_empty(self):
        """Test format_pg_options_for_display with empty list."""
        result = format_pg_options_for_display([])
        self.assertIn("No protecting group options", result)

    def test_multiple_site_selection_workflow(self):
        """Test workflow with multiple site selections."""
        smiles = "NCCO"  # Has both amine and alcohol
        sites = get_sites_for_selection(smiles)

        selections = []
        selected_pgs = []

        for site_info in sites[:2]:
            options = get_pg_options_for_site(site_info['site'],
                                              existing_pgs=selected_pgs)
            if options:
                selection = create_protection_selection(site_info, options[0])
                selections.append(selection)
                selected_pgs.append(
                    options[0]['suggestion'].protecting_group_id)

        plan = create_protection_plan(smiles, selections)

        self.assertEqual(len(plan.selections), len(selections))
        self.assertEqual(len(plan.get_all_pgs()), len(selections))
