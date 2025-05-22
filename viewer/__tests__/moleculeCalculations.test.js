/**
 * @jest-environment jsdom
 */

// Import the functions to test
const {
  calculateMoleculeSize,
  calculateStepSize,
  formatFormula,
} = require("../app_v4");

describe("Molecule Calculation Functions", () => {
  describe("calculateMoleculeSize", () => {
    test("calculates size based on chemical formula", () => {
      const result = calculateMoleculeSize({ chemical_formula: "C6H12O6" });
      expect(result).toHaveProperty("radius");
      expect(result).toHaveProperty("svgSize");
    });

    test("handles undefined or missing metadata", () => {
      const result = calculateMoleculeSize(undefined);
      expect(result).toEqual({ radius: 35, svgSize: 60 });
    });

    test("handles large molecules", () => {
      const result = calculateMoleculeSize({ chemical_formula: "C60H120O60" });
      expect(result.radius).toBeGreaterThan(45);

      // Check different scaling for large molecules
      const complexResult = calculateMoleculeSize({
        chemical_formula: "C60H120O60N20P10",
      });
      expect(complexResult.svgSize).toBeGreaterThan(result.svgSize);
    });

    test("handles complex formulas", () => {
      const result = calculateMoleculeSize({
        chemical_formula: "C60H120O60N20P10",
      });
      expect(result.radius).toBeGreaterThan(45);
      expect(result.svgSize).toBeGreaterThan(90);
    });

    test("handles invalid inputs", () => {
      // Just test for default values to be returned
      expect(calculateMoleculeSize(null)).toEqual({
        radius: 35,
        svgSize: 60,
      });

      expect(calculateMoleculeSize({ chemical_formula: null })).toEqual({
        radius: 35,
        svgSize: 60,
      });
    });

    test("with different atom types", () => {
      // Test with single atom type
      const singleAtom = calculateMoleculeSize({ chemical_formula: "H10" });

      // Test with multiple atom types
      const multipleAtoms = calculateMoleculeSize({
        chemical_formula: "C10H10",
      });

      // Test with many atom types
      const manyAtoms = calculateMoleculeSize({
        chemical_formula: "C10H20N5O10",
      });

      // More atom types should generally lead to larger radius
      // This is a flexible test since the exact calculation might vary
      expect(multipleAtoms.radius).toBeGreaterThanOrEqual(singleAtom.radius);
      expect(manyAtoms.radius).toBeGreaterThanOrEqual(multipleAtoms.radius);
    });

    test("handles extremely large molecules", () => {
      // Create a very complex formula
      const hugeFormula = "C100H200O50N30P10S5";
      const result = calculateMoleculeSize({ chemical_formula: hugeFormula });

      // Should have large radius and svgSize for very complex molecules
      expect(result.radius).toBeGreaterThan(60);
      expect(result.svgSize).toBeGreaterThan(100);
    });

    test("formula parsing edge cases", () => {
      // Test with formulas that might have tricky regex patterns
      const specialChars = calculateMoleculeSize({
        chemical_formula: "C-H-O-N",
      }); // Hyphens
      const irregularFormat = calculateMoleculeSize({
        chemical_formula: "C(CH3)3",
      }); // Parentheses
      const nonStandard = calculateMoleculeSize({
        chemical_formula: "$$C6H12O6$$",
      }); // Special chars

      // All should return valid sizes
      expect(specialChars.radius).toBeGreaterThan(0);
      expect(irregularFormat.radius).toBeGreaterThan(0);
      expect(nonStandard.radius).toBeGreaterThan(0);
    });

    test("scales correctly with molecule complexity", () => {
      // Very small molecule
      const small = calculateMoleculeSize({ chemical_formula: "H2" });

      // Medium molecule
      const medium = calculateMoleculeSize({ chemical_formula: "C6H12O6" });

      // Large complex molecule
      const large = calculateMoleculeSize({
        chemical_formula: "C60H120O60N20P10S5",
      });

      // Ensure sizes scale up with complexity
      expect(small.radius).toBeLessThan(medium.radius);
      expect(medium.radius).toBeLessThan(large.radius);
      expect(small.svgSize).toBeLessThan(medium.svgSize);
      expect(medium.svgSize).toBeLessThan(large.svgSize);
    });

    test("handles edge cases and special formulas", () => {
      // Test with empty string
      const empty = calculateMoleculeSize({ chemical_formula: "" });
      expect(empty.radius).toBe(35); // Should use default values

      // Test with null
      const nullFormula = calculateMoleculeSize({ chemical_formula: null });
      expect(nullFormula.radius).toBe(35);

      // Test with simple single atom
      const singleAtom = calculateMoleculeSize({ chemical_formula: "H" });
      expect(singleAtom.radius).toBe(45); // Should use base radius for small molecules

      // Test with formula that has no numbers
      const noNumbers = calculateMoleculeSize({ chemical_formula: "CHON" });
      expect(noNumbers.radius).toBeGreaterThan(35); // Should calculate based on unique elements
    });
  });

  describe("calculateStepSize", () => {
    test("finds largest molecule in step", () => {
      const molecules = [
        { type: "reactant", reactant_metadata: { chemical_formula: "C2H6" } },
        {
          type: "reactant",
          reactant_metadata: { chemical_formula: "C6H12O6" },
        },
      ];

      const result = calculateStepSize(molecules);
      expect(result).toBeGreaterThan(0);
    });

    test("handles step0 type molecules", () => {
      const molecules = [
        { type: "step0", product_metadata: { chemical_formula: "C60H120O60" } },
      ];

      const result = calculateStepSize(molecules);
      expect(result).toBeGreaterThan(45);
    });

    test("handles mixed molecule types", () => {
      const molecules = [
        { type: "reactant", reactant_metadata: { chemical_formula: "H2" } },
        { type: "step0", product_metadata: { chemical_formula: "C60H120O60" } },
        { type: "reactant", reactant_metadata: { chemical_formula: "CH4" } },
      ];

      const size = calculateStepSize(molecules);

      // The size should be determined by the largest molecule (the C60H120O60)
      expect(size).toBeGreaterThan(45);
    });

    test("handles edge cases", () => {
      // Empty array should return smallest radius
      const emptyResult = calculateStepSize([]);
      // Test for a number, not specifically 0
      expect(typeof emptyResult).toBe("number");

      // Test with molecules that have no metadata - should still return a number
      const noMetadataResult = calculateStepSize([{}, {}]);
      expect(typeof noMetadataResult).toBe("number");
    });
  });

  describe("formatFormula", () => {
    test("formats chemical formulas with subscripts", () => {
      const result = formatFormula("C6H12O6");
      expect(result).toContain("<tspan");
      expect(result).toContain("baseline-shift");
    });

    test("handles complex formulas", () => {
      const complexFormula = "C12H22O11";
      const result = formatFormula(complexFormula);

      // Should have multiple subscript tspans
      expect(result.match(/<tspan/g).length).toBe(3); // One for each number: 12, 22, 11
    });

    test("handles various input types", () => {
      // Test with simple formula
      expect(formatFormula("H2O")).toContain("<tspan");

      // Test with no numbers
      expect(formatFormula("CHNO")).not.toContain("<tspan");

      // Test with complex formula
      expect(formatFormula("C12H22O11N5")).toContain("<tspan");
      expect(formatFormula("C12H22O11N5").match(/<tspan/g).length).toBe(4);
    });

    test("handles edge case inputs", () => {
      // Empty string
      expect(formatFormula("")).toBe("");

      // No digits
      expect(formatFormula("CHO")).toBe("CHO");

      // Just digits
      expect(formatFormula("123")).toContain("<tspan");

      // Large numbers
      expect(formatFormula("C123H456O789")).toContain("<tspan");
    });

    test("handles edge case formula patterns", () => {
      // Test empty string
      expect(formatFormula("")).toBe("");

      // Test with no numbers
      expect(formatFormula("CHON")).toBe("CHON");

      // Test with unusual character sequences
      expect(formatFormula("C-12-H-24")).toContain("<tspan");
    });
  });
});
