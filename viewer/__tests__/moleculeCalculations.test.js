/**
 * @jest-environment jsdom
 *
 * Test Suite for Molecule Calculation Functions
 *
 * These tests verify the functionality of molecule-related calculations:
 * 1. calculateMoleculeSize: Determines display size based on chemical formula
 * 2. calculateStepSize: Finds the largest molecule in a reaction step
 * 3. formatFormula: Formats chemical formulas with subscript numbers
 */

const {
  calculateMoleculeSize,
  calculateStepSize,
  formatFormula,
} = require("../app_v4");

// Common test data
const SAMPLE_FORMULAS = {
  SMALL: "H2",
  MEDIUM: "C6H12O6",
  LARGE: "C60H120O60",
  COMPLEX: "C60H120O60N20P10S5",
  INVALID: "C-H-O-N",
  SPECIAL: "C(CH3)3",
  NON_STANDARD: "$$C6H12O6$$",
};

const SAMPLE_METADATA = {
  SMALL: { chemical_formula: SAMPLE_FORMULAS.SMALL },
  MEDIUM: { chemical_formula: SAMPLE_FORMULAS.MEDIUM },
  LARGE: { chemical_formula: SAMPLE_FORMULAS.LARGE },
  COMPLEX: { chemical_formula: SAMPLE_FORMULAS.COMPLEX },
};

const DEFAULT_SIZE = { radius: 35, svgSize: 60 };
const MIN_RADIUS = 45;

describe("Molecule Calculation Functions", () => {
  describe("calculateMoleculeSize", () => {
    // Basic Functionality Tests
    describe("Basic Size Calculations", () => {
      test("calculates size based on chemical formula", () => {
        const result = calculateMoleculeSize(SAMPLE_METADATA.MEDIUM);
        expect(result).toHaveProperty("radius");
        expect(result).toHaveProperty("svgSize");
      });

      test("handles undefined or missing metadata", () => {
        const result = calculateMoleculeSize(undefined);
        expect(result).toEqual(DEFAULT_SIZE);
      });

      test("scales correctly with molecule complexity", () => {
        const small = calculateMoleculeSize(SAMPLE_METADATA.SMALL);
        const medium = calculateMoleculeSize(SAMPLE_METADATA.MEDIUM);
        const large = calculateMoleculeSize(SAMPLE_METADATA.COMPLEX);

        // Verify size scaling
        expect(small.radius).toBeLessThan(medium.radius);
        expect(medium.radius).toBeLessThan(large.radius);
        expect(small.svgSize).toBeLessThan(medium.svgSize);
        expect(medium.svgSize).toBeLessThan(large.svgSize);
      });
    });

    // Edge Cases and Special Inputs
    describe("Edge Cases", () => {
      test("handles invalid inputs", () => {
        expect(calculateMoleculeSize(null)).toEqual(DEFAULT_SIZE);
        expect(calculateMoleculeSize({ chemical_formula: null })).toEqual(
          DEFAULT_SIZE
        );
      });

      test("handles special formula patterns", () => {
        const specialChars = calculateMoleculeSize({
          chemical_formula: SAMPLE_FORMULAS.INVALID,
        });
        const irregularFormat = calculateMoleculeSize({
          chemical_formula: SAMPLE_FORMULAS.SPECIAL,
        });
        const nonStandard = calculateMoleculeSize({
          chemical_formula: SAMPLE_FORMULAS.NON_STANDARD,
        });

        // Verify all return valid sizes
        expect(specialChars.radius).toBeGreaterThan(0);
        expect(irregularFormat.radius).toBeGreaterThan(0);
        expect(nonStandard.radius).toBeGreaterThan(0);
      });

      test("handles empty and null formulas", () => {
        expect(calculateMoleculeSize({ chemical_formula: "" })).toEqual(
          DEFAULT_SIZE
        );
        expect(calculateMoleculeSize({ chemical_formula: null })).toEqual(
          DEFAULT_SIZE
        );
        expect(calculateMoleculeSize({ chemical_formula: "H" }).radius).toBe(
          MIN_RADIUS
        );
        expect(
          calculateMoleculeSize({ chemical_formula: "CHON" }).radius
        ).toBeGreaterThan(DEFAULT_SIZE.radius);
      });
    });

    // Complex Molecule Tests
    describe("Complex Molecules", () => {
      test("handles large molecules", () => {
        const result = calculateMoleculeSize(SAMPLE_METADATA.LARGE);
        expect(result.radius).toBeGreaterThan(MIN_RADIUS);

        const complexResult = calculateMoleculeSize(SAMPLE_METADATA.COMPLEX);
        expect(complexResult.svgSize).toBeGreaterThan(result.svgSize);
      });

      test("with different atom types", () => {
        const singleAtom = calculateMoleculeSize({ chemical_formula: "H10" });
        const multipleAtoms = calculateMoleculeSize({
          chemical_formula: "C10H10",
        });
        const manyAtoms = calculateMoleculeSize({
          chemical_formula: "C10H20N5O10",
        });

        expect(multipleAtoms.radius).toBeGreaterThanOrEqual(singleAtom.radius);
        expect(manyAtoms.radius).toBeGreaterThanOrEqual(multipleAtoms.radius);
      });
    });
  });

  describe("calculateStepSize", () => {
    // Test Data
    const SAMPLE_MOLECULES = {
      SMALL: {
        type: "reactant",
        reactant_metadata: { chemical_formula: "H2" },
      },
      MEDIUM: {
        type: "reactant",
        reactant_metadata: { chemical_formula: "C6H12O6" },
      },
      LARGE: {
        type: "step0",
        product_metadata: { chemical_formula: "C60H120O60" },
      },
    };

    // Basic Functionality Tests
    describe("Basic Size Calculations", () => {
      test("finds largest molecule in step", () => {
        const molecules = [SAMPLE_MOLECULES.SMALL, SAMPLE_MOLECULES.MEDIUM];

        const result = calculateStepSize(molecules);
        expect(result).toBeGreaterThan(0);
      });

      test("handles step0 type molecules", () => {
        const molecules = [SAMPLE_MOLECULES.LARGE];
        const result = calculateStepSize(molecules);
        expect(result).toBeGreaterThan(MIN_RADIUS);
      });
    });

    // Edge Cases
    describe("Edge Cases", () => {
      test("handles mixed molecule types", () => {
        const molecules = [
          SAMPLE_MOLECULES.SMALL,
          SAMPLE_MOLECULES.LARGE,
          SAMPLE_MOLECULES.MEDIUM,
        ];

        const size = calculateStepSize(molecules);
        expect(size).toBeGreaterThan(MIN_RADIUS);
      });

      test("handles empty and invalid inputs", () => {
        expect(typeof calculateStepSize([])).toBe("number");
        expect(typeof calculateStepSize([{}, {}])).toBe("number");
      });
    });
  });

  describe("formatFormula", () => {
    // Test Data
    const SAMPLE_FORMULAS_WITH_NUMBERS = {
      SIMPLE: "H2O",
      MEDIUM: "C6H12O6",
      COMPLEX: "C12H22O11N5",
    };

    // Basic Functionality Tests
    describe("Basic Formatting", () => {
      test("formats chemical formulas with subscripts", () => {
        const result = formatFormula(SAMPLE_FORMULAS_WITH_NUMBERS.MEDIUM);
        expect(result).toContain("<tspan");
        expect(result).toContain("baseline-shift");
      });

      test("handles complex formulas", () => {
        const result = formatFormula(SAMPLE_FORMULAS_WITH_NUMBERS.COMPLEX);
        expect(result.match(/<tspan/g).length).toBe(4); // One for each number: 12, 22, 11, 5
      });
    });

    // Edge Cases
    describe("Edge Cases", () => {
      test("handles various input types", () => {
        expect(formatFormula(SAMPLE_FORMULAS_WITH_NUMBERS.SIMPLE)).toContain(
          "<tspan"
        );
        expect(formatFormula("CHNO")).not.toContain("<tspan");
        expect(formatFormula(SAMPLE_FORMULAS_WITH_NUMBERS.COMPLEX)).toContain(
          "<tspan"
        );
      });

      test("handles special cases", () => {
        expect(formatFormula("")).toBe("");
        expect(formatFormula("CHO")).toBe("CHO");
        expect(formatFormula("123")).toContain("<tspan");
        expect(formatFormula("C123H456O789")).toContain("<tspan");
        expect(formatFormula("C-12-H-24")).toContain("<tspan");
      });
    });
  });
});
