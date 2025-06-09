/**
 * @jest-environment jsdom
 *
 * Test Suite for processData Function
 *
 * These tests verify the data processing functionality for chemical reaction pathways.
 * The processData function is responsible for:
 * 1. Converting raw reaction data into a hierarchical tree structure
 * 2. Adding Step 0 as a root node
 * 3. Establishing parent-child relationships between steps
 * 4. Handling various edge cases and data formats
 */

// Import the functions to test
const { processData } = require("../app_v4");

// Common test data
const SAMPLE_PRODUCT = {
  smiles: "CCO",
  product_metadata: { chemical_formula: "C2H6O" },
};

const SAMPLE_REACTANT = {
  smiles: "C",
  reactant_metadata: { chemical_formula: "CH4" },
};

const SAMPLE_CONDITIONS = {
  temperature: "25C",
  pressure: "1 atm",
  solvent: "water",
  time: "1h",
};

describe("processData", () => {
  beforeEach(() => {
    // Reset console methods before each test
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };
  });

  // Basic Functionality Tests
  describe("Basic Data Processing", () => {
    test("handles input with steps as an array", () => {
      const testData = {
        steps: [
          {
            step: "1",
            products: [SAMPLE_PRODUCT],
            reactants: [SAMPLE_REACTANT],
            conditions: SAMPLE_CONDITIONS,
            reactionmetrics: [{ scalabilityindex: "10" }],
          },
        ],
        dependencies: { 1: [] },
      };

      const result = processData(testData);
      expect(result).toBeDefined();
      expect(console.log).toHaveBeenCalled();
    });

    test("handles complex dependencies", () => {
      const complexData = {
        steps: [
          {
            step: "1",
            products: [SAMPLE_PRODUCT],
            reactants: [SAMPLE_REACTANT],
          },
          {
            step: "2",
            products: [
              {
                smiles: "C",
                product_metadata: { chemical_formula: "CH4" },
              },
            ],
            reactants: [
              {
                smiles: "H2",
                reactant_metadata: { chemical_formula: "H2" },
              },
            ],
          },
        ],
        dependencies: {
          1: ["2"],
          2: [],
        },
      };

      const result = processData(complexData);
      expect(result).toBeDefined();
      expect(Object.keys(result)).toContain("0");
    });
  });

  // Edge Cases and Error Handling
  describe("Edge Cases", () => {
    test("handles empty or missing steps", () => {
      const emptyData = { steps: [] };
      const result = processData(emptyData);
      expect(result).toEqual({});
      expect(console.warn).toHaveBeenCalled();
    });

    test("handles step 1 with no products", () => {
      const noProductsData = {
        steps: [
          {
            step: "1",
            reactants: [SAMPLE_REACTANT],
            conditions: {},
            reactionmetrics: {},
          },
        ],
      };

      const result = processData(noProductsData);
      expect(result).toBeDefined();
      expect(console.warn).toHaveBeenCalledWith(
        expect.stringContaining("Step 1 has no products defined")
      );
    });

    test("handles step with empty or missing reactants", () => {
      const data = {
        steps: [
          {
            step: "1",
            products: [SAMPLE_PRODUCT],
            conditions: SAMPLE_CONDITIONS,
            reactionmetrics: [{ scalabilityindex: "10" }],
          },
        ],
        dependencies: { 1: [] },
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(Object.keys(result).length).toBe(1);
      expect(result["0"]).toBeDefined();
    });
  });

  // Dependency Handling Tests
  describe("Dependency Handling", () => {
    test("handles complex dependency chains", () => {
      const data = {
        steps: [
          { step: "1", products: [{ smiles: "CCO" }] },
          { step: "2", products: [{ smiles: "CC" }] },
          { step: "3", products: [{ smiles: "C" }] },
          { step: "4", products: [{ smiles: "O" }] },
        ],
        dependencies: {
          1: ["2"],
          2: ["3", "4"],
          3: [],
          4: [],
        },
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(Object.keys(result).length).toBeGreaterThan(0);
      expect(result["0"].children["1"]).toBeDefined();
    });

    test("handles null dependencies", () => {
      const data = {
        steps: [{ step: "1", products: [{ smiles: "CCO" }] }],
        dependencies: null,
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(Object.keys(result).length).toBe(1);
    });

    test("handles empty dependencies", () => {
      const data = {
        steps: [{ step: "1", products: [{ smiles: "CCO" }] }],
        dependencies: {},
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(Object.keys(result).length).toBe(1);
    });

    test("with missing dependencies for a step", () => {
      const data = {
        steps: [
          { step: "1", products: [{ smiles: "CCO" }] },
          { step: "2", products: [{ smiles: "CC" }] },
          { step: "3", products: [{ smiles: "C" }] }, // No dependencies defined
        ],
        dependencies: {
          1: ["2"],
          2: [],
        },
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(result["0"].children["1"].children["2"]).toBeDefined();
    });
  });

  // Data Format Handling Tests
  describe("Data Format Handling", () => {
    test("handles missing steps key", () => {
      const data = {
        dependencies: { 1: [] },
      };

      const result = processData(data);
      expect(result).toEqual({});
    });

    test("handles unusual step numbering", () => {
      const data = {
        steps: [
          { step: "1", products: [{ smiles: "CCO" }] }, // String step number
          { step: 2, products: [{ smiles: "CC" }] }, // Numeric step number
        ],
        dependencies: {
          1: ["2"],
          2: [],
        },
      };

      const result = processData(data);
      expect(result).toBeDefined();
      expect(Object.keys(result).length).toBeGreaterThan(0);
      expect(result["0"]).toBeDefined();
    });

    test("parent_id assignment logic", () => {
      const data = {
        steps: [
          { step: "1", products: [{ smiles: "CCO" }], parent_id: 0 }, // Numeric
          { step: "2", products: [{ smiles: "CC" }], parent_id: "1" }, // String
          { step: "3", products: [{ smiles: "C" }], parent_id: null }, // Null
          { step: "4", products: [{ smiles: "O" }] }, // Missing
        ],
        dependencies: {},
      };

      const result = processData(data);
      expect(result).toBeDefined();
    });
  });
});
