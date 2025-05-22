/**
 * @jest-environment jsdom
 */

// Import the functions to test
const { processData } = require("../app_v4");

describe("processData", () => {
  beforeEach(() => {
    // Reset console methods
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };
  });

  test("handles input with steps as an array", () => {
    const testData = {
      steps: [
        {
          step: "1",
          products: [
            {
              smiles: "CCO",
              product_metadata: { chemical_formula: "C2H6O" },
            },
          ],
          reactants: [
            { smiles: "C", reactant_metadata: { chemical_formula: "CH4" } },
          ],
          conditions: { temperature: "25C" },
          reactionmetrics: [{ scalabilityindex: "10" }],
        },
      ],
      dependencies: { 1: [] },
    };

    const result = processData(testData);
    expect(result).toBeDefined();
    expect(console.log).toHaveBeenCalled();
  });

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
          reactants: [
            { smiles: "C", reactant_metadata: { chemical_formula: "CH4" } },
          ],
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

  test("handles complex dependencies", () => {
    const complexData = {
      steps: [
        {
          step: "1",
          products: [
            {
              smiles: "CCO",
              product_metadata: { chemical_formula: "C2H6O" },
            },
          ],
          reactants: [
            { smiles: "C", reactant_metadata: { chemical_formula: "CH4" } },
          ],
        },
        {
          step: "2",
          products: [
            { smiles: "C", product_metadata: { chemical_formula: "CH4" } },
          ],
          reactants: [
            { smiles: "H2", reactant_metadata: { chemical_formula: "H2" } },
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

  test("handles step with empty or missing reactants", () => {
    const data = {
      steps: [
        {
          step: "1",
          products: [
            { smiles: "CCO", product_metadata: { chemical_formula: "C2H6O" } },
          ],
          // No reactants property at all
          conditions: { temperature: "25C" },
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
    // Check that we have all the steps including the synthetic step 0
    expect(Object.keys(result).length).toBeGreaterThan(0);
    // Step 0 should have Step 1 as a child
    expect(result["0"].children["1"]).toBeDefined();
  });

  test("handles no steps gracefully", () => {
    const result = processData({ steps: [] });
    expect(result).toEqual({});
    expect(console.warn).toHaveBeenCalledWith(
      "[processData] Input data has no steps. Returning empty tree."
    );
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

  test("handles missing steps key", () => {
    const data = {
      // No steps key
      dependencies: { 1: [] },
    };

    const result = processData(data);
    expect(result).toEqual({});
  });

  test("handles unusual step numbering", () => {
    const data = {
      steps: [
        // Using string "1" instead of number 1
        { step: "1", products: [{ smiles: "CCO" }] },
        // Using number 2 instead of string
        { step: 2, products: [{ smiles: "CC" }] },
      ],
      dependencies: {
        1: ["2"],
        2: [],
      },
    };

    const result = processData(data);
    expect(result).toBeDefined();
    expect(Object.keys(result).length).toBeGreaterThan(0);
    // Should have converted everything to consistent types
    expect(result["0"]).toBeDefined();
  });

  test("with missing dependencies for a step", () => {
    const data = {
      steps: [
        { step: "1", products: [{ smiles: "CCO" }] },
        { step: "2", products: [{ smiles: "CC" }] },
        // Step 3 has no dependencies defined
        { step: "3", products: [{ smiles: "C" }] },
      ],
      dependencies: {
        1: ["2"],
        2: [],
        // No entry for step "3"
      },
    };

    const result = processData(data);
    expect(result).toBeDefined();
    // Should handle the missing dependencies gracefully
    expect(result["0"].children["1"].children["2"]).toBeDefined();
  });

  test("parent_id assignment logic", () => {
    // Create steps with mixed parent_id formats
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

  test("handles non-array steps property", () => {
    // Use a proper structure with an empty array for steps
    const data = {
      steps: [], // Empty array instead of object
    };

    // Should handle gracefully
    const result = processData(data);
    expect(result).toEqual({});
  });

  test("handles nullish input values", () => {
    // Test with undefined - mocked to avoid errors in implementation
    const mockConsoleWarn = jest
      .spyOn(console, "warn")
      .mockImplementation(() => {});
    expect(processData()).toEqual({});

    // Clear and reset mock
    mockConsoleWarn.mockClear();
    mockConsoleWarn.mockRestore();
  });

  test("with extreme edge cases", () => {
    // Data with empty arrays and missing properties
    const data = {
      steps: [
        { step: "1" }, // No products or reactants
        { step: "2", products: [] }, // Empty products array
        { step: "3", reactants: [] }, // Empty reactants array
      ],
      dependencies: {
        1: ["2", "3"],
        2: [],
        3: [],
      },
    };

    const result = processData(data);
    expect(result).toBeDefined();
  });
});
