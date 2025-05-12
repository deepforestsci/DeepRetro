/**
 * @jest-environment jsdom
 */

// Import the functions to test
const {
  updatePathwayNumber,
  processData,
  calculateMoleculeSize,
  calculateStepSize,
  formatFormula,
  renderGraph,
} = require("./app_v4");

// Mock D3 with improved mocking for SVG operations
const mockAppendReturn = {
  attr: jest.fn().mockReturnThis(),
  style: jest.fn().mockReturnThis(),
  html: jest.fn().mockReturnThis(),
  text: jest.fn().mockReturnThis(),
  on: jest.fn().mockReturnThis(),
  append: jest.fn().mockReturnThis(),
  selectAll: jest.fn().mockReturnThis(),
  data: jest.fn().mockReturnThis(),
  enter: jest.fn().mockReturnThis(),
  call: jest.fn().mockReturnThis(),
  each: jest.fn().mockReturnThis(),
  select: jest.fn().mockReturnThis(),
  remove: jest.fn().mockReturnThis(),
};

// Create mock hierarchy node data for testing
const mockHierarchyData = [
  {
    data: {
      step: "0",
      products: [
        {
          smiles: "CCO",
          product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
        },
      ],
      reactionmetrics: [
        {
          scalabilityindex: "10",
          confidenceestimate: 0.9,
          closestliterature: "",
        },
      ],
      conditions: {
        temperature: "25C",
        pressure: "1 atm",
        solvent: "water",
        time: "1h",
      },
    },
    x: 0,
    y: 0,
  },
  {
    data: {
      step: "1",
      reactants: [
        {
          smiles: "CC=O",
          reactant_metadata: { chemical_formula: "C2H4O", mass: 44.05 },
        },
      ],
      reactionmetrics: [
        {
          scalabilityindex: "8",
          confidenceestimate: 0.8,
          closestliterature: "",
        },
      ],
      conditions: {
        temperature: "25C",
        pressure: "1 atm",
        solvent: "water",
        time: "1h",
      },
    },
    x: 10,
    y: 10,
  },
];

// Create a more complete mock for hierarchyRoot
const mockHierarchyRoot = {
  x: 0,
  y: 0,
  data: {
    step: "0",
    products: [
      {
        smiles: "CCO",
        product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
      },
    ],
  },
  each: jest.fn().mockImplementation((callback) => {
    mockHierarchyData.forEach((node) => callback(node));
    return mockHierarchyRoot;
  }),
  links: jest.fn().mockReturnValue([]),
  descendants: jest.fn().mockReturnValue(mockHierarchyData),
};

global.d3 = {
  select: jest.fn().mockReturnValue({
    selectAll: jest.fn().mockReturnValue({
      remove: jest.fn(),
    }),
    append: jest.fn().mockReturnValue(mockAppendReturn),
    style: jest.fn().mockReturnThis(),
    attr: jest.fn().mockReturnThis(),
    call: jest.fn().mockReturnThis(),
  }),
  hierarchy: jest.fn().mockReturnValue(mockHierarchyRoot),
  tree: jest.fn().mockReturnValue({
    nodeSize: jest.fn().mockReturnValue({
      separation: jest.fn().mockReturnValue(jest.fn()),
    }),
  }),
  zoom: jest.fn().mockReturnValue({
    scaleExtent: jest.fn().mockReturnValue({
      on: jest.fn().mockReturnValue({}),
    }),
    transform: jest.fn(),
    scaleBy: jest.fn(),
  }),
  zoomIdentity: {},
};

// Mock OpenChemLib
global.OCL = {
  Molecule: {
    fromSmiles: jest.fn().mockReturnValue({
      toSVG: jest.fn().mockReturnValue("<svg><g></g></svg>"),
      getAllAtoms: jest.fn().mockReturnValue(10),
    }),
  },
};

// Mock DOMParser
global.DOMParser = class {
  parseFromString() {
    return {
      documentElement: {
        innerHTML: "<g></g>",
      },
      getElementsByTagName: jest.fn().mockReturnValue([]),
    };
  }
};

// Create a few mock elements
document.body.innerHTML = `
  <div id="graph"></div>
  <div id="pathway-number"><span id="current-pathway">-</span></div>
`;

describe("App_v4.js Tests", () => {
  beforeEach(() => {
    // Reset all mocks before each test
    jest.clearAllMocks();

    // Reset DOM elements for each test
    document.body.innerHTML = `
      <div id="graph"></div>
      <div id="pathway-number" style="display: none;"><span id="current-pathway">-</span></div>
    `;

    // Reset console methods
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };
  });

  // Test updatePathwayNumber function (line 7-27)
  describe("updatePathwayNumber", () => {
    test("updates pathway number correctly", () => {
      updatePathwayNumber("1");
      expect(document.getElementById("current-pathway").textContent).toBe("1");
      expect(document.getElementById("pathway-number").style.display).toBe(
        "block"
      );
    });

    test("handles missing elements gracefully", () => {
      document.body.innerHTML = `<div></div>`;
      updatePathwayNumber("2");
      expect(console.error).toHaveBeenCalled();
    });
  });

  // Test processData function (lines 39-204)
  describe("processData", () => {
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
  });

  // Test calculateMoleculeSize function (lines 206-229)
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
  });

  // Test calculateStepSize function (lines 231-242)
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
  });

  // Test formatFormula function (line 244-246)
  describe("formatFormula", () => {
    test("formats chemical formulas with subscripts", () => {
      const result = formatFormula("C6H12O6");
      expect(result).toContain("<tspan");
      expect(result).toContain("baseline-shift");
    });
  });

  // Test renderGraph function (lines 248-896)
  describe("renderGraph", () => {
    test("renders graph with root step", () => {
      const mockRootStep = {
        step: {
          step: "0",
          products: [
            {
              smiles: "CCO",
              product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: "10",
              confidenceestimate: 0.9,
              closestliterature: "",
            },
          ],
          conditions: {
            temperature: "25C",
            pressure: "1 atm",
            solvent: "water",
            time: "1h",
          },
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(d3.select).toHaveBeenCalledWith("#graph");
    });

    test("handles errors in molecule rendering", () => {
      // Mock OCL to throw an error
      global.OCL.Molecule.fromSmiles.mockImplementationOnce(() => {
        throw new Error("SMILES parsing error");
      });

      const mockRootStep = {
        step: {
          step: "0",
          products: [
            {
              smiles: "InvalidSMILES",
              product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: "10",
              confidenceestimate: 0.9,
              closestliterature: "",
            },
          ],
          conditions: {
            temperature: "25C",
            pressure: "1 atm",
            solvent: "water",
            time: "1h",
          },
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(console.error).toHaveBeenCalled();
    });

    test("handles invalid SMILES", () => {
      // Mock getAllAtoms to return 0 (invalid molecule)
      global.OCL.Molecule.fromSmiles.mockImplementationOnce(() => ({
        toSVG: jest.fn(),
        getAllAtoms: jest.fn().mockReturnValue(0),
      }));

      const mockRootStep = {
        step: {
          step: "0",
          products: [
            {
              smiles: "Invalid",
              product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: "10",
              confidenceestimate: 0.9,
              closestliterature: "",
            },
          ],
          conditions: {
            temperature: "25C",
            pressure: "1 atm",
            solvent: "water",
            time: "1h",
          },
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(console.error).toHaveBeenCalled();
    });

    test("handles failed SVG parsing", () => {
      // Mock DOMParser to return parser error
      global.DOMParser = class {
        parseFromString() {
          return {
            documentElement: { innerHTML: "" },
            getElementsByTagName: jest
              .fn()
              .mockImplementation((tag) =>
                tag === "parsererror"
                  ? [{ textContent: "SVG parsing error" }]
                  : []
              ),
          };
        }
      };

      const mockRootStep = {
        step: {
          step: "0",
          products: [
            {
              smiles: "CCO",
              product_metadata: { chemical_formula: "C2H6O", mass: 46.07 },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: "10",
              confidenceestimate: 0.9,
              closestliterature: "",
            },
          ],
          conditions: {
            temperature: "25C",
            pressure: "1 atm",
            solvent: "water",
            time: "1h",
          },
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(console.error).toHaveBeenCalled();
    });

    test("handles missing SMILES in molecules", () => {
      const mockRootStep = {
        step: {
          step: "0",
          products: [
            { product_metadata: { chemical_formula: "C2H6O", mass: 46.07 } },
          ], // No SMILES
          reactionmetrics: [
            {
              scalabilityindex: "10",
              confidenceestimate: 0.9,
              closestliterature: "",
            },
          ],
          conditions: {
            temperature: "25C",
            pressure: "1 atm",
            solvent: "water",
            time: "1h",
          },
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(console.warn).toHaveBeenCalled();
    });
  });
});
