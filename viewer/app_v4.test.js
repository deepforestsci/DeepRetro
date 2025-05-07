/**
 * @jest-environment jsdom
 */

// Mock global objects first
global.d3 = {
  select: jest.fn().mockReturnThis(),
  selectAll: jest.fn().mockReturnThis(),
  append: jest.fn().mockReturnThis(),
  attr: jest.fn().mockReturnThis(),
  style: jest.fn().mockReturnThis(),
  text: jest.fn().mockReturnThis(),
  html: jest.fn().mockReturnThis(),
  data: jest.fn().mockReturnThis(),
  enter: jest.fn().mockReturnThis(),
  on: jest.fn().mockReturnThis(),
  call: jest.fn().mockReturnThis(),
  transition: jest.fn().mockReturnThis(),
  duration: jest.fn().mockReturnThis(),
  each: jest.fn().mockImplementation(function (callback) {
    if (callback) callback(); // Execute the callback for testing
    return this;
  }),
  hierarchy: jest.fn().mockReturnValue({
    descendants: jest.fn().mockReturnValue([]),
    links: jest.fn().mockReturnValue([]),
    each: jest.fn().mockImplementation((callback) => {
      if (callback) callback({ data: { step: "0", products: [] }, x: 0, y: 0 });
      return this;
    }),
  }),
  tree: jest.fn().mockReturnValue(function () {
    return {
      nodeSize: jest.fn().mockReturnThis(),
      separation: jest.fn().mockReturnThis(),
    };
  }),
  zoom: jest.fn().mockReturnValue({
    scaleExtent: jest.fn().mockReturnThis(),
    on: jest.fn().mockReturnThis(),
    scaleBy: jest.fn(),
    transform: jest.fn(),
  }),
  zoomIdentity: {},
};

// Mock OpenChemLib
global.OCL = {
  Molecule: {
    fromSmiles: jest.fn().mockReturnValue({
      toSVG: jest.fn().mockReturnValue("<svg><g></g></svg>"),
      getAllAtoms: jest.fn().mockReturnValue(1),
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

// Mock window.alert
global.alert = jest.fn();
global.window = {
  ...global.window,
  alert: jest.fn(),
};

// Mock File constructor
global.File = class File {
  constructor(content, name, options = {}) {
    this.name = name;
    this.size = content.length;
    this.type = options.type || "";
    this._content = content;
  }
};

// Import real implementations for testing and coverage
const app = require("./app_v4");
global.updatePathwayNumber = app.updatePathwayNumber;
global.handleFileSelect = app.handleFileSelect;
global.processData = app.processData;
global.calculateMoleculeSize = app.calculateMoleculeSize;
global.calculateStepSize = app.calculateStepSize;
global.formatFormula = app.formatFormula;
global.renderGraph = app.renderGraph;

describe("app_v4.js", () => {
  // Setup common test data
  let mockData;
  let mockStep;

  beforeEach(() => {
    // Reset mocks
    jest.clearAllMocks();

    // Reset global variables that might be modified by tests
    global.fileData = undefined;
    global.fileHandlerInitialized = false;

    // Setup mock DOM elements
    document.body.innerHTML = `
      <div id="graph"></div>
      <div id="pathway-number"><span id="current-pathway">-</span></div>
      <div id="file-result" style="display: none;"></div>
      <div id="file-json-toggle" style="display: none;"></div>
      <div id="file-json-arrow"></div>
      <input type="file" id="file-input" />
    `;

    // Create mock data structure
    mockData = {
      steps: [
        {
          step: 1,
          reactants: [
            {
              smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
              reactant_metadata: {
                chemical_formula: "C9H8O4",
                mass: 180.16,
              },
            },
          ],
          products: [
            {
              smiles: "CC(=O)O",
              product_metadata: {
                chemical_formula: "C2H4O2",
                mass: 60.05,
              },
            },
            {
              smiles: "OC1=CC=CC=C1C(=O)O",
              product_metadata: {
                chemical_formula: "C7H6O3",
                mass: 138.12,
              },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: 8.5,
              confidenceestimate: 0.94,
              closestliterature: "Patent XYZ",
            },
          ],
          conditions: {
            temperature: "25°C",
            pressure: "1 atm",
            solvent: "H2O",
            time: "2h",
          },
        },
        {
          step: 2,
          reactants: [
            {
              smiles: "OC1=CC=CC=C1C(=O)O",
              reactant_metadata: {
                chemical_formula: "C7H6O3",
                mass: 138.12,
              },
            },
          ],
          products: [
            {
              smiles: "OC1=CC=CC=C1",
              product_metadata: {
                chemical_formula: "C6H6O",
                mass: 94.11,
              },
            },
          ],
          parent_id: 1,
          reactionmetrics: [
            {
              scalabilityindex: 7.2,
              confidenceestimate: 0.88,
              closestliterature: "Journal ABC",
            },
          ],
          conditions: {
            temperature: "100°C",
            pressure: "1 atm",
            solvent: "Toluene",
            time: "4h",
          },
        },
      ],
      dependencies: {
        1: ["2"],
      },
    };

    // Mock step structure for renderGraph
    mockStep = {
      step: {
        step: "0",
        products: [
          {
            smiles: "CC(=O)O",
            product_metadata: {
              chemical_formula: "C2H4O2",
              mass: 60.05,
            },
          },
        ],
        reactionmetrics: [
          {
            scalabilityindex: 9.0,
            confidenceestimate: 0.96,
            closestliterature: "Reference ABC",
          },
        ],
        conditions: {
          temperature: "25°C",
          pressure: "1 atm",
          solvent: "H2O",
          time: "1h",
        },
      },
      children: {
        1: {
          step: {
            step: "1",
            reactants: [
              {
                smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
                reactant_metadata: {
                  chemical_formula: "C9H8O4",
                  mass: 180.16,
                },
              },
            ],
            reactionmetrics: [
              {
                scalabilityindex: 8.5,
                confidenceestimate: 0.94,
                closestliterature: "Patent XYZ",
              },
            ],
            conditions: {
              temperature: "25°C",
              pressure: "1 atm",
              solvent: "H2O",
              time: "2h",
            },
          },
          children: {},
        },
      },
    };
  });

  describe("updatePathwayNumber", () => {
    test("should update the current pathway number", () => {
      updatePathwayNumber(3);
      expect(document.getElementById("current-pathway").textContent).toBe("3");
      expect(document.getElementById("pathway-number").style.display).toBe(
        "block"
      );
    });

    test("should handle error when elements are not found", () => {
      // Remove the elements to simulate error
      document.body.innerHTML = "";

      // Mock console.error
      const consoleErrorSpy = jest.spyOn(console, "error").mockImplementation();

      updatePathwayNumber(5);

      expect(consoleErrorSpy).toHaveBeenCalledWith(
        "Pathway display elements not found"
      );

      consoleErrorSpy.mockRestore();
    });
  });

  describe("handleFileSelect", () => {
    test("should handle error when parsing invalid JSON", () => {
      // Mock FileReader to return invalid JSON
      global.FileReader = class {
        constructor() {
          this.onload = null;
        }
        readAsText() {
          setTimeout(() => {
            if (this.onload) {
              this.onload({ target: { result: "invalid json" } });
            }
          }, 0);
        }
      };

      // Mock alert
      global.alert.mockClear();

      // Create mock event
      const mockEvent = {
        target: {
          files: [
            new File(["invalid json"], "test.json", { type: "text/plain" }),
          ],
        },
      };

      // Call the function
      handleFileSelect(mockEvent);

      // Use setTimeout to give FileReader.onload a chance to run
      return new Promise((resolve) => setTimeout(resolve, 10)).then(() => {
        // Verify alert was called
        expect(global.alert).toHaveBeenCalledWith(
          expect.stringContaining("Error parsing JSON")
        );
      });
    });
  });

  describe("processData", () => {
    test("should process data and construct a tree", () => {
      const result = processData(mockData);

      // Verify we have a root node (Step 0)
      expect(result).toHaveProperty("0");

      // Verify Step 0 was created with correct data
      const step0 = result["0"].step;
      expect(step0.step).toBe("0");
      expect(step0.products).toHaveLength(1);
      expect(step0.products[0]).toEqual(mockData.steps[0].products[0]);

      // Verify Step 1 is a child of Step 0
      expect(result["0"].children).toHaveProperty("1");
      expect(result["0"].children["1"].step.step).toBe("1");
    });

    test("should handle empty input data", () => {
      // Test with null data
      let result = processData(null);
      expect(result).toEqual({});

      // Test with empty steps array
      result = processData({ steps: [] });
      expect(result).toEqual({});
    });
  });

  describe("calculateMoleculeSize", () => {
    test("should calculate correct size based on formula", () => {
      const metadata = {
        chemical_formula: "C6H12O6",
        mass: 180.16,
      };

      const result = calculateMoleculeSize(metadata);

      expect(result).toHaveProperty("radius");
      expect(result).toHaveProperty("svgSize");
      expect(result.radius).toBeGreaterThan(45); // Base radius
      expect(result.svgSize).toBeGreaterThan(result.radius); // SVG size should be larger than radius
    });

    test("should return default size for missing metadata", () => {
      const result = calculateMoleculeSize(null);

      expect(result).toEqual({
        radius: 35,
        svgSize: 60,
      });
    });
  });

  describe("formatFormula", () => {
    test("should format chemical formula with subscripts", () => {
      const result = formatFormula("C6H12O6");

      expect(result).toBe(
        'C<tspan baseline-shift="sub">6</tspan>H<tspan baseline-shift="sub">12</tspan>O<tspan baseline-shift="sub">6</tspan>'
      );
    });
  });
});
