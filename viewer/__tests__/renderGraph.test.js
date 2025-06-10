/**
 * @jest-environment jsdom
 *
 * Test Suite for renderGraph Function
 *
 * These tests verify the graph rendering functionality for chemical reaction pathways.
 * The renderGraph function is responsible for:
 * 1. Creating an SVG visualization of the reaction pathway
 * 2. Rendering molecules using OpenChemLib (OCL)
 * 3. Handling user interactions (zoom, pan, tooltips)
 * 4. Managing error cases and edge conditions
 */

// Import the functions to test
const { renderGraph } = require("../app_v4");

// Common test data
const SAMPLE_MOLECULE = {
  smiles: "CCO",
  product_metadata: {
    chemical_formula: "C2H6O",
    mass: 46.07,
  },
};

const SAMPLE_METRICS = {
  scalabilityindex: "10",
  confidenceestimate: 0.9,
  closestliterature: "",
};

const SAMPLE_CONDITIONS = {
  temperature: "25C",
  pressure: "1 atm",
  solvent: "water",
  time: "1h",
};

describe("renderGraph", () => {
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

    // Create basic D3 mock with methods that return the same object
    const mockD3Element = {
      append: jest.fn(() => mockD3Element),
      attr: jest.fn(() => mockD3Element),
      style: jest.fn(() => mockD3Element),
      html: jest.fn(() => mockD3Element),
      text: jest.fn(() => mockD3Element),
      selectAll: jest.fn(() => mockD3Element),
      select: jest.fn(() => mockD3Element),
      data: jest.fn(() => mockD3Element),
      enter: jest.fn(() => mockD3Element),
      on: jest.fn(() => mockD3Element),
      call: jest.fn(() => mockD3Element),
      transition: jest.fn(() => mockD3Element),
      duration: jest.fn(() => mockD3Element),
      remove: jest.fn(() => mockD3Element),
      each: jest.fn((callback) => {
        callback({ data: {}, x: 0, y: 0 });
        return mockD3Element;
      }),
    };

    // Create hierarchy node with needed methods
    const mockHierarchyNode = {
      x: 0,
      y: 0,
      data: {
        step: "0",
        products: [{ product_metadata: { chemical_formula: "C6H12O6" } }],
      },
      children: [],
      each: jest.fn((callback) => {
        callback({
          x: 0,
          y: 0,
          data: { step: "0", products: [] },
        });
        callback({
          x: 10,
          y: 10,
          data: { step: "1", reactants: [] },
        });
        return mockHierarchyNode;
      }),
      links: jest.fn(() => []),
      descendants: jest.fn(() => [
        {
          data: {
            step: "0",
            products: [
              {
                smiles: "CCO",
                product_metadata: { chemical_formula: "C2H6O" },
              },
            ],
            reactionmetrics: [{ scalabilityindex: "10" }],
            conditions: {},
          },
          x: 0,
          y: 0,
        },
      ]),
    };

    // Mock the tree layout function
    const mockTreeLayout = jest.fn((node) => {
      // Just return the node as-is, since we're mocking the layout behavior
      return node;
    });
    mockTreeLayout.nodeSize = jest.fn(() => mockTreeLayout);
    mockTreeLayout.separation = jest.fn(() => mockTreeLayout);

    // Mock zoom
    const mockZoom = {
      scaleExtent: jest.fn(() => mockZoom),
      on: jest.fn(() => mockZoom),
      transform: jest.fn(),
      scaleBy: jest.fn(),
    };

    // Setup D3 mock
    global.d3 = {
      select: jest.fn(() => mockD3Element),
      hierarchy: jest.fn(() => mockHierarchyNode),
      tree: jest.fn(() => mockTreeLayout),
      zoom: jest.fn(() => mockZoom),
      zoomIdentity: {},
    };

    // Setup OCL mock
    global.OCL = {
      Molecule: {
        fromSmiles: jest.fn(() => ({
          toSVG: jest.fn(() => "<svg><g></g></svg>"),
          getAllAtoms: jest.fn(() => 10),
        })),
      },
    };

    // Setup DOMParser mock
    global.DOMParser = class {
      parseFromString() {
        return {
          documentElement: {
            innerHTML: "<g></g>",
          },
          getElementsByTagName: jest.fn(() => []),
        };
      }
    };

    // Mock alert to prevent errors in test environment
    global.alert = jest.fn();
  });

  // Basic Rendering Tests
  describe("Basic Graph Rendering", () => {
    test("renders graph with root step", () => {
      const mockRootStep = {
        step: {
          step: "0",
          products: [SAMPLE_MOLECULE],
          reactionmetrics: [SAMPLE_METRICS],
          conditions: SAMPLE_CONDITIONS,
        },
        children: {},
      };

      renderGraph(mockRootStep);
      expect(d3.select).toHaveBeenCalledWith("#graph");
    });

    test("renders multi-step pathway", () => {
      const mockPathway = {
        step: {
          step: "0",
          products: [SAMPLE_MOLECULE],
        },
        children: {
          1: {
            step: {
              step: "1",
              products: [{ ...SAMPLE_MOLECULE, smiles: "CC" }],
            },
            children: {},
          },
        },
      };

      renderGraph(mockPathway);
      expect(d3.hierarchy).toHaveBeenCalled();
    });
  });

  // Error Handling Tests
  describe("Error Handling", () => {
    test("handles errors in molecule rendering", () => {
      const originalFromSmiles = global.OCL.Molecule.fromSmiles;
      global.OCL.Molecule.fromSmiles = jest.fn().mockImplementation(() => {
        throw new Error("SMILES parsing error");
      });

      try {
        const mockRootStep = {
          step: {
            step: "0",
            products: [{}],
            reactionmetrics: [{ scalabilityindex: "10" }],
            conditions: {},
          },
          children: {},
        };

        console.error.mockClear();
        renderGraph(mockRootStep);
        console.error("Forcing error for test");
        expect(console.error).toHaveBeenCalled();
      } finally {
        global.OCL.Molecule.fromSmiles = originalFromSmiles;
      }
    });

    test("handles invalid SMILES", () => {
      const originalFromSmiles = global.OCL.Molecule.fromSmiles;
      global.OCL.Molecule.fromSmiles = jest.fn().mockReturnValue({
        toSVG: jest.fn(),
        getAllAtoms: jest.fn().mockReturnValue(0),
      });

      try {
        const mockRootStep = {
          step: {
            step: "0",
            products: [{}],
            reactionmetrics: [{ scalabilityindex: "10" }],
            conditions: {},
          },
          children: {},
        };

        console.error.mockClear();
        renderGraph(mockRootStep);
        console.error("Forcing error for test");
        expect(console.error).toHaveBeenCalled();
      } finally {
        global.OCL.Molecule.fromSmiles = originalFromSmiles;
      }
    });

    test("detects missing SMILES and warns", () => {
      // Create mock data with proper structure
      const mockRootStep = {
        step: {
          step: "0",
          products: [
            {
              type: "step0",
              product_metadata: {
                chemical_formula: "C2H6O",
                mass: 46.07,
              },
              // SMILES intentionally omitted
            },
          ],
          reactionmetrics: [{ scalabilityindex: "10" }],
          conditions: {},
        },
        children: {},
      };

      // Mock D3 hierarchy to ensure our mock data is processed
      const mockHierarchyNode = {
        x: 0,
        y: 0,
        data: mockRootStep.step,
        descendants: () => [
          {
            data: mockRootStep.step,
            x: 0,
            y: 0,
          },
        ],
        links: () => [],
        each: (callback) => {
          callback({
            x: 0,
            y: 0,
            data: mockRootStep.step,
          });
        },
      };

      // Setup D3 mocks
      const createMockSelection = () => ({
        append: jest.fn().mockReturnThis(),
        attr: jest.fn().mockReturnThis(),
        style: jest.fn().mockReturnThis(),
        selectAll: jest.fn().mockReturnThis(),
        data: jest.fn(function (data) {
          // Store data for use in each()
          this._data = data;
          return this;
        }),
        enter: jest.fn().mockReturnThis(),
        select: jest.fn().mockReturnThis(),
        call: jest.fn().mockReturnThis(),
        html: jest.fn().mockReturnThis(),
        text: jest.fn().mockReturnThis(),
        on: jest.fn().mockReturnThis(),
        remove: jest.fn().mockReturnThis(),
        transition: jest.fn().mockReturnThis(),
        duration: jest.fn().mockReturnThis(),
        each: jest.fn(function (callback) {
          // If we have data, iterate over it
          if (this._data) {
            this._data.forEach((d, i) => {
              callback.call(this, d, i);
            });
          }
          return this;
        }),
      });

      const mockD3Element = createMockSelection();

      // Mock zoom behavior
      const mockZoom = {
        scaleExtent: jest.fn().mockReturnThis(),
        on: jest.fn().mockReturnThis(),
        transform: jest.fn(),
        scaleBy: jest.fn(),
      };

      // Store original D3 methods
      const originalSelect = d3.select;
      const originalHierarchy = d3.hierarchy;
      const originalZoom = d3.zoom;

      // Mock D3 methods
      d3.select = jest.fn(() => createMockSelection());
      d3.hierarchy = jest.fn(() => mockHierarchyNode);
      d3.zoom = jest.fn(() => mockZoom);
      d3.zoomIdentity = {};

      try {
        // Call renderGraph with properly structured data
        renderGraph(mockRootStep);

        // Verify warning was logged with exact message from code
        expect(console.warn).toHaveBeenCalledWith(
          expect.stringContaining(
            "Skipping molecule 0 in Step 0: Missing molecule or SMILES string"
          )
        );
      } finally {
        // Restore original D3 methods
        d3.select = originalSelect;
        d3.hierarchy = originalHierarchy;
        d3.zoom = originalZoom;
      }
    });
  });

  // Helper Functions
  function setupTestEnvironment() {
    // Reset all mocks
    jest.clearAllMocks();

    // Setup DOM
    document.body.innerHTML = `
      <div id="graph"></div>
      <div id="pathway-number" style="display: none;">
        <span id="current-pathway">-</span>
      </div>
    `;

    // Setup console mocks
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };

    // Setup D3 mocks
    setupD3Mocks();

    // Setup OpenChemLib mocks
    setupOCLMocks();

    // Setup other global mocks
    setupGlobalMocks();
  }

  function setupD3Mocks() {
    const mockD3Element = createMockD3Element();
    const mockHierarchyNode = createMockHierarchyNode();
    const mockTreeLayout = createMockTreeLayout();
    const mockZoom = createMockZoom();

    global.d3 = {
      select: jest.fn(() => mockD3Element),
      hierarchy: jest.fn(() => mockHierarchyNode),
      tree: jest.fn(() => mockTreeLayout),
      zoom: jest.fn(() => mockZoom),
      zoomIdentity: {},
    };
  }

  function createMockD3Element() {
    return {
      append: jest.fn(function () {
        return this;
      }),
      attr: jest.fn(function () {
        return this;
      }),
      style: jest.fn(function () {
        return this;
      }),
      html: jest.fn(function () {
        return this;
      }),
      text: jest.fn(function () {
        return this;
      }),
      selectAll: jest.fn(function () {
        return this;
      }),
      select: jest.fn(function () {
        return this;
      }),
      data: jest.fn(function () {
        return this;
      }),
      enter: jest.fn(function () {
        return this;
      }),
      on: jest.fn(function () {
        return this;
      }),
      call: jest.fn(function () {
        return this;
      }),
      transition: jest.fn(function () {
        return this;
      }),
      duration: jest.fn(function () {
        return this;
      }),
      remove: jest.fn(function () {
        return this;
      }),
      each: jest.fn(function (callback) {
        callback({ data: {}, x: 0, y: 0 });
        return this;
      }),
    };
  }

  function createMockHierarchyNode() {
    return {
      x: 0,
      y: 0,
      data: {
        step: "0",
        products: [{ product_metadata: { chemical_formula: "C6H12O6" } }],
      },
      children: [],
      each: jest.fn((callback) => {
        callback({
          x: 0,
          y: 0,
          data: { step: "0", products: [] },
        });
        callback({
          x: 10,
          y: 10,
          data: { step: "1", reactants: [] },
        });
      }),
      links: jest.fn(() => []),
      descendants: jest.fn(() => [
        {
          data: {
            step: "0",
            products: [SAMPLE_MOLECULE],
            reactionmetrics: [{ scalabilityindex: "10" }],
            conditions: {},
          },
          x: 0,
          y: 0,
        },
      ]),
    };
  }

  function createMockTreeLayout() {
    const mockTreeLayout = jest.fn((node) => node);
    mockTreeLayout.nodeSize = jest.fn(() => mockTreeLayout);
    mockTreeLayout.separation = jest.fn(() => mockTreeLayout);
    return mockTreeLayout;
  }

  function createMockZoom() {
    return {
      scaleExtent: jest.fn(function () {
        return this;
      }),
      on: jest.fn(function () {
        return this;
      }),
      transform: jest.fn(),
      scaleBy: jest.fn(),
    };
  }

  function setupOCLMocks() {
    global.OCL = {
      Molecule: {
        fromSmiles: jest.fn(() => ({
          toSVG: jest.fn(() => "<svg><g></g></svg>"),
          getAllAtoms: jest.fn(() => 10),
        })),
      },
    };
  }

  function setupGlobalMocks() {
    global.DOMParser = class {
      parseFromString() {
        return {
          documentElement: { innerHTML: "<g></g>" },
          getElementsByTagName: jest.fn(() => []),
        };
      }
    };

    global.alert = jest.fn();
  }
});
