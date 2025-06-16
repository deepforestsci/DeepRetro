/**
 * Test Setup Configuration
 *
 * This module provides mock implementations of external dependencies used in tests.
 * It includes mocks for:
 * - D3.js visualization library
 * - OpenChemLib (OCL) molecular structure library
 * - DOM Parser and manipulation
 * - Console and alert functions
 */

/**
 * Creates mock implementations of D3.js functionality
 * @returns {Object} Object containing all D3 mock objects
 */
function setupD3Mocks() {
  // Basic D3 selection mock that chains all methods
  const mockD3Element = {
    // Visualization methods
    append: jest.fn(() => mockD3Element),
    attr: jest.fn(() => mockD3Element),
    style: jest.fn(() => mockD3Element),
    html: jest.fn(() => mockD3Element),
    text: jest.fn(() => mockD3Element),

    // Selection methods
    selectAll: jest.fn(() => mockD3Element),
    select: jest.fn(() => mockD3Element),
    data: jest.fn(() => mockD3Element),
    enter: jest.fn(() => mockD3Element),

    // Event and animation methods
    on: jest.fn(() => mockD3Element),
    call: jest.fn(() => mockD3Element),
    transition: jest.fn(() => mockD3Element),
    duration: jest.fn(() => mockD3Element),
    remove: jest.fn(() => mockD3Element),

    // Data iteration
    each: jest.fn((callback) => {
      callback({ data: {}, x: 0, y: 0 });
      return mockD3Element;
    }),
  };

  // Hierarchy node mock with test data
  const mockHierarchyNode = {
    // Node position
    x: 0,
    y: 0,

    // Test molecule data
    data: {
      step: "0",
      products: [
        {
          product_metadata: {
            chemical_formula: "C6H12O6",
          },
        },
      ],
    },
    children: [],

    // Tree traversal methods
    each: jest.fn((callback) => {
      // Mock root node
      callback({
        x: 0,
        y: 0,
        data: { step: "0", products: [] },
      });
      // Mock child node
      callback({
        x: 10,
        y: 10,
        data: { step: "1", reactants: [] },
      });
      return mockHierarchyNode;
    }),

    // Tree structure methods
    links: jest.fn(() => []),
    descendants: jest.fn(() => [
      {
        data: {
          step: "0",
          products: [
            {
              smiles: "CCO",
              product_metadata: {
                chemical_formula: "C2H6O",
              },
            },
          ],
          reactionmetrics: [
            {
              scalabilityindex: "10",
            },
          ],
          conditions: {},
        },
        x: 0,
        y: 0,
      },
    ]),
  };

  // Tree layout mock that preserves node structure
  const mockTreeLayout = jest.fn((node) => node);
  mockTreeLayout.nodeSize = jest.fn(() => mockTreeLayout);
  mockTreeLayout.separation = jest.fn(() => mockTreeLayout);

  // Zoom behavior mock
  const mockZoom = {
    scaleExtent: jest.fn(() => mockZoom),
    on: jest.fn(() => mockZoom),
    transform: jest.fn(),
    scaleBy: jest.fn(),
  };

  // Global D3 object with all mock functionality
  global.d3 = {
    select: jest.fn(() => mockD3Element),
    hierarchy: jest.fn(() => mockHierarchyNode),
    tree: jest.fn(() => mockTreeLayout),
    zoom: jest.fn(() => mockZoom),
    zoomIdentity: {},
  };

  return {
    mockD3Element,
    mockHierarchyNode,
    mockTreeLayout,
    mockZoom,
  };
}

/**
 * Creates mock implementation of OpenChemLib functionality
 * Provides basic molecule parsing and SVG generation capabilities
 */
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

/**
 * Creates mock implementation of DOMParser
 * Returns a simple SVG structure for molecule rendering
 */
function setupDOMParserMocks() {
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
}

/**
 * Sets up DOM elements and console mocks
 * Creates necessary HTML structure and replaces console methods with Jest spies
 */
function setupDOMAndConsole() {
  // Create minimal DOM structure needed for tests
  document.body.innerHTML = `
    <div id="graph"></div>
    <div id="pathway-number" style="display: none;">
      <span id="current-pathway">-</span>
    </div>
  `;

  // Replace console methods with Jest spies
  global.console = {
    log: jest.fn(),
    error: jest.fn(),
    warn: jest.fn(),
  };

  // Replace alert with Jest spy
  global.alert = jest.fn();
}

// Export setup functions for use in tests
module.exports = {
  setupD3Mocks,
  setupOCLMocks,
  setupDOMParserMocks,
  setupDOMAndConsole,
};
