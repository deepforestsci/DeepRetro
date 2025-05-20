/**
 * Common setup for tests
 */

// Mock D3
function setupD3Mocks() {
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

  return {
    mockD3Element,
    mockHierarchyNode,
    mockTreeLayout,
    mockZoom,
  };
}

// Mock OCL
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

// Mock DOMParser
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

// Reset DOM and console
function setupDOMAndConsole() {
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

  // Mock alert
  global.alert = jest.fn();
}

module.exports = {
  setupD3Mocks,
  setupOCLMocks,
  setupDOMParserMocks,
  setupDOMAndConsole,
};
