// Mock browser globals that might be missing in jest-environment-jsdom
global.URL.createObjectURL = jest.fn();
global.URL.revokeObjectURL = jest.fn();

// Mock fetch
global.fetch = jest.fn();

// Mock D3.js
global.d3 = {
  v6: {
    min: {},
  },
};

// Mock graph-related functions
global.processData = jest.fn().mockReturnValue({ 0: { children: [] } });
global.renderGraph = jest.fn();
global.updatePathwayNumber = jest.fn();

// Common test data
global.testData = {
  smiles: "CC(=O)Oc1ccccc1C(=O)O", // Aspirin
  mockPathway: {
    steps: {
      1: {
        reactant: "CC(=O)Oc1ccccc1C(=O)O",
        products: [{ smiles: "CC(=O)O" }],
        reaction_type: "Hydrolysis",
      },
    },
  },
};

// Mock configuration
global.config = {
  instances: [
    {
      url: "http://test-server-1",
      defaults: {
        model_type: "claude3",
        advanced_prompt: false,
        model_version: "USPTO",
        stability_flag: false,
        hallucination_check: false,
      },
    },
    {
      url: "http://test-server-2",
      defaults: {
        model_type: "claude37",
        advanced_prompt: true,
        model_version: "Pistachio_25",
        stability_flag: true,
        hallucination_check: false,
      },
    },
  ],
  endpoints: {
    retrosynthesis: "/api/retrosynthesis",
    rerun: "/api/rerun",
    partial_rerun: "/api/partial_rerun",
  },
};

// Reset all mocks before each test
beforeEach(() => {
  jest.clearAllMocks();
});
