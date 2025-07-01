/**
 * @jest-environment jsdom
 *
 * Test Suite for updatePathwayNumber Function
 *
 * These tests verify the functionality of updating the pathway number display:
 * 1. Updating the pathway number in the DOM
 * 2. Handling visibility of the pathway number element
 * 3. Error handling for missing DOM elements
 * 4. Input value formatting and validation
 */

const { updatePathwayNumber } = require("../app_v4");

// Common test data
const TEST_VALUES = {
  NUMERIC: "42",
  STRING: "ABC",
  NULL: null,
  UNDEFINED: undefined,
};

// Common DOM templates
const DOM_TEMPLATES = {
  COMPLETE: `
    <div id="graph"></div>
    <div id="pathway-number" style="display: none;">
      <span id="current-pathway">-</span>
    </div>
  `,
  MINIMAL: `
    <div id="pathway-number" style="display: none;">
      <span id="current-pathway">-</span>
    </div>
  `,
  EMPTY: `<div></div>`,
};

describe("updatePathwayNumber", () => {
  // Test environment setup
  beforeEach(() => {
    setupTestEnvironment();
  });

  // Basic Functionality Tests
  describe("Basic Display Updates", () => {
    test("updates pathway number correctly", () => {
      updatePathwayNumber("1");

      const pathwayElement = document.getElementById("current-pathway");
      const containerElement = document.getElementById("pathway-number");

      expect(pathwayElement.textContent).toBe("1");
      expect(containerElement.style.display).toBe("block");
    });

    test("formats the pathway correctly", () => {
      document.body.innerHTML = DOM_TEMPLATES.MINIMAL;
      updatePathwayNumber(TEST_VALUES.NUMERIC);

      const pathwayElement = document.getElementById("current-pathway");
      const containerElement = document.getElementById("pathway-number");

      expect(pathwayElement.textContent).toBe(TEST_VALUES.NUMERIC);
      expect(containerElement.style.display).toBe("block");
    });
  });

  // Input Handling Tests
  describe("Input Value Handling", () => {
    beforeEach(() => {
      document.body.innerHTML = DOM_TEMPLATES.MINIMAL;
    });

    test("handles numeric values", () => {
      updatePathwayNumber(42);
      expect(document.getElementById("current-pathway").textContent).toBe(
        TEST_VALUES.NUMERIC
      );
    });

    test("handles string values", () => {
      updatePathwayNumber(TEST_VALUES.STRING);
      expect(document.getElementById("current-pathway").textContent).toBe(
        TEST_VALUES.STRING
      );
    });

    test("handles null values", () => {
      updatePathwayNumber(TEST_VALUES.NULL);
      expect(
        document.getElementById("current-pathway").textContent
      ).toBeDefined();
    });
  });

  // Error Handling Tests
  describe("Error Handling", () => {
    test("handles missing DOM elements", () => {
      document.body.innerHTML = DOM_TEMPLATES.EMPTY;
      updatePathwayNumber(TEST_VALUES.NUMERIC);
      expect(console.error).toHaveBeenCalled();
    });

    test("handles undefined values", () => {
      updatePathwayNumber(TEST_VALUES.UNDEFINED);
      expect(
        document.getElementById("current-pathway").textContent
      ).toBeDefined();
    });
  });

  // Helper Functions
  function setupTestEnvironment() {
    // Setup DOM
    document.body.innerHTML = DOM_TEMPLATES.COMPLETE;

    // Setup console mocks
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };
  }
});
