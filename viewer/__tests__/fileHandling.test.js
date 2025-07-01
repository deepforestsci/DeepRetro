/**
 * @jest-environment jsdom
 *
 * Test Suite for File Handling Functions
 *
 * These tests verify the functionality of file handling operations:
 * 1. Processing uploaded JSON files
 * 2. Parsing file contents and updating the pathway
 * 3. Error handling for invalid files and JSON parsing
 * 4. DOM updates based on file contents
 */

const { handleFileSelect } = require("../app_v4");

// Common test data
const SAMPLE_DATA = {
  VALID_JSON: {
    steps: [
      {
        step: "1",
        products: [{ smiles: "CCO" }],
      },
    ],
  },
  INVALID_JSON: "This is not valid JSON",
};

// Common DOM templates
const DOM_TEMPLATES = {
  COMPLETE: `
    <div id="graph"></div>
    <div id="pathway-number" style="display: none;">
      <span id="current-pathway">-</span>
    </div>
  `,
};

describe("File Handling Functions", () => {
  // Test environment setup
  beforeEach(() => {
    setupTestEnvironment();
  });

  // Basic Functionality Tests
  describe("Valid File Processing", () => {
    test("processes valid JSON file input", () => {
      // Setup test data
      const { mockEvent, reader } = createMockFileEvent(
        SAMPLE_DATA.VALID_JSON,
        "application/json"
      );

      // Setup file reader mock
      setupFileReaderMock(reader, SAMPLE_DATA.VALID_JSON);

      // Execute file handling
      handleFileSelect(mockEvent);

      // Verify pathway update
      expect(console.log).toHaveBeenCalledWith("Updated pathway number to:", 1);
    });
  });

  // Error Handling Tests
  describe("Error Handling", () => {
    test("handles JSON parsing errors", () => {
      // Setup test data
      const { mockEvent, reader } = createMockFileEvent(
        SAMPLE_DATA.INVALID_JSON,
        "text/plain"
      );

      // Setup file reader mock
      setupFileReaderMock(reader, SAMPLE_DATA.INVALID_JSON);

      // Execute file handling
      handleFileSelect(mockEvent);

      // Verify error handling
      expect(window.alert).toHaveBeenCalledWith(
        expect.stringContaining("Error parsing JSON")
      );
    });

    test("handles missing file input", () => {
      const mockEvent = { target: {} };
      handleFileSelect(mockEvent);
      expect(console.error).toHaveBeenCalledWith(
        "[handleFileSelect] Invalid event or missing files list"
      );
    });

    test("handles null file list", () => {
      const mockEvent = { target: { files: null } };
      handleFileSelect(mockEvent);
      expect(console.error).toHaveBeenCalledWith(
        "[handleFileSelect] Invalid event or missing files list"
      );
    });

    test("handles empty file list", () => {
      const mockEvent = { target: { files: [] } };
      handleFileSelect(mockEvent);
      expect(console.error).toHaveBeenCalledWith(
        "[handleFileSelect] No file selected"
      );
    });
  });

  // Helper Functions
  function setupTestEnvironment() {
    // Reset all mocks
    jest.clearAllMocks();

    // Setup DOM
    document.body.innerHTML = DOM_TEMPLATES.COMPLETE;

    // Setup console mocks
    global.console = {
      log: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
    };

    // Setup alert mock
    global.alert = jest.fn();
  }

  function createMockFileEvent(content, type) {
    // Create mock file
    const mockFile = new Blob(
      [typeof content === "string" ? content : JSON.stringify(content)],
      { type }
    );

    // Create mock event
    const mockEvent = {
      target: {
        files: [mockFile],
      },
    };

    // Setup file input element
    const fileInput = document.createElement("input");
    fileInput.type = "file";
    document.body.appendChild(fileInput);

    // Setup FileReader
    const reader = new FileReader();
    jest.spyOn(window, "FileReader").mockImplementation(() => reader);

    return { mockEvent, reader };
  }

  function setupFileReaderMock(reader, content) {
    jest.spyOn(reader, "readAsText").mockImplementation(function () {
      this.onload({
        target: {
          result:
            typeof content === "string" ? content : JSON.stringify(content),
        },
      });
    });
  }
});
