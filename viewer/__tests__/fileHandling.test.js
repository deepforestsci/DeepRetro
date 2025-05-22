/**
 * @jest-environment jsdom
 */

// Import the functions to test
const { handleFileSelect } = require("../app_v4");

describe("File Handling Functions", () => {
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

    // Mock alert to prevent errors in test environment
    global.alert = jest.fn();
  });

  test("handleFileSelect processes valid file input", () => {
    const mockFile = new Blob(
      [
        JSON.stringify({
          steps: [{ step: "1", products: [{ smiles: "CCO" }] }],
        }),
      ],
      { type: "application/json" }
    );
    const mockEvent = { target: { files: [mockFile] } };

    const fileInput = document.createElement("input");
    fileInput.type = "file";
    document.body.appendChild(fileInput);

    const reader = new FileReader();
    jest.spyOn(window, "FileReader").mockImplementation(() => reader);
    jest.spyOn(reader, "readAsText").mockImplementation(function () {
      this.onload({
        target: {
          result: JSON.stringify({
            steps: [{ step: "1", products: [{ smiles: "CCO" }] }],
          }),
        },
      });
    });

    handleFileSelect(mockEvent);

    expect(console.log).toHaveBeenCalledWith("Updated pathway number to:", 1);
  });

  test("handles JSON parsing errors", () => {
    // Spy on alert to verify it's called
    jest.spyOn(window, "alert").mockImplementation(() => {});

    // Create an invalid JSON file
    const invalidJsonFile = new Blob(["This is not valid JSON"], {
      type: "text/plain",
    });

    const mockEvent = { target: { files: [invalidJsonFile] } };

    // Setup DOM element for the test
    const fileInput = document.createElement("input");
    fileInput.type = "file";
    document.body.appendChild(fileInput);

    // Mock FileReader
    const reader = new FileReader();
    jest.spyOn(window, "FileReader").mockImplementation(() => reader);
    jest.spyOn(reader, "readAsText").mockImplementation(function () {
      this.onload({ target: { result: "This is not valid JSON" } });
    });

    // Execute the function
    handleFileSelect(mockEvent);

    // Verify alert was called with error message
    expect(window.alert).toHaveBeenCalledWith(
      expect.stringContaining("Error parsing JSON")
    );
  });
});
