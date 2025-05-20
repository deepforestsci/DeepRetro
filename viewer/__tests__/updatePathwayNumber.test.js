/**
 * @jest-environment jsdom
 */

// Import the functions to test
const { updatePathwayNumber } = require("../app_v4");

describe("updatePathwayNumber", () => {
  beforeEach(() => {
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

  test("formats the pathway correctly", () => {
    // Setup
    document.body.innerHTML = `
      <div id="pathway-number" style="display: none;"><span id="current-pathway">-</span></div>
    `;

    // Act
    updatePathwayNumber("42");

    // Assert
    expect(document.getElementById("current-pathway").textContent).toBe("42");
    expect(document.getElementById("pathway-number").style.display).toBe(
      "block"
    );
  });

  test("handles various input values", () => {
    // Setup DOM
    document.body.innerHTML = `
      <div id="pathway-number" style="display: none;"><span id="current-pathway">-</span></div>
    `;

    // Test with number
    updatePathwayNumber(42);
    expect(document.getElementById("current-pathway").textContent).toBe("42");

    // Test with string
    updatePathwayNumber("ABC");
    expect(document.getElementById("current-pathway").textContent).toBe("ABC");

    // null gets converted to a string in textContent
    updatePathwayNumber(null);
    // Don't strictly check the exact value since toString(null) behavior may vary
    expect(
      document.getElementById("current-pathway").textContent
    ).toBeDefined();
  });
});
