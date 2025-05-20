/** @type {import('jest').Config} */
const config = {
  // The root directory where Jest should look for tests
  rootDir: ".",

  // The test environment to use
  testEnvironment: "jsdom",

  // The directories where Jest should look for tests
  testMatch: ["**/__tests__/**/*.test.js", "**/?(*.)+(spec|test).js"],

  // Automatically clear mock calls and instances between tests
  clearMocks: true,

  // Indicates whether the coverage information should be collected
  collectCoverage: true,

  // The directory where Jest should output its coverage files
  coverageDirectory: "coverage",

  // Coverage reporters
  coverageReporters: ["text", "html"],

  // Files to collect coverage from
  collectCoverageFrom: ["app*.js", "!**/node_modules/**"],

  // An array of regexp pattern strings used to skip coverage collection
  coveragePathIgnorePatterns: [
    "/node_modules/",
    "/coverage/",
    "/dist/",
    "/__tests__/",
    "/jest.config.js",
    "/jest.setup.js",
  ],

  // Setup files to run before each test
  setupFilesAfterEnv: ["<rootDir>/jest.setup.js"],

  // A list of paths to modules that run some code to configure testing framework
  // before each test
  setupFiles: [],

  // Transform files with babel-jest
  transform: {},

  // An array of regexp pattern strings that are matched against all test paths,
  // matched tests are skipped
  testPathIgnorePatterns: ["/node_modules/"],

  // Enable verbose output
  verbose: true,
};

module.exports = config;
