# DeepRetro Viewer

A web-based interface for visualizing retrosynthesis pathways generated with large language models.

## Testing

The DeepRetro Viewer uses Jest with JSDOM for testing the UI components and functionality.

### Running Tests

To run the tests:

```bash
# Run tests with standard configuration
npm test

# Run tests without using watchman
npm run test:no-watchman
```

### Test Coverage

The test suite includes coverage for:

- View toggling between API and File views
- Advanced settings display
- JSON data display
- SMILES input handling
- JSON file upload handling
- JSON editing
- Partial rerun functionality

### Adding Tests

When adding new tests, follow these patterns:

1. Import necessary dependencies at the top of the test file
2. Group related tests using `describe` blocks
3. Set up test environment in `beforeEach` function
4. Use meaningful test descriptions that explain the expected behavior

### Test Setup

The `jest.setup.js` file contains common mocks and test data used across all tests. You can add additional test data or mocks to this file when needed.
