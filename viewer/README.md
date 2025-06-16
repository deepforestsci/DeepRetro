# DeepRetro Viewer

A web-based interface for visualizing retrosynthesis pathways generated with large language models.

## Getting Started

### Prerequisites

- Git
- Modern web browser (Chrome, Firefox, Safari)
- Text editor for configuration

### Installation

1. Clone the repository:

```bash
git clone https://github.com/your-username/deepretro-viewer.git
cd deepretro-viewer
```

2. Configure environment variables:
   - Create a `.env` file in the root directory
   - Add your API key:
   ```
   API_KEY=your_api_key_here
   ```

### Local Setup

1. Navigate to the viewer directory:

```bash
cd viewer
```

2. Open the application:

   - Mac: Press `Cmd + O` in your browser
   - Windows: Press `Ctrl + O` in your browser
   - Linux: Press `Ctrl + O` in your browser
   - Navigate to the cloned repository
   - Select `viewer/index.html`

   Alternatively, you can also:

   - Drag and drop the `index.html` file directly into your browser
   - Double-click the `index.html` file to open in your default browser

## Usage

### Smart Retrosynthesis

1. Navigate to the "Smart Retrosynthesis" tab
2. Enter the SMILES notation of your target molecule
3. Click "Generate" to analyze the retrosynthesis pathway
4. The pathway will be displayed as an interactive visualization

#### Partial Re-runs

- You can edit any step in the generated pathway
- Click "Partial Rerun Analysis" to regenerate the pathway from that point

### View Pathway

1. Navigate to the "View Pathway" tab
2. Upload previously generated pathway files (.json format)
3. View and analyze existing pathways
4. Use partial re-run functionality to modify specific steps

## Features

- Interactive pathway visualization
- SMILES input support
- JSON pathway editing
- Partial pathway re-generation
- File upload for existing pathways
- Advanced settings configuration

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
