// app.js

// Wait for the DOM to be fully loaded
document.addEventListener('DOMContentLoaded', () => {
	const fileInput = document.getElementById('fileInput');

	fileInput.addEventListener('change', handleFileSelect, false);

	function handleFileSelect(event) {
		const file = event.target.files[0];
		if (file) {
			const reader = new FileReader();
			reader.onload = function (e) {
				try {
					const data = JSON.parse(e.target.result);
					renderGraph(data);
				} catch (error) {
					alert('Error parsing JSON: ' + error.message);
				}
			};
			reader.readAsText(file);
		}
	}
});

function buildHierarchy(stepId, dependencies, stepMap) {
	const node = {
		id: stepId,
		...stepMap[stepId],
		children: []
	};

	const deps = dependencies[stepId] || [];
	node.children = deps.map(depId => buildHierarchy(depId, dependencies, stepMap));

	return node;
}



function renderGraph(data) {
	// Clear any existing SVG (if re-rendering)
	d3.select('#graph').selectAll('*').remove();
	d3.select('body').selectAll('.tooltip').remove(); // Remove any existing tooltips

	// Create a map from step number to step data
	const stepMap = {};
	data.steps.forEach(step => {
		stepMap[step.step] = step;
	});

	const width = 1400; // Increase width to accommodate the tree

	// Build the hierarchical data starting from Step 1
	const rootData = buildHierarchy('1', data.dependencies, stepMap);

	// Convert the hierarchical data into a D3 hierarchy object
	const root = d3.hierarchy(rootData);

	// Prepare an object to store circle radii for each node
	const circleRadii = {};

	// First pass to calculate circle sizes
	root.each(node => {
		// Determine if the node is a terminal step (no children)
		const isTerminal = !node.children || node.children.length === 0;

		// Determine if the step has reagents
		const hasReagents = node.data.reagents && node.data.reagents.length > 0;

		// Prepare array of molecules to display
		let moleculesToDisplay = [];

		if (hasReagents) {
			// Add reagents to the array
			node.data.reagents.forEach(reagent => {
				moleculesToDisplay.push({ smiles: reagent.smiles });
			});
		}

		if (isTerminal) {
			// Add reactants to the array
			node.data.reactants.forEach(reactant => {
				moleculesToDisplay.push({ smiles: reactant.smiles });
			});
		}

		// Add the products
		node.data.products.forEach(product => {
			moleculesToDisplay.push({ smiles: product.smiles });
		});

		// Calculate total height of molecules
		const moleculeCount = moleculesToDisplay.length;
		const moleculeHeight = 100; // Height of each molecule image
		const moleculeSpacing = 10; // Space between molecules
		const totalHeight = moleculeCount * moleculeHeight + (moleculeCount - 1) * moleculeSpacing;

		// Calculate the circle radius to encompass all molecules
		const circleRadius = Math.max(60, totalHeight / 2 + 20);

		// Store the circle radius
		circleRadii[node.data.id] = circleRadius;

		// Store the molecules to display for use later
		node.data.moleculesToDisplay = moleculesToDisplay;
		node.data.circleRadius = circleRadius;
	});

	// Set the tree layout with nodeSize
	const treeLayout = d3.tree()
		.nodeSize([0, 250]); // Horizontal spacing is fixed; vertical spacing will be calculated

	// Compute the new tree layout
	const treeData = treeLayout(root);

	// Reverse the y-coordinates to orient the tree from right to left
	treeData.descendants().forEach(d => {
		d.x = d.x; // vertical position
		d.y = -d.y; // invert horizontal position
	});

	// Adjust vertical positions to prevent overlap based on circle sizes
	let currentY = 0;

	treeData.each(node => {
		const radius = circleRadii[node.data.id];
		node.x = currentY + radius;
		currentY += radius * 2 + 30; // Add extra spacing between nodes
	});

	const height = currentY;

	const svg = d3.select('#graph').append('svg')
		.attr('width', width + 200) // Add extra space for labels
		.attr('height', height + 100);

	const g = svg.append('g')
		.attr('transform', `translate(${width - 100}, 50)`); // Start from the right side

	// Add links between nodes
	const linkGenerator = d3.linkHorizontal()
		.x(d => d.y)
		.y(d => d.x);

	const link = g.selectAll('.link')
		.data(treeData.links())
		.enter().append('path')
		.attr('class', 'link')
		.attr('d', linkGenerator)
		.attr('fill', 'none')
		.attr('stroke', '#555')
		.attr('stroke-width', 2);

	// Add nodes
	const node = g.selectAll('.node')
		.data(treeData.descendants())
		.enter().append('g')
		.attr('class', 'node')
		.attr('transform', d => `translate(${d.y},${d.x})`);

	// Adjust the node rendering
	node.each(function (d) {
		const group = d3.select(this);

		// Retrieve stored data
		const moleculesToDisplay = d.data.moleculesToDisplay;
		const circleRadius = d.data.circleRadius;

		// Add a circle behind the molecules
		group.append('circle')
			.attr('r', circleRadius)
			.attr('fill', '#f0f8ff')
			.attr('stroke', '#69b3a2')
			.attr('stroke-width', 2);

		// Start position to center the molecules vertically
		const moleculeCount = moleculesToDisplay.length;
		const moleculeHeight = 100; // Height of each molecule image
		const moleculeSpacing = 10; // Space between molecules
		const totalHeight = moleculeCount * moleculeHeight + (moleculeCount - 1) * moleculeSpacing;
		let yOffset = -totalHeight / 2;

		// Remove any existing content in the group (if any)
		group.selectAll('.molecule-group').remove();

		moleculesToDisplay.forEach(molecule => {
			try {
				const mol = OCL.Molecule.fromSmiles(molecule.smiles);
				const svgOptions = {
					suppressChiralText: true
				};
				let moleculeSVG = mol.toSVG(100, 100, 'node-molecule', svgOptions);

				// Parse the SVG string
				const parser = new DOMParser();
				const svgDoc = parser.parseFromString(moleculeSVG, 'image/svg+xml');
				const svgElement = svgDoc.documentElement;

				// Extract inner SVG content
				const innerContent = svgElement.innerHTML;

				// Append the inner content to the group
				group.append('g')
					.attr('class', 'molecule-group')
					.html(innerContent)
					.attr('transform', `translate(-50, ${yOffset})`); // Center horizontally, adjust vertically

			} catch (error) {
				// Handle error in rendering molecule
				group.append('text')
					.attr('x', 0)
					.attr('y', yOffset + 50)
					.attr('text-anchor', 'middle')
					.text('Error rendering molecule');
			}

			yOffset += moleculeHeight + moleculeSpacing;
		});
	});

	// Tooltip for node details
	const tooltip = d3.select('body').append('div')
		.attr('class', 'tooltip')
		.style('opacity', 0);

	node.on('mouseover', (event, d) => {
		tooltip.transition()
			.duration(200)
			.style('opacity', .9);
		tooltip.html(getStepDetails(d.data))
			.style('left', (event.pageX + 15) + 'px')
			.style('top', (event.pageY - 28) + 'px');
	})
		.on('mouseout', () => {
			tooltip.transition()
				.duration(500)
				.style('opacity', 0);
		});
}

// ... [Rest of the code remains unchanged] ...

function getStepDetails(step) {
	let html = `<strong>Step ${step.id}</strong><br/>`;

	html += `<em>Scalability Index:</em> ${step.reactionmetrics[0].scalabilityindex}<br/>`;
	html += `<em>Confidence Estimate:</em> ${step.reactionmetrics[0].confidenceestimate}<br/>`;

	return html;
}

