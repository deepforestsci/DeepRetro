// app.js

// Wait for the DOM to be fully loaded
document.addEventListener('DOMContentLoaded', () => {
	const fileInput = document.getElementById('fileInput');
  
	fileInput.addEventListener('change', handleFileSelect, false);
  
	function handleFileSelect(event) {
	  const file = event.target.files[0];
	  if (file) {
		const reader = new FileReader();
		reader.onload = function(e) {
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
  
	// Build the hierarchical data starting from Step 1
	const root = buildHierarchy('1', data.dependencies, stepMap);
  
	const width = 1000;
	const height = 800;
  
	const svg = d3.select('#graph').append('svg')
	  .attr('width', width)
	  .attr('height', height);
  
	const g = svg.append('g')
	  .attr('transform', 'translate(100,50)'); // Adjust margins as needed
  
	// Create a tree layout
	const treeLayout = d3.tree()
	  .size([height - 100, width - 200]); // Adjust size if needed
  
	// Assigns the x and y position for the nodes
	const treeData = treeLayout(d3.hierarchy(root));
  
	// Reverse the y-coordinates to orient the tree from right to left
	treeData.descendants().forEach(d => {
	  d.y = width - d.y;
	});
  
	// Add links between nodes
	const link = g.selectAll('.link')
	  .data(treeData.links())
	  .enter().append('path')
	  .attr('class', 'link')
	  .attr('d', d3.linkHorizontal()
		.x(d => d.y)
		.y(d => d.x))
	  .attr('fill', 'none')
	  .attr('stroke', '#555')
	  .attr('stroke-width', 2);
  
	// Add nodes
	const node = g.selectAll('.node')
	  .data(treeData.descendants())
	  .enter().append('g')
	  .attr('class', 'node')
	  .attr('transform', d => `translate(${d.y},${d.x})`);
  
	// Render the molecule images on the nodes
	node.each(function(d) {
	  const group = d3.select(this);
  
	  // Add a circle behind the molecule for better visuals
	  group.append('circle')
		.attr('r', 55)
		.attr('fill', '#f0f8ff')
		.attr('stroke', '#69b3a2')
		.attr('stroke-width', 2);
  
	  // Try to render the product molecule
	  try {
		const product = d.data.products[0]; // Assuming the first product is the main one
		console.log(product);
		const mol = OCL.Molecule.fromSmiles(product.smiles);
  
		// Remove chirality labels
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
		  .html(innerContent)
		  .attr('transform', 'translate(-50, -50)'); // Center the molecule
  
	  } catch (error) {
		// Handle error in rendering molecule
		group.append('text')
		  .attr('x', 0)
		  .attr('y', 0)
		  .text('Error rendering molecule');
	  }
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
  
  function getStepDetails(step) {
	let html = `<strong>Step ${step.id}</strong><br/>`;
  
	html += `<em>Scalability Index:</em> ${step.reactionmetrics[0].scalabilityindex}<br/>`;
	html += `<em>Confidence Estimate:</em> ${step.reactionmetrics[0].confidenceestimate}<br/>`;
  
	return html;
  }
  