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
					const rootStep = data['1']; // Start with step 1 as the root
					renderGraph(rootStep);
				} catch (error) {
					alert('Error parsing JSON: ' + error.message);
				}
			};
			reader.readAsText(file);
		}
	}
});

function renderGraph(rootStep) {
	// Clear any existing SVG (if re-rendering)
	d3.select('#graph').selectAll('*').remove();
	d3.select('body').selectAll('.tooltip').remove();

	// Build the tree structure using the new nested data format
	const buildTree = (step) => {
		const node = {
			...step.step,
			children: Object.values(step.children).map(childStep => buildTree(childStep))
		};
		return node;
	};

	const root = buildTree(rootStep);

	const graphContainer = d3.select('#graph')
		.style('overflow', 'auto')
		.style('position', 'relative');

	const width = 1600;
	const height = 1000;
	const circlePadding = 200;

	const svg = graphContainer.append('svg')
		.attr('width', width)
		.attr('height', height)
		.style('display', 'block')
		.style('margin', 'auto');

	const g = svg.append('g')
		.attr('transform', `translate(${circlePadding},${circlePadding})`);

	const treeLayout = d3.tree()
		.size([height - (2 * circlePadding), width - (2 * circlePadding)]);

	const treeData = treeLayout(d3.hierarchy(root));

	// Add reaction arrows first (so they appear behind nodes)
	const link = g.selectAll('.link')
		.data(treeData.links())
		.enter().append('g')
		.attr('class', 'link');

	// Draw paths between matching molecules
	link.append('path')
		.attr('d', d => {
			const parentStep = d.source.data;
			const childStep = d.target.data;

			// Find matching molecules between parent products and child reactants
			const matchingPairs = parentStep.products.flatMap(product =>
				childStep.reactants.map(reactant => ({
					product,
					reactant,
					matches: reactant.smiles === product.smiles
				}))
			).filter(pair => pair.matches);

			if (matchingPairs.length > 0) {
				return d3.linkHorizontal()
					.x(d => d.y)
					.y(d => d.x)({
						source: [d.source.y, d.source.x + (parentStep.products.length - 1) * 35 / 2],
						target: [d.target.y, d.target.x + (childStep.reactants.length - 1) * 35 / 2]
					});
			}
		})
		.attr('fill', 'none')
		.attr('stroke', '#555')
		.attr('stroke-width', 2)
		.attr('marker-end', 'url(#arrow)');

	// Add arrowhead marker
	svg.append('defs').append('marker')
		.attr('id', 'arrow')
		.attr('viewBox', '0 -5 10 10')
		.attr('refX', 8)
		.attr('refY', 0)
		.attr('markerWidth', 6)
		.attr('markerHeight', 6)
		.attr('orient', 'auto')
		.append('path')
		.attr('d', 'M0,-5L10,0L0,5')
		.attr('fill', '#555');

	// Add arrows between parent and child nodes
	link.append('line')
		.attr('x1', d => d.source.y)
		.attr('y1', d => d.source.x + (d.source.data.reactants.length - 1) * 35 / 2)
		.attr('x2', d => d.target.y)
		.attr('y2', d => d.target.x + (d.target.data.reactants.length - 1) * 35 / 2)
		.attr('stroke', '#888')
		.attr('stroke-width', 1)
		.attr('marker-end', 'url(#arrow)');

	// Add step labels
	link.append('text')
		.attr('x', d => (d.source.y + d.target.y) / 2)
		.attr('y', d => (d.source.x + d.target.x) / 2 - 10)
		.attr('text-anchor', 'middle')
		.text(d => `Step ${d.target.data.step}`);

	// Tooltip for reaction metrics
	const tooltip = d3.select('body').append('div')
		.attr('class', 'tooltip')
		.style('opacity', 0);

	link.on('mouseover', (event, d) => {
		tooltip.transition()
			.duration(200)
			.style('opacity', .9);
		tooltip.html(`
		<strong>Step ${d.target.data.step} Metrics</strong><br/>
		<em>Scalability Index:</em> ${d.target.data.reactionmetrics[0].scalabilityindex}<br/>
		<em>Confidence Estimate:</em> ${d.target.data.reactionmetrics[0].confidenceestimate}
	  `)
			.style('left', (event.pageX + 15) + 'px')
			.style('top', (event.pageY - 28) + 'px');
	})
		.on('mouseout', () => {
			tooltip.transition()
				.duration(500)
				.style('opacity', 0);
		});

	// Add nodes
	const node = g.selectAll('.node')
		.data(treeData.descendants())
		.enter().append('g')
		.attr('class', 'node')
		.attr('transform', d => `translate(${d.y},${d.x})`);

	// Render molecules
	node.each(function (d) {
		const group = d3.select(this);
		try {
			const molecules = [];

			if (d.data.step === '1') {
				// For step 1, show products first then reactants horizontally
				molecules.push(...d.data.products.map(m => ({ ...m, type: 'product' })));
				molecules.push(...d.data.reactants.map(m => ({ ...m, type: 'reactant' })));

				const xOffset = 100;
				molecules.forEach((molecule, i) => {
					const x = i * xOffset;

					const molGroup = group.append('g')
						.attr('transform', `translate(${x}, 0)`);

					molGroup.append('circle')
						.attr('r', 35)
						.attr('fill', molecule.type === 'product' ? '#e8f5e9' : '#f0f8ff')
						.attr('stroke', molecule.type === 'product' ? '#66bb6a' : '#69b3a2')
						.attr('stroke-width', 2);

					const mol = OCL.Molecule.fromSmiles(molecule.smiles);
					let molSVG = mol.toSVG(60, 60, 'molecule', { suppressChiralText: true });

					const parser = new DOMParser();
					const svgDoc = parser.parseFromString(molSVG, 'image/svg+xml');

					molGroup.append('g')
						.html(svgDoc.documentElement.innerHTML)
						.attr('transform', 'translate(-30, -30)');
				});
			} else {
				// For other steps, show reactants only vertically
				molecules.push(...d.data.reactants.map(m => ({ ...m, type: 'reactant' })));

				const yOffset = 70;
				molecules.forEach((molecule, i) => {
					const y = i * yOffset;

					const molGroup = group.append('g')
						.attr('transform', `translate(0, ${y})`);

					molGroup.append('circle')
						.attr('r', 35)
						.attr('fill', molecule.type === 'product' ? '#e8f5e9' : '#f0f8ff')
						.attr('stroke', molecule.type === 'product' ? '#66bb6a' : '#69b3a2')
						.attr('stroke-width', 2);

					const mol = OCL.Molecule.fromSmiles(molecule.smiles);
					let molSVG = mol.toSVG(60, 60, 'molecule', { suppressChiralText: true });

					const parser = new DOMParser();
					const svgDoc = parser.parseFromString(molSVG, 'image/svg+xml');

					molGroup.append('g')
						.html(svgDoc.documentElement.innerHTML)
						.attr('transform', 'translate(-30, -30)');
				});
			}
		} catch (error) {
			group.append('text')
				.attr('x', 0)
				.attr('y', 0)
				.text('Error rendering molecule');
		}
	});
}
