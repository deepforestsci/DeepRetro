// app.js

document.addEventListener('DOMContentLoaded', () => {
    const fileInput = document.getElementById('fileInput');
    fileInput.addEventListener('change', handleFileSelect, false);
});

function handleFileSelect(event) {
    const file = event.target.files[0];
    if (file) {
        const reader = new FileReader();
        reader.onload = function(e) {
            try {
                const data = JSON.parse(e.target.result);
                const processedTree = processData(data);
                const rootStep = processedTree['1'];
                renderGraph(rootStep);
            } catch (error) {
                alert('Error parsing JSON: ' + error.message);
            }
        };
        reader.readAsText(file);
    }
}

function processData(data) {
    const dependencyLength = Object.keys(data.dependencies).length;
    let parent_id_list = new Array(dependencyLength + 2).fill(-1);
    parent_id_list[0] = null;
    parent_id_list[1] = 0;

    for (const step in data.dependencies) {
        for (const child of data.dependencies[step]) {
            parent_id_list[parseInt(child)] = parseInt(step);
        }
    }

    let data_steps = data.steps;
    for (const step in data.dependencies) {
        data_steps[parseInt(step) - 1]['child_id'] = data.dependencies[step];
        data_steps[parseInt(step) - 1]['parent_id'] = parent_id_list[parseInt(step)];
    }

    function buildTree(data_steps, parent_id) {
        let tree = {};
        for (const step of data_steps) {
            if (step.parent_id === parent_id) {
                tree[parseInt(step.step)] = {
                    step: step,
                    children: buildTree(data_steps, parseInt(step.step))
                };
            }
        }
        return tree;
    }

    return buildTree(data_steps, 0);
}

function formatFormula(formula) {
    return formula.replace(/(\d+)/g, '<tspan baseline-shift="sub">$1</tspan>');
}

function renderGraph(rootStep) {
    d3.select('#graph').selectAll('*').remove();
    d3.select('body').selectAll('.tooltip').remove();

    const buildTree = (step) => {
        const node = {
            ...step.step,
            children: Object.values(step.children).map(childStep => buildTree(childStep))
        };
        return node;
    };

    const root = buildTree(rootStep);

    const styles = document.createElement('style');
    styles.textContent = `
        .tooltip {
            position: absolute;
            padding: 12px;
            background: rgba(255, 255, 255, 0.98);
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            pointer-events: none;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
            font-size: 14px;
            border: 1px solid #f0f0f0;
            max-width: 300px;
            z-index: 1000;
        }
        .molecule-node {
            transition: all 0.2s ease;
        }
        .molecule-node:hover circle {
            stroke-width: 3px;
            filter: brightness(0.98);
        }
        .link path, .link line {
            transition: all 0.2s ease;
        }
        .link:hover path, .link:hover line {
            stroke: #2196F3;
            stroke-width: 2.5px;
        }
        .mol-info text {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
            font-size: 12px;
            fill: #666;
        }
    `;
    document.head.appendChild(styles);

    const width = 1200;
    const height = 800;
    const circlePadding = 100;

    const svg = d3.select('#graph')
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .style('display', 'block')
        .style('margin', 'auto')
        .style('background', '#ffffff');

    const g = svg.append('g')
        .attr('transform', `translate(${circlePadding},${circlePadding})`);

    const defs = svg.append('defs');
    
    defs.append('linearGradient')
        .attr('id', 'productGradient')
        .attr('x1', '0%').attr('y1', '0%')
        .attr('x2', '100%').attr('y2', '100%')
        .selectAll('stop')
        .data([
            {offset: '0%', color: '#e8f5e9'},
            {offset: '100%', color: '#c8e6c9'}
        ])
        .enter().append('stop')
        .attr('offset', d => d.offset)
        .attr('stop-color', d => d.color);

    defs.append('linearGradient')
        .attr('id', 'reactantGradient')
        .attr('x1', '0%').attr('y1', '0%')
        .attr('x2', '100%').attr('y2', '100%')
        .selectAll('stop')
        .data([
            {offset: '0%', color: '#e3f2fd'},
            {offset: '100%', color: '#bbdefb'}
        ])
        .enter().append('stop')
        .attr('offset', d => d.offset)
        .attr('stop-color', d => d.color);

    const treeLayout = d3.tree()
        .size([height - (2 * circlePadding), width - (2 * circlePadding)]);

    const treeData = treeLayout(d3.hierarchy(root));

    const link = g.selectAll('.link')
        .data(treeData.links())
        .enter().append('g')
        .attr('class', 'link');

    link.append('path')
        .attr('d', d => {
            const parentStep = d.source.data;
            const childStep = d.target.data;

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

    defs.append('marker')
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

    link.append('line')
        .attr('x1', d => d.source.y)
        .attr('y1', d => d.source.x + (d.source.data.reactants.length - 1) * 35 / 2)
        .attr('x2', d => d.target.y)
        .attr('y2', d => d.target.x + (d.target.data.reactants.length - 1) * 35 / 2)
        .attr('stroke', '#888')
        .attr('stroke-width', 1)
        .attr('marker-end', 'url(#arrow)');

    link.append('text')
        .attr('x', d => (d.source.y + d.target.y) / 2)
        .attr('y', d => (d.source.x + d.target.x) / 2 - 10)
        .attr('text-anchor', 'middle')
        .style('font-family', '-apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif')
        .style('font-weight', '500')
        .text(d => `Step ${d.target.data.step}`);

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

    const node = g.selectAll('.node')
        .data(treeData.descendants())
        .enter().append('g')
        .attr('class', 'node')
        .attr('transform', d => `translate(${d.y},${d.x})`);

    node.each(function(d) {
        const group = d3.select(this);
        try {
            const molecules = [];

            group.append('text')
                .attr('x', 0)
                .attr('y', -50)
                .attr('text-anchor', 'middle')
                .style('font-family', '-apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif')
                .style('font-weight', '500')
                .text(`Step ${d.data.step}`);

            if (d.data.step === '1') {
                molecules.push(...d.data.products.map(m => ({ ...m, type: 'product' })));
                molecules.push(...d.data.reactants.map(m => ({ ...m, type: 'reactant' })));

                const xOffset = 80;
                molecules.forEach((molecule, i) => {
                    const x = i * xOffset;

                    const molGroup = group.append('g')
                        .attr('class', 'molecule-node')
                        .attr('transform', `translate(${x}, 0)`);

                    molGroup.append('circle')
                        .attr('r', 35)
                        .attr('fill', `url(#${molecule.type}Gradient)`)
                        .attr('stroke', molecule.type === 'product' ? '#66bb6a' : '#42a5f5')
                        .attr('stroke-width', 2)
                        .style('filter', 'drop-shadow(0px 2px 3px rgba(0,0,0,0.1))');

                    const mol = OCL.Molecule.fromSmiles(molecule.smiles);
                    let molSVG = mol.toSVG(60, 60, 'molecule', { suppressChiralText: true });

                    const parser = new DOMParser();
                    const svgDoc = parser.parseFromString(molSVG, 'image/svg+xml');

                    molGroup.append('g')
                        .html(svgDoc.documentElement.innerHTML)
                        .attr('transform', 'translate(-30, -30)');

                    // Add molecular info group
                    const infoGroup = molGroup.append('g')
                        .attr('class', 'mol-info');

                    const metadata = molecule.type === 'product' ? 
                        molecule.product_metadata : molecule.reactant_metadata;

                    infoGroup.append('text')
                        .attr('x', 0)
                        .attr('y', 45)
                        .attr('text-anchor', 'middle')
                        .text(`${metadata.mass.toFixed(1)} g/mol`);

                    infoGroup.append('text')
                        .attr('x', 0)
                        .attr('y', 60)
                        .attr('text-anchor', 'middle')
                        .html(formatFormula(metadata.chemical_formula));
                });
            } else {
                molecules.push(...d.data.reactants.map(m => ({ ...m, type: 'reactant' })));

                const yOffset = 60;
                molecules.forEach((molecule, i) => {
                    const y = i * yOffset;

                    const molGroup = group.append('g')
                        .attr('class', 'molecule-node')
                        .attr('transform', `translate(0, ${y})`);

                    molGroup.append('circle')
                        .attr('r', 35)
                        .attr('fill', `url(#${molecule.type}Gradient)`)
                        .attr('stroke', molecule.type === 'product' ? '#66bb6a' : '#42a5f5')
                        .attr('stroke-width', 2)
                        .style('filter', 'drop-shadow(0px 2px 3px rgba(0,0,0,0.1))');

                    const mol = OCL.Molecule.fromSmiles(molecule.smiles);
                    let molSVG = mol.toSVG(60, 60, 'molecule', { suppressChiralText: true });

                    const parser = new DOMParser();
                    const svgDoc = parser.parseFromString(molSVG, 'image/svg+xml');

                    molGroup.append('g')
                        .html(svgDoc.documentElement.innerHTML)
                        .attr('transform', 'translate(-30, -30)');

                    // Add molecular info group
                    const infoGroup = molGroup.append('g')
                        .attr('class', 'mol-info');

                    // Add molecular weight
                    infoGroup.append('text')
                        .attr('x', 0)
                        .attr('y', 45)
                        .attr('text-anchor', 'middle')
                        .text(`${molecule.reactant_metadata.mass.toFixed(1)} g/mol`);

                    // Add chemical formula
                    infoGroup.append('text')
                        .attr('x', 0)
                        .attr('y', 60)
                        .attr('text-anchor', 'middle')
                        .html(formatFormula(molecule.reactant_metadata.chemical_formula));
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