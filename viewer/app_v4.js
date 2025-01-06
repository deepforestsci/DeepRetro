// app.js

document.addEventListener('DOMContentLoaded', () => {
    const fileInput = document.getElementById('fileInput');
    fileInput.addEventListener('change', handleFileSelect, false);
});

function handleFileSelect(event) {
    const file = event.target.files[0];
    if (file) {
        const reader = new FileReader();
        reader.onload = function (e) {
            try {
                const data = JSON.parse(e.target.result);
                const processedTree = processData(data);
                const rootStep = processedTree['0'];
                renderGraph(rootStep);
            } catch (error) {
                alert('Error parsing JSON: ' + error.message);
            }
        };
        reader.readAsText(file);
    }
}

function processData(data) {
    // Create step 0 from the green product in step 1
    const step1Data = data.steps[0];
    const greenProduct = step1Data.products[0];  // Assuming first product is green
    const step0 = {
        step: '0',
        products: [greenProduct],
        reactants: [],
        conditions: step1Data.conditions,
        reactionmetrics: step1Data.reactionmetrics
    };
    
    // Modify step 1 to only include blue molecules
    const step1Modified = {
        ...step1Data,
        step: '1',
        products: step1Data.products.slice(1),  // Take all products except first one
        parent_id: 0
    };

    // Create new steps array with step 0
    const newSteps = [step0, step1Modified, ...data.steps.slice(1)];

    // Update dependencies to start from step 0
    const newDependencies = {
        '0': ['1'],
        ...Object.entries(data.dependencies).reduce((acc, [key, value]) => {
            acc[key] = value;
            return acc;
        }, {})
    };

    // Process parent-child relationships
    const dependencyLength = Object.keys(newDependencies).length;
    let parent_id_list = new Array(dependencyLength + 2).fill(-1);
    parent_id_list[0] = null;
    parent_id_list[1] = 0;  // Step 1 now has Step 0 as parent

    for (const step in newDependencies) {
        for (const child of newDependencies[step]) {
            parent_id_list[parseInt(child)] = parseInt(step);
        }
    }

    // Update steps with parent-child information
    let data_steps = newSteps;
    for (const step in newDependencies) {
        const stepIndex = parseInt(step);
        if (stepIndex < data_steps.length) {
            data_steps[stepIndex]['child_id'] = newDependencies[step];
            data_steps[stepIndex]['parent_id'] = parent_id_list[stepIndex];
        }
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

    return buildTree(data_steps, null);
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
            z-index: 9999;
            transition: opacity 0.2s ease-in-out;
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

    const width = 1800;
    const height = 1200;
    const circlePadding = 200;
    const nodeSpacing = 300;
    const moleculeSpacing = 120;

    const svg = d3.select('#graph')
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .style('display', 'block')
        .style('margin', 'auto')
        .style('background', '#ffffff');

    const g = svg.append('g')
        .attr('transform', `translate(${circlePadding},${height/2})`);

    const defs = svg.append('defs');

    defs.append('linearGradient')
        .attr('id', 'step0Gradient')
        .attr('x1', '0%').attr('y1', '0%')
        .attr('x2', '100%').attr('y2', '100%')
        .selectAll('stop')
        .data([
            { offset: '0%', color: '#e8f5e9' },
            { offset: '100%', color: '#c8e6c9' }
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
            { offset: '0%', color: '#e3f2fd' },
            { offset: '100%', color: '#bbdefb' }
        ])
        .enter().append('stop')
        .attr('offset', d => d.offset)
        .attr('stop-color', d => d.color);

    const treeLayout = d3.tree()
        .size([height/3, width - (3 * circlePadding)])
        .nodeSize([moleculeSpacing * 2, nodeSpacing])
        .separation((a, b) => {
            return (a.parent === b.parent ? 1.5 : 2);
        });

    const treeData = treeLayout(d3.hierarchy(root));

    const tooltip = d3.select('body').append('div')
        .attr('class', 'tooltip')
        .style('opacity', 0);

    const link = g.selectAll('.link')
        .data(treeData.links())
        .enter().append('g')
        .attr('class', 'link');

    link.append('path')
        .attr('d', d => {
            const sourceX = d.source.x;
            const sourceY = d.source.y;
            const targetX = d.target.x;
            const targetY = d.target.y;
            
            // Calculate control points for the curved path
            const midY = (sourceY + targetY) / 2;
            
            return `M ${sourceY} ${sourceX}
                    C ${midY} ${sourceX},
                      ${midY} ${targetX},
                      ${targetY} ${targetX}`;
        })
        .attr('fill', 'none')
        .attr('stroke', '#999')
        .attr('stroke-width', 1.5)
        .attr('marker-end', 'url(#arrow)');

    link.on('mouseover', function(event, d) {
        const path = d3.select(this).select('path');
        path.attr('stroke', '#2196F3')
            .attr('stroke-width', 2.5);

        // Ensure we have the metrics data
        if (d.target.data && d.target.data.reactionmetrics && d.target.data.reactionmetrics[0]) {
            tooltip.style('opacity', 1)
                .html(`
                    <strong>Step ${d.target.data.step} Metrics</strong><br/>
                    <em>Scalability Index:</em> ${d.target.data.reactionmetrics[0].scalabilityindex}<br/>
                    <em>Confidence Estimate:</em> ${d.target.data.reactionmetrics[0].confidenceestimate}<br/>
                    <em>Closest Literature:</em> ${d.target.data.reactionmetrics[0].closestliterature}<br/>
                    <em>Reaction Conditions:</em> <br/>
                    <ul style="margin: 5px 0; padding-left: 20px;"> 
                        <li><em>Temperature:</em> ${d.target.data.conditions.temperature} </li>
                        <li><em>Pressure:</em> ${d.target.data.conditions.pressure} </li>
                        <li><em>Solvent:</em> ${d.target.data.conditions.solvent}</li>
                        <li><em>Time:</em> ${d.target.data.conditions.time} </li>
                    </ul>
                `)
                .style('left', (event.pageX + 15) + 'px')
                .style('top', (event.pageY - 28) + 'px');
        }
    })
    .on('mouseout', function() {
        const path = d3.select(this).select('path');
        path.attr('stroke', '#999')
            .attr('stroke-width', 1.5);
            
        tooltip.style('opacity', 0);
    });

    defs.append('marker')
        .attr('id', 'arrow')
        .attr('viewBox', '0 -5 10 10')
        .attr('refX', 20)
        .attr('refY', 0)
        .attr('markerWidth', 6)
        .attr('markerHeight', 6)
        .attr('orient', 'auto')
        .append('path')
        .attr('d', 'M0,-5L10,0L0,5')
        .attr('fill', '#999');

    const node = g.selectAll('.node')
        .data(treeData.descendants())
        .enter().append('g')
        .attr('class', 'node')
        .attr('transform', d => {
            const yOffset = d.data.step === '0' ? 0 : 
                          (d.data.reactants ? (d.data.reactants.length - 1) * moleculeSpacing / 2 : 0);
            return `translate(${d.y},${d.x - yOffset})`;
        });

    node.each(function(d) {
        const group = d3.select(this);
        try {
            const molecules = [];

            group.append('text')
                .attr('x', 0)
                .attr('y', -70)
                .attr('text-anchor', 'middle')
                .style('font-family', '-apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif')
                .style('font-weight', '500')
                .text(`Step ${d.data.step}`);

            if (d.data.step === '0' && d.data.products) {
                molecules.push(...d.data.products.map(m => ({ ...m, type: 'step0' })));
            } else if (d.data.reactants) {
                molecules.push(...d.data.reactants.map(m => ({ ...m, type: 'reactant' })));
            }

            molecules.forEach((molecule, i) => {
                if (!molecule || !molecule.smiles) return;

                const molGroup = group.append('g')
                    .attr('class', 'molecule-node')
                    .attr('transform', `translate(0, ${i * moleculeSpacing})`);

                molGroup.append('circle')
                    .attr('r', 35)
                    .attr('fill', `url(#${molecule.type}Gradient)`)
                    .attr('stroke', molecule.type === 'step0' ? '#66bb6a' : '#42a5f5')
                    .attr('stroke-width', 2)
                    .style('filter', 'drop-shadow(0px 2px 3px rgba(0,0,0,0.1))');

                molGroup.on('mouseover', function(event) {
                    d3.select(this).select('circle')
                        .attr('stroke-width', 3)
                        .style('filter', 'brightness(0.98)');

                    const metadata = molecule.type === 'step0' ? 
                        molecule.product_metadata : molecule.reactant_metadata;

                    // Get step metadata from parent node data
                    const stepMetadata = d.data;

                    if (metadata) {
                        tooltip.style('opacity', 1)
                            .html(`
                                <div style="padding: 8px;">
                                    <strong>Molecule Information</strong><br/>
                                    <em>Formula:</em> ${metadata.chemical_formula}<br/>
                                    <em>Mass:</em> ${metadata.mass.toFixed(1)} g/mol<br/>
                                    ${metadata.smiles ? `<em>SMILES:</em> ${metadata.smiles}` : ''}<br/>
                                    ${metadata.inchi ? `<em>InChI:</em> ${metadata.inchi}<br/>` : ''}
                                    <br/>
                                    <strong>Step ${d.data.step} Metrics</strong><br/>
                                    <em>Scalability Index:</em> ${stepMetadata.reactionmetrics[0].scalabilityindex}<br/>
                                    <em>Confidence Estimate:</em> ${stepMetadata.reactionmetrics[0].confidenceestimate}<br/>
                                    <em>Closest Literature:</em> ${stepMetadata.reactionmetrics[0].closestliterature}<br/>
                                    <em>Reaction Conditions:</em><br/>
                                    <ul style="margin: 5px 0; padding-left: 20px;"> 
                                        <li><em>Temperature:</em> ${stepMetadata.conditions.temperature}</li>
                                        <li><em>Pressure:</em> ${stepMetadata.conditions.pressure}</li>
                                        <li><em>Solvent:</em> ${stepMetadata.conditions.solvent}</li>
                                        <li><em>Time:</em> ${stepMetadata.conditions.time}</li>
                                    </ul>
                                </div>
                            `)
                            .style('left', (event.pageX + 15) + 'px')
                            .style('top', (event.pageY - 28) + 'px');
                    }
                })
                .on('mouseout', function() {
                    d3.select(this).select('circle')
                        .attr('stroke-width', 2)
                        .style('filter', 'drop-shadow(0px 2px 3px rgba(0,0,0,0.1))');
                        
                    tooltip.style('opacity', 0);
                });

                const mol = OCL.Molecule.fromSmiles(molecule.smiles);
                let molSVG = mol.toSVG(60, 60, 'molecule', { suppressChiralText: true });

                const parser = new DOMParser();
                const svgDoc = parser.parseFromString(molSVG, 'image/svg+xml');

                molGroup.append('g')
                    .html(svgDoc.documentElement.innerHTML)
                    .attr('transform', 'translate(-30, -30)');

                const infoGroup = molGroup.append('g')
                    .attr('class', 'mol-info');

                const metadata = molecule.type === 'step0' ?
                    molecule.product_metadata : molecule.reactant_metadata;

                if (metadata) {
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
                }
            });
        } catch (error) {
            console.error('Error rendering molecule:', error);
            group.append('text')
                .attr('x', 0)
                .attr('y', 0)
                .text('Error rendering molecule');
        }
    });
}