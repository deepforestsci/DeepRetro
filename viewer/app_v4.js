/**
 * DeepRetro Viewer
 *
 * This application provides an interactive visualization of chemical reaction pathways.
 * It allows users to:
 * 1. Upload JSON files containing reaction pathway data
 * 2. View molecular structures and reaction steps in a tree layout
 * 3. Interact with molecules to see detailed information
 * 4. Navigate through different reaction pathways
 *
 * The visualization uses D3.js for rendering the reaction tree and OpenChemLib for
 * molecular structure visualization.
 */

// app.js

document.addEventListener("DOMContentLoaded", () => {
  const fileInput = document.getElementById("fileInput");
  fileInput.addEventListener("change", handleFileSelect, false);
});

/**
 * Updates the pathway number display in the UI
 * @param {number} pathwayNum - The pathway number to display
 */
function updatePathwayNumber(pathwayNum) {
  try {
    // Get DOM elements for pathway display
    const pathwayDisplay = document.getElementById("pathway-number");
    const currentPathway = document.getElementById("current-pathway");

    // Validate DOM elements exist
    if (!pathwayDisplay || !currentPathway) {
      console.error("Pathway display elements not found");
      return;
    }

    // Show pathway number and ensure it's visible
    pathwayDisplay.style.display = "block";
    currentPathway.textContent = pathwayNum;

    // Ensure display is on top of other elements
    pathwayDisplay.style.position = "relative";
    pathwayDisplay.style.zIndex = "1000";

    console.log("Updated pathway number to:", pathwayNum);
  } catch (error) {
    console.error("Error updating pathway number:", error);
  }
}

/**
 * Handles file selection and initiates visualization
 * @param {Event} event - The file input change event
 */
function handleFileSelect(event) {
  // Validate event and files
  if (!event || !event.target || !event.target.files) {
    console.error("[handleFileSelect] Invalid event or missing files list");
    return;
  }

  const file = event.target.files[0];
  if (file) {
    const reader = new FileReader();
    reader.onload = function (e) {
      try {
        // Parse uploaded JSON file
        const data = JSON.parse(e.target.result);

        // Process data into hierarchical structure
        const processedTree = processData(data);
        const rootStep = processedTree["0"];

        // Always start with pathway 1 for new files
        updatePathwayNumber(1);

        // Render the reaction pathway graph
        renderGraph(rootStep);
      } catch (error) {
        alert("Error parsing JSON: " + error.message);
      }
    };
    reader.readAsText(file);
  } else {
    console.error("[handleFileSelect] No file selected");
  }
}

/**
 * Processes raw reaction data into a hierarchical structure
 * Adds Step 0 as root node and establishes parent-child relationships
 * @param {Object} data - Raw reaction pathway data
 * @returns {Object} Processed hierarchical data structure
 */
function processData(data) {
  // Log incoming data for debugging
  console.log("[processData] Received data:", JSON.stringify(data, null, 2));

  // Validate input data
  if (!data || !data.steps || data.steps.length === 0) {
    console.warn(
      "[processData] Input data has no steps. Returning empty tree."
    );
    return {};
  }

  // Helper function to build tree structure from steps
  function buildTree(steps, parent_id) {
    let tree = {};
    for (const step of steps) {
      // Handle null/undefined parent_id cases
      const stepParentId = step.parent_id === undefined ? null : step.parent_id;
      if (stepParentId === parent_id) {
        // Create tree node for this step and recursively process children
        const stepId = parseInt(step.step);
        tree[stepId] = {
          step: step,
          children: buildTree(steps, stepId),
        };

        // Log tree construction progress
        if (parent_id === null) {
          console.log(`[buildTree] Adding root step ${stepId}`);
        } else {
          console.log(
            `[buildTree] Adding step ${stepId} as child of ${parent_id}`
          );
        }
      }
    }
    return tree;
  }

  // Extract data from first step to construct Step 0
  const step1Data = data.steps[0];

  // Validate Step 1 has products
  if (!step1Data.products || step1Data.products.length === 0) {
    console.warn("[processData] Step 1 has no products defined");
  }

  // Get first product from Step 1 to use as Step 0's product
  const greenProduct =
    step1Data.products && step1Data.products.length > 0
      ? step1Data.products[0]
      : null;

  // Construct Step 0 as root node
  const step0 = {
    step: "0",
    products: greenProduct ? [greenProduct] : [],
    reactants: [],
    conditions: step1Data.conditions || {},
    reactionmetrics: step1Data.reactionmetrics || {},
  };

  // Modify Step 1 to connect with Step 0
  const step1Modified = {
    ...step1Data,
    step: "1",
    products: step1Data.products ? step1Data.products.slice(1) : [],
    parent_id: 0,
  };

  // Create new steps array with Step 0 and modified Step 1
  const remainingSteps = data.steps
    .slice(1)
    .filter((s) => String(s.step) !== "1");
  const initialSteps = [step0, step1Modified, ...remainingSteps];

  // Update dependencies to include Step 0
  const originalDependencies = data.dependencies || {};
  const newDependencies = {
    0: ["1"], // Step 0 always points to Step 1
  };

  // Preserve existing dependencies
  Object.keys(originalDependencies).forEach((key) => {
    newDependencies[key] = originalDependencies[key] || [];
  });

  // Calculate parent-child relationships
  let parent_id_map = {};
  parent_id_map[0] = null;
  parent_id_map[1] = 0;

  // Build parent-child map from dependencies
  for (const parentStep in newDependencies) {
    const children = Array.isArray(newDependencies[parentStep])
      ? newDependencies[parentStep]
      : [];
    for (const childStep of children) {
      const childNum = parseInt(childStep);
      const parentNum = parseInt(parentStep);
      if (childNum !== 1) {
        parent_id_map[childNum] = parentNum;
      }
    }
  }

  // Update steps with parent-child information
  const data_steps = initialSteps.map((step) => {
    const stepNum = parseInt(step.step);
    const parentId = parent_id_map[stepNum];
    const childIds = newDependencies[String(stepNum)] || [];

    return {
      ...step,
      parent_id: parentId !== undefined ? parentId : null,
      child_id: childIds,
    };
  });

  // Build and return final tree structure
  return buildTree(data_steps, null);
}

/**
 * Calculates appropriate size for molecule visualization based on complexity
 * @param {Object} metadata - Molecule metadata containing chemical formula
 * @returns {Object} Object containing radius and SVG size
 */
function calculateMoleculeSize(metadata) {
  // Return default size if no metadata available
  if (!metadata || !metadata.chemical_formula) {
    return {
      radius: 35,
      svgSize: 60,
    };
  }

  const formula = metadata.chemical_formula;

  // Count total atoms including repeating ones
  const atomMatches = formula.match(/[A-Z][a-z]?\d*/g) || [];
  let totalAtoms = 0;
  atomMatches.forEach((match) => {
    const count = match.match(/\d+/);
    totalAtoms += count ? parseInt(count[0]) : 1;
  });

  // Count unique elements for complexity
  const uniqueElements = new Set(formula.replace(/[0-9]/g, "").split("")).size;

  // Calculate radius based on molecular complexity
  const baseRadius = 45; // Base size for all molecules
  const complexityFactor = Math.log(totalAtoms) * 20; // Additional size based on atom count
  const radius = Math.max(baseRadius, baseRadius + complexityFactor);

  // Calculate SVG size with extra space for complex molecules
  const svgSize = totalAtoms > 50 ? radius * 2 : radius * 1.8;

  return {
    radius,
    svgSize,
  };
}

/**
 * Calculates the overall size needed for a reaction step
 * @param {Array} molecules - Array of molecules in the step
 * @returns {number} Maximum radius needed for the step
 */
function calculateStepSize(molecules) {
  // Find the largest molecule in the step to determine minimum space needed
  let maxRadius = 0;
  molecules.forEach((molecule) => {
    // Get appropriate metadata based on molecule type
    const metadata =
      molecule.type === "step0"
        ? molecule.product_metadata
        : molecule.reactant_metadata;

    // Calculate size and update maximum if needed
    const { radius } = calculateMoleculeSize(metadata);
    maxRadius = Math.max(maxRadius, radius);
  });
  return maxRadius;
}

/**
 * Formats chemical formulas with subscripts for proper display
 * @param {string} formula - Chemical formula to format
 * @returns {string} HTML formatted formula with subscripts
 */
function formatFormula(formula) {
  // Convert numbers in chemical formulas to subscripts using tspan
  return formula.replace(/(\d+)/g, '<tspan baseline-shift="sub">$1</tspan>');
}

/**
 * Main function for rendering the reaction pathway graph
 * Sets up the SVG canvas, renders nodes and edges, and adds interactivity
 * @param {Object} rootStep - Root node of the reaction pathway tree
 */
function renderGraph(rootStep) {
  // Clear existing graph and tooltips
  d3.select("#graph").selectAll("*").remove();
  d3.select("body").selectAll(".tooltip").remove();

  // Add container styles
  const containerStyles = document.createElement("style");
  containerStyles.textContent = `
    #graph {
      width: 100%;
      height: 800px;
      position: relative;
      overflow: hidden;
      border: 1px solid #ddd;
    }
    
    #graph button {
      width: 30px;
      height: 30px;
      margin: 2px;
      color: black;
      font-weight: bold;
      border: 1px solid #ddd;
      border-radius: 4px;
      background: white;
      cursor: pointer;
    }
    #graph button:hover {
      background: #f5f5f5;
    }
  `;
  document.head.appendChild(containerStyles);

  // Helper function to build D3 hierarchy
  const buildTree = (step) => {
    console.log(`[renderGraph.buildTree] Processing step ${step.step.step}`);
    const node = {
      ...step.step,
      children: Object.values(step.children).map((childStep) =>
        buildTree(childStep)
      ),
    };
    return node;
  };

  // Create D3 hierarchy from data
  const root = buildTree(rootStep);
  let hierarchyRoot = d3.hierarchy(root);

  // Calculate spacing based on molecule sizes
  function calculateSpacingParams(hierarchyRoot) {
    let maxStepSize = 0;
    hierarchyRoot.each((d) => {
      const molecules = [];
      // Collect molecules from step
      if (d.data.step === "0" && d.data.products) {
        molecules.push(
          ...d.data.products.map((m) => ({ ...m, type: "step0" }))
        );
      } else if (d.data.reactants) {
        molecules.push(
          ...d.data.reactants.map((m) => ({ ...m, type: "reactant" }))
        );
      }
      const stepSize = calculateStepSize(molecules);
      maxStepSize = Math.max(maxStepSize, stepSize);
    });

    return {
      nodeSpacing: maxStepSize * 3,
      moleculeSpacing: maxStepSize * 2,
    };
  }

  // Set up SVG and zoom behavior
  const svg = d3
    .select("#graph")
    .append("svg")
    .style("display", "block")
    .style("margin", "auto")
    .style("background", "#ffffff");

  const g = svg.append("g");

  // Add gradient definitions for molecule styling
  const defs = svg.append("defs");

  defs
    .append("linearGradient")
    .attr("id", "step0Gradient")
    .attr("x1", "0%")
    .attr("y1", "0%")
    .attr("x2", "100%")
    .attr("y2", "100%")
    .selectAll("stop")
    .data([
      { offset: "0%", color: "#e8f5e9" },
      { offset: "100%", color: "#c8e6c9" },
    ])
    .enter()
    .append("stop")
    .attr("offset", (d) => d.offset)
    .attr("stop-color", (d) => d.color);

  defs
    .append("linearGradient")
    .attr("id", "reactantGradient")
    .attr("x1", "0%")
    .attr("y1", "0%")
    .attr("x2", "100%")
    .attr("y2", "100%")
    .selectAll("stop")
    .data([
      { offset: "0%", color: "#e3f2fd" },
      { offset: "100%", color: "#bbdefb" },
    ])
    .enter()
    .append("stop")
    .attr("offset", (d) => d.offset)
    .attr("stop-color", (d) => d.color);

  // Create tree layout with dynamic spacing
  const spacing = calculateSpacingParams(hierarchyRoot);
  const treeLayout = d3
    .tree()
    .nodeSize([spacing.moleculeSpacing * 2, spacing.nodeSpacing])
    .separation((a, b) => {
      const aReactants = a.data.reactants ? a.data.reactants.length : 0;
      const bReactants = b.data.reactants ? b.data.reactants.length : 0;
      const maxReactants = Math.max(aReactants, bReactants);
      return (a.parent === b.parent ? 2 : 2.5) * (1 + maxReactants * 0.2);
    });

  // Apply layout and calculate bounds
  hierarchyRoot = treeLayout(hierarchyRoot);

  // Calculate SVG bounds based on node positions
  function calculateBounds(hierarchyRoot) {
    let minX = Infinity,
      maxX = -Infinity;
    let minY = Infinity,
      maxY = -Infinity;

    hierarchyRoot.each((d) => {
      minX = Math.min(minX, d.x);
      maxX = Math.max(maxX, d.x);
      minY = Math.min(minY, d.y);
      maxY = Math.max(maxY, d.y);
    });

    const padding = 100;
    return {
      width: maxY - minY + padding * 2,
      height: 4 * (maxX - minX + padding * 2),
      minX: minX - padding,
      minY: minY - padding,
    };
  }

  // Set SVG viewBox based on content bounds
  const bounds = calculateBounds(hierarchyRoot);
  svg.attr(
    "viewBox",
    `${bounds.minY} ${bounds.minX} ${bounds.width} ${bounds.height}`
  );

  // Add zoom behavior
  const zoom = d3
    .zoom()
    .scaleExtent([0.1, 4])
    .on("zoom", (event) => {
      g.attr("transform", event.transform);
    });

  svg.call(zoom);

  // Add zoom controls
  const zoomControls = d3
    .select("#graph")
    .append("div")
    .style("position", "absolute")
    .style("top", "10px")
    .style("right", "10px")
    .style("background", "white")
    .style("border", "1px solid #ddd")
    .style("border-radius", "4px")
    .style("padding", "5px");

  zoomControls
    .append("button")
    .text("+")
    .on("click", () => {
      svg.transition().duration(300).call(zoom.scaleBy, 1.5);
    });

  zoomControls
    .append("button")
    .text("-")
    .on("click", () => {
      svg.transition().duration(300).call(zoom.scaleBy, 0.75);
    });

  zoomControls
    .append("button")
    .text("⟲")
    .on("click", () => {
      svg.transition().duration(300).call(zoom.transform, d3.zoomIdentity);
    });

  const tooltip = d3
    .select("body")
    .append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  // Render links between nodes
  const link = g
    .selectAll(".link")
    .data(hierarchyRoot.links())
    .enter()
    .append("g")
    .attr("class", "link");

  // Add curved paths for links
  link
    .append("path")
    .attr("d", (d) => {
      const sourceX = d.source.x;
      const sourceY = d.source.y;
      const targetX = d.target.x;
      const targetY = d.target.y;
      const midY = (sourceY + targetY) / 2;
      return `M ${sourceY} ${sourceX}
              C ${midY} ${sourceX},
                ${midY} ${targetX},
                ${targetY} ${targetX}`;
    })
    .attr("fill", "none")
    .attr("stroke", "#999")
    .attr("stroke-width", 1.5)
    .attr("marker-end", "url(#arrow)");

  // Add hover effects and tooltips to links
  link
    .on("mouseover", function (event, d) {
      const path = d3.select(this).select("path");
      path.attr("stroke", "#2196F3").attr("stroke-width", 2.5);

      // Ensure we have the metrics data
      if (
        d.target.data &&
        d.target.data.reactionmetrics &&
        d.target.data.reactionmetrics[0]
      ) {
        tooltip
          .style("opacity", 1)
          .html(
            `
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
                `
          )
          .style("left", event.pageX + 15 + "px")
          .style("top", event.pageY - 28 + "px");
      }
    })
    .on("mouseout", function () {
      const path = d3.select(this).select("path");
      path.attr("stroke", "#999").attr("stroke-width", 1.5);

      tooltip.style("opacity", 0);
    });

  defs
    .append("marker")
    .attr("id", "arrow")
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 20)
    .attr("refY", 0)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0,-5L10,0L0,5")
    .attr("fill", "#999");

  // Render nodes
  const node = g
    .selectAll(".node")
    .data(hierarchyRoot.descendants())
    .enter()
    .append("g")
    .attr("class", "node")
    .attr("transform", (d) => {
      const yOffset =
        d.data.step === "0"
          ? 0
          : d.data.reactants
          ? ((d.data.reactants.length - 1) * spacing.moleculeSpacing) / 2
          : 0;
      return `translate(${d.y},${d.x - yOffset})`;
    });

  // Render molecules within nodes
  node.each(function (d) {
    const group = d3.select(this);
    try {
      // Collect molecules for this node
      const molecules = [];
      if (d.data.step === "0" && d.data.products) {
        molecules.push(
          ...d.data.products.map((m) => ({ ...m, type: "step0" }))
        );
      } else if (d.data.reactants) {
        molecules.push(
          ...d.data.reactants.map((m) => ({ ...m, type: "reactant" }))
        );
      }

      // Calculate step size and add step number
      const stepSize = calculateStepSize(molecules);
      group
        .append("text")
        .attr("x", 0)
        .attr("y", -stepSize * 1.5)
        .attr("text-anchor", "middle")
        .style(
          "font-family",
          '-apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif'
        )
        .style("font-weight", "500")
        .text(`Step ${d.data.step}`);

      // Render each molecule
      molecules.forEach((molecule, i) => {
        // Wrap individual molecule rendering in try...catch
        let molGroup; // Define molGroup here to be accessible in catch block
        try {
          if (!molecule || !molecule.smiles) {
            console.warn(
              `Skipping molecule ${i} in Step ${d.data.step}: Missing molecule or SMILES string.`
            );
            return; // Skip this iteration
          }

          molGroup = group
            .append("g") // Assign here
            .attr("class", "molecule-node")
            .attr("transform", `translate(0, ${i * spacing.moleculeSpacing})`);

          const metadata =
            molecule.type === "step0"
              ? molecule.product_metadata
              : molecule.reactant_metadata;

          const { radius, svgSize } = calculateMoleculeSize(metadata);

          molGroup
            .append("circle")
            .attr("r", radius)
            .attr("fill", `url(#${molecule.type}Gradient)`)
            .attr("stroke", molecule.type === "step0" ? "#66bb6a" : "#42a5f5")
            .attr("stroke-width", 2)
            .style("filter", "drop-shadow(0px 2px 3px rgba(0,0,0,0.1))");

          molGroup
            .on("mouseover", function (event) {
              d3.select(this)
                .select("circle")
                .attr("stroke-width", 3)
                .style("filter", "brightness(0.98)");

              if (metadata) {
                tooltip
                  .style("opacity", 1)
                  .html(
                    `
                                    <div style="padding: 8px;">
                                        <strong>Molecule Information</strong><br/>
                                        <em>Formula:</em> ${
                                          metadata.chemical_formula
                                        }<br/>
                                        <em>Mass:</em> ${metadata.mass.toFixed(
                                          1
                                        )} g/mol<br/>
                                        ${
                                          metadata.smiles
                                            ? `<em>SMILES:</em> ${metadata.smiles}`
                                            : ""
                                        }<br/>
                                        ${
                                          metadata.inchi
                                            ? `<em>InChI:</em> ${metadata.inchi}<br/>`
                                            : ""
                                        }
                                        <br/>
                                        <strong>Step ${
                                          d.data.step
                                        } Metrics</strong><br/>
                                        <em>Scalability Index:</em> ${
                                          d.data.reactionmetrics[0]
                                            .scalabilityindex
                                        }<br/>
                                        <em>Confidence Estimate:</em> ${
                                          d.data.reactionmetrics[0]
                                            .confidenceestimate
                                        }<br/>
                                        <em>Closest Literature:</em> ${
                                          d.data.reactionmetrics[0]
                                            .closestliterature
                                        }<br/>
                                        <em>Reaction Conditions:</em><br/>
                                        <ul style="margin: 5px 0; padding-left: 20px;"> 
                                            <li><em>Temperature:</em> ${
                                              d.data.conditions.temperature
                                            }</li>
                                            <li><em>Pressure:</em> ${
                                              d.data.conditions.pressure
                                            }</li>
                                            <li><em>Solvent:</em> ${
                                              d.data.conditions.solvent
                                            }</li>
                                            <li><em>Time:</em> ${
                                              d.data.conditions.time
                                            }</li>
                                        </ul>
                                    </div>
                                `
                  )
                  .style("left", event.pageX + 15 + "px")
                  .style("top", event.pageY - 28 + "px");
              }
            })
            .on("mouseout", function () {
              d3.select(this)
                .select("circle")
                .attr("stroke-width", 2)
                .style("filter", "drop-shadow(0px 2px 3px rgba(0,0,0,0.1))");

              tooltip.style("opacity", 0);
            });

          // --- Add log here to check SMILES before parsing ---
          console.log(
            `[renderGraph] Step ${d.data.step}, Molecule ${i}: Attempting to parse SMILES: '${molecule.smiles}'`
          );
          // --- Log the full object for the problematic case ---
          if (String(d.data.step) === "1" && i === 0) {
            try {
              // Use structured clone for a deep copy independent of original object reference
              const moleculeCopy = structuredClone(molecule);
              console.log(
                `[renderGraph] Full molecule object for Step 1, Molecule 0:`,
                moleculeCopy
              );
              // Also log the direct reference just in case
              console.log(
                `[renderGraph] Direct molecule object reference for Step 1, Molecule 0:`,
                molecule
              );
            } catch (cloneError) {
              console.warn(
                "[renderGraph] Could not structuredClone molecule object, logging directly:",
                molecule
              );
            }
          }
          // ---------------------------------------------------
          // ----------------------------------------------------

          // --- Force creation of a new string primitive to isolate from object context ---
          const smilesToParse = String(molecule.smiles); // Explicitly create a new string
          console.log(
            `[renderGraph] Step ${d.data.step}, Molecule ${i}: Explicit smilesToParse variable: '${smilesToParse}'`
          );
          // ---------------------------------------------------------------------------

          const mol = OCL.Molecule.fromSmiles(smilesToParse); // Use the isolated string
          // Check if molecule parsing was successful before generating SVG
          if (!mol || mol.getAllAtoms() === 0) {
            throw new Error(`OCL could not parse SMILES: ${smilesToParse}`); // Update error message too
          }
          let molSVG = mol.toSVG(svgSize, svgSize, "molecule", {
            suppressChiralText: true,
          });

          const parser = new DOMParser();
          const svgDoc = parser.parseFromString(molSVG, "image/svg+xml");

          // Check if SVG parsing was successful
          if (
            !svgDoc ||
            svgDoc.getElementsByTagName("parsererror").length > 0
          ) {
            const parserError = svgDoc.getElementsByTagName("parsererror")[0];
            console.error(
              "SVG Parser Error:",
              parserError
                ? parserError.textContent
                : "Unknown SVG parsing error"
            );
            throw new Error(
              `Could not parse generated SVG for SMILES: ${smilesToParse}`
            );
          }
          // --- End OCL Rendering ---

          molGroup
            .append("g")
            .html(svgDoc.documentElement.innerHTML)
            .attr("transform", `translate(-${svgSize / 2}, -${svgSize / 2})`);

          const infoGroup = molGroup.append("g").attr("class", "mol-info");

          if (metadata) {
            infoGroup
              .append("text")
              .attr("x", 0)
              .attr("y", radius + 10)
              .attr("text-anchor", "middle")
              .text(`${metadata.mass.toFixed(1)} g/mol`);

            infoGroup
              .append("text")
              .attr("x", 0)
              .attr("y", radius + 25)
              .attr("text-anchor", "middle")
              .html(formatFormula(metadata.chemical_formula));
          }
        } catch (molError) {
          console.error(
            `Error rendering molecule ${i} in Step ${d.data.step} (SMILES: ${
              molecule ? molecule.smiles : "N/A"
            }):`,
            molError
          );

          // If molGroup was created before the error, add error text to it
          if (molGroup) {
            molGroup
              .append("text")
              .attr("x", 0)
              .attr("y", 0) // Position roughly in the center of where the molecule would be
              .attr("text-anchor", "middle")
              .attr("fill", "red")
              .style("font-size", "10px")
              .style("font-weight", "bold")
              .text("[Render Error]");
          } else {
            // If molGroup wasn't even created, add error text to the main step group
            group
              .append("text")
              .attr("x", 0)
              .attr("y", i * spacing.moleculeSpacing) // Position at the molecule's y-offset
              .attr("text-anchor", "middle")
              .attr("fill", "red")
              .style("font-size", "10px")
              .style("font-weight", "bold")
              .text(`[Error Molecule ${i}]`);
          }
        }
      });
    } catch (error) {
      // This outer catch handles errors before the molecule loop (e.g., calculating step size)
      console.error(`Error processing step ${d.data.step}:`, error);
      group
        .append("text")
        .attr("x", 0)
        .attr("y", 0)
        .text("Error processing step");
    }
  });
}

// Export functions for testing and coverage
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    updatePathwayNumber,
    handleFileSelect,
    processData,
    calculateMoleculeSize,
    calculateStepSize,
    formatFormula,
    renderGraph,
  };
}
