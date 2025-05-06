// app.js

document.addEventListener("DOMContentLoaded", () => {
  const fileInput = document.getElementById("fileInput");
  fileInput.addEventListener("change", handleFileSelect, false);
});

function updatePathwayNumber(pathwayNum) {
  try {
    const pathwayDisplay = document.getElementById("pathway-number");
    const currentPathway = document.getElementById("current-pathway");

    if (!pathwayDisplay || !currentPathway) {
      console.error("Pathway display elements not found");
      return;
    }

    // Always show the pathway number
    pathwayDisplay.style.display = "block";
    currentPathway.textContent = pathwayNum;

    // Make sure it's visible by bringing it to front
    pathwayDisplay.style.position = "relative";
    pathwayDisplay.style.zIndex = "1000";

    console.log("Updated pathway number to:", pathwayNum);
  } catch (error) {
    console.error("Error updating pathway number:", error);
  }
}

function handleFileSelect(event) {
  const file = event.target.files[0];
  if (file) {
    const reader = new FileReader();
    reader.onload = function (e) {
      try {
        const data = JSON.parse(e.target.result);
        const processedTree = processData(data);
        const rootStep = processedTree["0"];

        // Always show pathway 1 for new file
        updatePathwayNumber(1);

        renderGraph(rootStep);
      } catch (error) {
        alert("Error parsing JSON: " + error.message);
      }
    };
    reader.readAsText(file);
  }
}

function processData(data) {
  // --- Add log here ---
  console.log("[processData] Received data:", JSON.stringify(data, null, 2));
  // --------------------

  // Check if the data already contains an explicit Step 0 (handle string or number)
  // const hasExplicitStep0 = data.steps && data.steps.some(s => String(s.step) === '0'); // REMOVED CHECK
  // --- Add log here ---
  // console.log(`[processData] Does it have explicit Step 0? ${hasExplicitStep0}`); // REMOVED LOG
  // --------------------

  let data_steps;
  let newDependencies;

  // Define buildTree helper function (can be used by both branches)
  function buildTree(steps, parent_id) {
    let tree = {};
    for (const step of steps) {
      // Ensure parent_id comparison handles null/undefined/0 correctly
      const stepParentId = step.parent_id === undefined ? null : step.parent_id;
      if (stepParentId === parent_id) {
        // Process this step and its children
        const stepId = parseInt(step.step);
        tree[stepId] = {
          step: step, // Use the whole step object
          children: buildTree(steps, stepId),
        };

        // Log tree building information for debugging
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

  // REMOVED if (hasExplicitStep0) { ... } else { ... } structure
  // Always run the reconstruction logic:

  // Construct Step 0 and process dependencies
  console.log("[processData] Constructing Step 0 and processing dependencies.");

  // --- Ensure we handle cases where data.steps might be empty or null ---
  if (!data || !data.steps || data.steps.length === 0) {
    console.warn(
      "[processData] Input data has no steps. Returning empty tree."
    );
    return {}; // Return an empty object or handle as appropriate
  }
  // --- End handling empty steps ---

  const step1Data = data.steps[0]; // Assumes at least one step exists after check

  // --- Handle case where step 1 might not have products ---
  const greenProduct =
    step1Data.products && step1Data.products.length > 0
      ? step1Data.products[0]
      : null;
  if (!greenProduct) {
    console.warn(
      "[processData] Step 1 has no products defined. Step 0 will be empty."
    );
    // Decide how to handle this - maybe create a placeholder Step 0?
    // For now, we'll proceed but Step 0's product will be null.
  }
  // --- End handling no products in step 1 ---

  const step0 = {
    step: "0",
    products: greenProduct ? [greenProduct] : [], // Use the extracted product or empty array
    reactants: [],
    // --- Safely access potentially missing nested properties ---
    conditions: step1Data.conditions || {},
    reactionmetrics: step1Data.reactionmetrics || {},
    // --- End safe access ---
  };

  const step1Modified = {
    ...step1Data,
    step: "1",
    products: step1Data.products ? step1Data.products.slice(1) : [], // Take all products except first one
    parent_id: 0, // Explicitly set parent_id for the modified step 1
  };
  // --- Add log here ---
  console.log(
    "[processData] Constructed step0:",
    JSON.stringify(step0, null, 2)
  );
  console.log(
    "[processData] Constructed step1Modified:",
    JSON.stringify(step1Modified, null, 2)
  );
  // --------------------

  // Create new steps array with step 0
  // Filter out the original step 1 before adding the modified one
  const remainingSteps = data.steps
    .slice(1)
    .filter((s) => String(s.step) !== "1");
  const initialSteps = [step0, step1Modified, ...remainingSteps];

  // Update dependencies to start from step 0
  // Make sure data.dependencies exists before trying to access it
  const originalDependencies = data.dependencies || {};
  newDependencies = {
    0: ["1"], // Step 0 always points to step 1
  };

  // Copy all existing dependencies to preserve the tree structure
  Object.keys(originalDependencies).forEach((key) => {
    // Keep all dependencies except the original step 1's if any
    if (key !== "1") {
      newDependencies[key] = originalDependencies[key];
    }
  });

  // If step 1 had dependencies to other steps, preserve them
  if (originalDependencies["1"] && Array.isArray(originalDependencies["1"])) {
    newDependencies["1"] = originalDependencies["1"];
  }

  // --- Log new dependencies ---
  console.log(
    "[processData] Constructed newDependencies:",
    JSON.stringify(newDependencies, null, 2)
  );
  // --------------------------

  // Calculate parent_id list based on reconstructed dependencies
  let parent_id_map = {}; // Use a map for easier lookup { childId: parentId }
  parent_id_map[0] = null; // Step 0 has no parent
  parent_id_map[1] = 0; // Step 1's parent is 0

  for (const parentStep in newDependencies) {
    // Ensure the value is an array before iterating
    const children = Array.isArray(newDependencies[parentStep])
      ? newDependencies[parentStep]
      : [];
    for (const childStep of children) {
      // Ensure keys are treated as numbers for lookups if necessary
      const childNum = parseInt(childStep);
      const parentNum = parseInt(parentStep);
      // Avoid overwriting Step 1's parent if defined elsewhere
      if (childNum !== 1) {
        parent_id_map[childNum] = parentNum;
      }
    }
  }
  // --- Log parent map ---
  console.log(
    "[processData] Constructed parent_id_map:",
    JSON.stringify(parent_id_map, null, 2)
  );
  // ---------------------

  // Update steps array with correct parent_id and child_id information
  data_steps = initialSteps.map((step) => {
    const stepNum = parseInt(step.step);
    const parentId = parent_id_map[stepNum];

    // Calculate child_id based on dependencies
    const childIds = [];
    if (
      newDependencies[String(stepNum)] &&
      Array.isArray(newDependencies[String(stepNum)])
    ) {
      childIds.push(...newDependencies[String(stepNum)]);
    }

    return {
      ...step,
      parent_id: parentId !== undefined ? parentId : null, // Assign from map
      child_id: childIds, // Use calculated child IDs
    };
  });
  // --- Log final data_steps ---
  console.log(
    "[processData] Final data_steps before buildTree:",
    JSON.stringify(data_steps, null, 2)
  );
  // ---------------------------

  return buildTree(data_steps, null); // Build tree from newly processed steps
}

function calculateMoleculeSize(metadata) {
  if (!metadata || !metadata.chemical_formula) {
    return {
      radius: 35,
      svgSize: 60,
    };
  }

  const formula = metadata.chemical_formula;

  // Count total atoms including repeats
  const atomMatches = formula.match(/[A-Z][a-z]?\d*/g) || [];
  let totalAtoms = 0;
  atomMatches.forEach((match) => {
    const count = match.match(/\d+/);
    totalAtoms += count ? parseInt(count[0]) : 1;
  });

  // Count unique elements
  const uniqueElements = new Set(formula.replace(/[0-9]/g, "").split("")).size;

  // New calculation that better accounts for molecular size
  const baseRadius = 45; // Increased base radius
  const complexityFactor = Math.log(totalAtoms) * 20; // Steeper scaling
  const radius = Math.max(baseRadius, baseRadius + complexityFactor);

  // SVG size proportional to radius with larger multiplier for complex molecules
  const svgSize = totalAtoms > 50 ? radius * 2 : radius * 1.8;

  return {
    radius,
    svgSize,
  };
}

function calculateStepSize(molecules) {
  // Find largest molecule in step
  let maxRadius = 0;
  molecules.forEach((molecule) => {
    const metadata =
      molecule.type === "step0"
        ? molecule.product_metadata
        : molecule.reactant_metadata;
    const { radius } = calculateMoleculeSize(metadata);
    maxRadius = Math.max(maxRadius, radius);
  });
  return maxRadius;
}

function formatFormula(formula) {
  return formula.replace(/(\d+)/g, '<tspan baseline-shift="sub">$1</tspan>');
}

function renderGraph(rootStep) {
  d3.select("#graph").selectAll("*").remove();
  d3.select("body").selectAll(".tooltip").remove();

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

  const buildTree = (step) => {
    // Add diagnostic logging
    console.log(
      `[renderGraph.buildTree] Processing step ${step.step.step} with ${
        Object.keys(step.children).length
      } children`
    );

    const node = {
      ...step.step,
      children: Object.values(step.children).map((childStep) =>
        buildTree(childStep)
      ),
    };
    return node;
  };

  const root = buildTree(rootStep);
  console.log(
    `[renderGraph] Built tree with root step ${root.step}, children count: ${
      root.children ? root.children.length : 0
    }`
  );

  const styles = document.createElement("style");
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

  // Create hierarchy and calculate spacing
  let hierarchyRoot = d3.hierarchy(root);

  function calculateSpacingParams(hierarchyRoot) {
    let maxStepSize = 0;
    hierarchyRoot.each((d) => {
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
      const stepSize = calculateStepSize(molecules);
      maxStepSize = Math.max(maxStepSize, stepSize);
    });

    return {
      nodeSpacing: maxStepSize * 3,
      moleculeSpacing: maxStepSize * 2,
    };
  }

  const spacing = calculateSpacingParams(hierarchyRoot);

  const svg = d3
    .select("#graph")
    .append("svg")
    .style("display", "block")
    .style("margin", "auto")
    .style("background", "#ffffff");

  const g = svg.append("g");

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
  const treeLayout = d3
    .tree()
    .nodeSize([spacing.moleculeSpacing * 2, spacing.nodeSpacing])
    .separation((a, b) => {
      const aReactants = a.data.reactants ? a.data.reactants.length : 0;
      const bReactants = b.data.reactants ? b.data.reactants.length : 0;
      const maxReactants = Math.max(aReactants, bReactants);
      return (a.parent === b.parent ? 2 : 2.5) * (1 + maxReactants * 0.2);
    });

  // Apply layout to hierarchy
  hierarchyRoot = treeLayout(hierarchyRoot);

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

  const link = g
    .selectAll(".link")
    .data(hierarchyRoot.links())
    .enter()
    .append("g")
    .attr("class", "link");

  link
    .append("path")
    .attr("d", (d) => {
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
    .attr("fill", "none")
    .attr("stroke", "#999")
    .attr("stroke-width", 1.5)
    .attr("marker-end", "url(#arrow)");

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

  node.each(function (d) {
    const group = d3.select(this);
    try {
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

// Export main functions for testing and coverage
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
