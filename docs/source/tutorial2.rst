.. title:: Tutorial 2

Tutorial 2
==========

Using the DeepRetro GUI for Retrosynthesis Analysis and Pathway Visualization

This tutorial guides you through analyzing and visualizing the retrosynthetic pathway for a bioactive heterocyclic compound using the DeepRetro GUI.

Getting Started
---------------

Before you begin, ensure you have:

- The DeepRetro backend server running
- An API key for backend access
- A modern web browser (Chrome/Firefox/Safari)

To start, navigate to `viewer/index.html`, enter your API key when prompted, and choose your analysis mode.

Entering the Target Molecule
---------------------------
.. image:: _static/tut2.png
   :alt: Example Pathway Visualization
   :align: center
   :width: 800px

Once the GUI is open, enter the SMILES notation for your target molecule in the input field:

  Example: ``O=C1N(CC2=C(F)C=C(C3=CC=CC4=NN(C)C=C43)C=C2F)CC5=NC=CC=C51``

Click "Analyze" to begin the retrosynthesis process. The system will validate the SMILES and start generating possible synthetic pathways.

Running the Analysis
--------------------

After clicking "Analyze," the backend will process your molecule and, if successful, display the retrosynthetic pathway as an interactive graph.

Pathway Visualization and Chemical Context
-----------------------------------------

Once the analysis is complete, you will see the retrosynthetic pathway for your molecule displayed as an interactive graph.

- **Explore the pathway:**
  - Click on molecule nodes to view their structure and details.
  - Hover over reaction arrows to see reaction conditions and success metrics.
  - Use navigation tools to switch between alternative routes if available.

.. image:: _static/tut2_pathway.png
   :alt: Example Pathway Visualization
   :align: center
   :width: 800px

This example demonstrates the synthesis of 6-(2,6-Difluoro-4-(2-methyl-2H-indazol-4-yl)benzyl)-6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-5-one, a bioactive heterocyclic compound with potential therapeutic applications, including as a kinase inhibitor.

**Two-Step Synthesis Approach**

  - **Step 1 (Suzuki Coupling):** Attach the 2-methylindazole moiety to the difluorobenzene ring via a Suzuki coupling reaction.
  - **Step 2 (SN2 Reaction):** Perform a simple SN2 reaction on the secondary amine group of 6,7-dihydro-5H-pyrrolo[3,4-b]pyridin-5-one to obtain the final product.

**Applications**

  - Potential kinase inhibitor
  - Therapeutic agent for various diseases
  - Useful scaffold for medicinal chemistry research

Exploring and Exporting Results
-------------------------------
.. image:: _static/tut2_info.png
   :alt: Example Pathway Visualization
   :align: center
   :width: 800px
- Click on molecule nodes to view detailed information (structure, formula, SMILES, etc.).
- Hover over reaction arrows to inspect conditions and success metrics.
- Export the pathway as JSON for further analysis or record-keeping.

Further Features
----------------

For advanced editing, troubleshooting, or more detailed features, refer to :doc:`Tutorial 1 <tutorial>`. 