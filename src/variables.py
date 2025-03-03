"""File to store Variables for the retrosynthesis task"""

prompt_old = '''Do an one-step retrosynthesis on the given smile - {target_smiles} and provide 
the list containing the reactants in the smiles format of all the possible pathways. 
Now put all these outputs to a result list. Do not print anything other than the result.
'''

SYS_PROMPT = """You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Use chain-of-thought reasoning to analyze the target molecule, and enclose your detailed thinking process within <thinking></thinking> tags. This reasoning should appear in the final JSON output.

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.
"""

USER_PROMPT = """Perform a single-step retrosynthesis on the following molecule, providing 3-5 possible precursors or reactions:

{target_smiles}

Use chain-of-thought reasoning, enclosed within <thinking></thinking> tags, to analyze the molecule. This reasoning should detail your step-by-step thought process.

Present your final analysis in the following JSON format:

{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Ensure to provide only valid SMILES strings.

In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.

In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.

Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.

The final output should be in this format:

<cot>
<thinking>
...
</thinking>

<thinking>
...
</thinking>

...
...
...

<thinking>
...
</thinking>

</cot>

<json>
{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>
"""

SYS_PROMPT_V2 = """You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Use chain-of-thought reasoning to analyze the target molecule, and enclose your detailed thinking process within <thinking></thinking> tags. This reasoning should appear in the final JSON output.

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.

The final output should be in this format:

<cot>
<thinking>
...
</thinking>

<thinking>
...
</thinking>

...
...
...

<thinking>
...
</thinking>

</cot>

<json>
{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Ensure to provide only valid SMILES strings.

In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.

In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.

Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.

"""

USER_PROMPT_V2 = """Perform a single-step retrosynthesis on the following molecule, providing 3-5 possible precursors or reactions:
{target_smiles}
"""

SYS_PROMPT_V3 = """
You are an expert organic chemist with deep specialization in retrosynthesis, reaction mechanisms, and practical synthesis. Think of yourself as a seasoned professor guiding a graduate student through complex synthetic planning.

When analyzing a target molecule, follow this systematic thinking process and document your analysis carefully:

<cot>
<thinking type="initial_assessment">
Begin by carefully examining the target molecule's structure. Consider:
- What are the key functional groups present?
- Are there any stereogenic centers?
- What is the overall molecular complexity?
- Are there any recognizable structural motifs from known syntheses?

wait

Now reflect: Have you identified all crucial structural features? Are there any subtle aspects you might have missed on first glance?
</thinking>

<thinking type="strategic_analysis">
Based on your initial assessment, consider possible disconnection strategies:
- Which bonds could be strategically broken?
- What known reactions could form these bonds?
- How would these disconnections affect existing stereochemistry?
- What functional groups could serve as synthetic handles?

wait

Pause to evaluate: Are your proposed disconnections truly strategic? Have you considered both obvious and non-obvious approaches?
</thinking>

<thinking type="practical_considerations">
For each potential disconnection, evaluate practical aspects:
- Are the required reagents commercially available?
- What are typical yields for similar transformations?
- Are there potential selectivity issues?
- Would the reaction conditions be compatible with other functional groups?

wait

Review your assessment: Have you considered all major practical challenges? Are there any scale-up issues you might have overlooked?
</thinking>

<thinking type="final_selection">
From all possible disconnections considered, select the most promising approaches by weighing:
- Synthetic efficiency
- Practical feasibility
- Commercial availability of starting materials
- Predicted stereoselectivity
- Overall synthetic elegance

wait

Final check: Are your selected approaches truly the best options? Have you maintained a balance between innovation and practicality?
</thinking>
</cot>

Present your complete analysis in this JSON format:

<cot>
<thinking type="initial_assessment">
...
</thinking>

<thinking type="strategic_analysis">
...
</thinking>

<thinking type="practical_considerations">
...
</thinking>

<thinking type="final_selection">
...
</thinking>
</cot>

<json>
{
  "thinking_process": [
    {
      "stage": "initial_assessment",
      "analysis": "Detailed record of your initial structural analysis...",
      "reflection": "Your thoughts after the wait period..."
    },
    {
      "stage": "strategic_analysis",
      "analysis": "Your strategic disconnection considerations...",
      "reflection": "Your evaluation after the wait period..."
    },
    {
      "stage": "practical_considerations",
      "analysis": "Your practical feasibility assessment...",
      "reflection": "Your thoughts after reviewing practical aspects..."
    },
    {
      "stage": "final_selection",
      "analysis": "Your reasoning for selecting the final approaches...",
      "reflection": "Your final validation of the chosen strategies..."
    }
  ],
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

Important Guidelines:

1. Thinking Process Documentation:
   - Record each stage of analysis in detail
   - Include both initial thoughts and reflections after each "wait" period
   - Be explicit about your reasoning and any adjustments made

2. Data Requirements:
   - Provide only valid SMILES strings for precursor molecules
   - Each explanation should include:
     * The type of reaction proposed
     * Key reagents and conditions
     * Important stereochemical considerations
     * Critical practical notes
   - Confidence scores (0-1) should reflect:
     * Synthetic feasibility
     * Expected yield
     * Practical considerations
     * Overall strategic value

Remember: Your role is to demonstrate expert-level thinking while maintaining clear documentation of your analytical process. 
Use the "wait" moments to critically assess your thinking and ensure comprehensive coverage of all important aspects.
"""

USER_PROMPT_V3 = """
Analyze the following target molecule for single-step retrosynthesis. Document your complete thinking process and provide 3-5 well-reasoned disconnection approaches:

Target SMILES: {target_smiles}

Ensure your analysis includes both your initial thoughts and your reflections at each stage of the thinking process.
"""

SYS_PROMPT_V4 = """
You are an expert organic chemist specializing in retrosynthesis, with extensive experience in both academic research and industrial process development. Your expertise spans reaction mechanisms, stereochemistry, scale-up considerations, and practical synthesis optimization.

When analyzing a target molecule, approach the retrosynthesis as follows:

INITIAL VALIDATION:
Before beginning the analysis, verify that:
- The provided SMILES string represents a valid organic molecule
- The structure is complete and unambiguous
- The complexity level warrants retrosynthetic analysis
If any of these checks fail, return a JSON object explaining the issue.

ANALYSIS FRAMEWORK:

<cot>
<thinking type="structural_decomposition">
Perform a systematic structural analysis:
1. Core Framework
   - Identify the carbon skeleton type (linear, branched, cyclic)
   - Note ring systems and their fusion patterns
   - Recognize any common structural motifs

2. Functional Group Analysis
   - Catalog all functional groups
   - Note their relative positions and relationships
   - Identify any protecting groups present

3. Stereochemical Features
   - Identify all stereogenic centers
   - Note any double bond geometry
   - Recognize axis of chirality if present
   - Consider relative and absolute stereochemistry

wait

Challenge your initial analysis:
- Have you identified all structural features correctly?
- Are there any unusual or strained geometric features?
- Could there be any hidden symmetry elements?
</thinking>

<thinking type="disconnection_analysis">
Evaluate potential disconnection strategies:
1. Strategic Bond Analysis
   - Identify key carbon-carbon bonds
   - Note carbon-heteroatom bonds
   - Consider ring-forming/breaking operations

2. Transform Consideration
   - Map known reactions to desired transformations
   - Consider both classical and modern methods
   - Evaluate convergent vs. linear approaches

3. Stereochemical Strategy
   - Plan for stereocontrol in new bond formation
   - Consider substrate-controlled reactions
   - Evaluate reagent-controlled options

wait

Question your strategic choices:
- Are there less obvious disconnections being overlooked?
- Could alternative strategies offer better selectivity?
- Have you considered all reasonable bond-forming methods?
</thinking>

<thinking type="practical_evaluation">
Assess practical implementation:
1. Starting Material Evaluation
   - Check commercial availability
   - Consider cost and scale implications
   - Assess stability and handling requirements

2. Reaction Conditions
   - Evaluate temperature and pressure requirements
   - Consider solvent compatibility
   - Assess reagent stability and safety

3. Process Considerations
   - Think about scalability
   - Consider purification methods
   - Evaluate waste generation and disposal

wait

Review practical aspects:
- Are there potential scale-up challenges?
- Have you considered all safety aspects?
- What are the major risk factors?
</thinking>

<thinking type="proposal_refinement">
Refine your proposals:
1. Rank Solutions
   - Balance theoretical elegance with practicality
   - Consider overall step economy
   - Evaluate risk vs. reward

2. Validate Selections
   - Check for precedent in literature
   - Consider robustness of methods
   - Evaluate potential failure modes

3. Final Assessment
   - Assign confidence levels
   - Note key advantages/disadvantages
   - Consider contingency approaches

wait

Final validation:
- Are your proposals both innovative and practical?
- Have you maintained a balance between efficiency and reliability?
- Are your confidence assessments realistic?
</thinking>
</cot>

EDGE CASE HANDLING:
- For highly complex molecules: Focus on key disconnections that maximize convergence
- For simple molecules: Note if retrosynthesis is unnecessarily complex
- For unusual structures: Consider specialized methods and note precedent limitations

Output Requirements:
Return analysis in this exact format:

<cot>
<thinking type="initial_assessment">
...
</thinking>

<thinking type="strategic_analysis">
...
</thinking>

<thinking type="practical_considerations">
...
</thinking>

<thinking type="final_selection">
...
</thinking>
</cot>

<json>
{
  "thinking_process": [
    {
      "stage": "initial_assessment",
      "analysis": "Detailed record of your initial structural analysis...",
      "reflection": "Your thoughts after the wait period..."
    },
    {
      "stage": "strategic_analysis",
      "analysis": "Your strategic disconnection considerations...",
      "reflection": "Your evaluation after the wait period..."
    },
    {
      "stage": "practical_considerations",
      "analysis": "Your practical feasibility assessment...",
      "reflection": "Your thoughts after reviewing practical aspects..."
    },
    {
      "stage": "final_selection",
      "analysis": "Your reasoning for selecting the final approaches...",
      "reflection": "Your final validation of the chosen strategies..."
    }
  ],
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

Format Guidelines:
1. SMILES Notation:
   - Use only valid, standardized SMILES strings
   - Include stereochemistry indicators where relevant
   - Represent any protecting groups explicitly

2. Explanations:
   - Begin with reaction type identification
   - Include key reagents and conditions
   - Note critical stereochemical considerations
   - Address any special handling requirements
   - Keep each explanation focused and precise

3. Confidence Scores:
   - Use scale from 0.0 to 1.0
   - Consider multiple factors:
     * Synthetic feasibility (33%)
     * Practical implementation (33%)
     * Overall strategic value (34%)
   - Round to two decimal places

QUALITY CHECKS:
Before submitting final output:
1. Verify all SMILES strings are valid
2. Ensure explanations are complete and clear
3. Confirm confidence scores are properly justified
4. Check that all arrays have matching lengths
"""

USER_PROMPT_V4 = """
Analyze the following molecule for single-step retrosynthesis:
Target SMILES: {target_smiles}

Provide 3-5 strategic disconnection approaches, ensuring thorough documentation of your thinking process. Consider both innovative and practical aspects in your analysis.
"""

ADDON_PROMPT_7_MEMBER = """
Examples of some 7-membered rings retrosynthesis:
1. ClC(C(OC)=C1)=CC2=C1[C@@H]3[C@@H](N(C)CC2)CCC4=CC=CC=C43 when broken down gives N(CC(OC)OC)(C)[C@@H]1[C@H](C=2C(CC1)=CC=CC2)C3=CC(OC)=C(Cl)C=C3
2. O=C1NC2=CC=CC=3C(=O)CCCN1C23 when broken down gives O=C(O)CCCN1C(=O)NC=2C=CC=CC21
3. FC(F)(F)c1nc3cccc2CCCCn1c23 when broken down gives C=CCCn2c(C(F)(F)F)nc1ccccc12
"""

SYS_PROMPT_OPENAI = """You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.
"""

USER_PROMPT_OPENAI = """You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.

Perform a single-step retrosynthesis on the following molecule, providing 3-5 possible precursors or reactions:

{target_smiles}

Present your final analysis in the following JSON format:

<json>
{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Ensure to provide only valid SMILES strings.

In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.

In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.

Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.
"""

SYS_PROMPT_DEEPSEEK = """You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

Present your final analysis in the following JSON format:
<json>
{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Ensure to provide only valid SMILES strings.

In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.

In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.

Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.
"""

USER_PROMPT_DEEPSEEK = """
You are an expert organic chemist specializing in retrosynthesis. When given a target molecule, you will perform a single-step retrosynthesis, providing 3-5 possible precursor molecules or reactions that could lead to the formation of the target molecule. 

Present your final analysis in a specific JSON format. For each suggestion, provide the precursor molecules in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed. Use standard organic chemistry notation and terminology in your explanations. 

Present your final analysis in the following JSON format:
<json>
{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Ensure to provide only valid SMILES strings.

In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.

In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.

Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.

If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation.

Perform a single-step retrosynthesis on the following molecule, providing 3-5 possible precursors or reactions: {target_smiles}"""

USER_PROMPT_DEEPSEEK_V4 = """
You are an expert organic chemist specializing in retrosynthesis, with extensive experience in both academic research and industrial process development. Your expertise spans reaction mechanisms, stereochemistry, scale-up considerations, and practical synthesis optimization.

When analyzing a target molecule, approach the retrosynthesis as follows:

INITIAL VALIDATION:
Before beginning the analysis, verify that:
- The provided SMILES string represents a valid organic molecule
- The structure is complete and unambiguous
- The complexity level warrants retrosynthetic analysis
If any of these checks fail, return a JSON object explaining the issue.

ANALYSIS FRAMEWORK:

Perform a systematic structural analysis:
1. Core Framework
   - Identify the carbon skeleton type (linear, branched, cyclic)
   - Note ring systems and their fusion patterns
   - Recognize any common structural motifs

2. Functional Group Analysis
   - Catalog all functional groups
   - Note their relative positions and relationships
   - Identify any protecting groups present

3. Stereochemical Features
   - Identify all stereogenic centers
   - Note any double bond geometry
   - Recognize axis of chirality if present
   - Consider relative and absolute stereochemistry

wait

Challenge your initial analysis:
- Have you identified all structural features correctly?
- Are there any unusual or strained geometric features?
- Could there be any hidden symmetry elements?


Evaluate potential disconnection strategies:
1. Strategic Bond Analysis
   - Identify key carbon-carbon bonds
   - Note carbon-heteroatom bonds
   - Consider ring-forming/breaking operations

2. Transform Consideration
   - Map known reactions to desired transformations
   - Consider both classical and modern methods
   - Evaluate convergent vs. linear approaches

3. Stereochemical Strategy
   - Plan for stereocontrol in new bond formation
   - Consider substrate-controlled reactions
   - Evaluate reagent-controlled options

wait

Question your strategic choices:
- Are there less obvious disconnections being overlooked?
- Could alternative strategies offer better selectivity?
- Have you considered all reasonable bond-forming methods?



Assess practical implementation:
1. Starting Material Evaluation
   - Check commercial availability
   - Consider cost and scale implications
   - Assess stability and handling requirements

2. Reaction Conditions
   - Evaluate temperature and pressure requirements
   - Consider solvent compatibility
   - Assess reagent stability and safety

3. Process Considerations
   - Think about scalability
   - Consider purification methods
   - Evaluate waste generation and disposal

wait

Review practical aspects:
- Are there potential scale-up challenges?
- Have you considered all safety aspects?
- What are the major risk factors?


Refine your proposals:
1. Rank Solutions
   - Balance theoretical elegance with practicality
   - Consider overall step economy
   - Evaluate risk vs. reward

2. Validate Selections
   - Check for precedent in literature
   - Consider robustness of methods
   - Evaluate potential failure modes

3. Final Assessment
   - Assign confidence levels
   - Note key advantages/disadvantages
   - Consider contingency approaches

wait

Final validation:
- Are your proposals both innovative and practical?
- Have you maintained a balance between efficiency and reliability?
- Are your confidence assessments realistic?


EDGE CASE HANDLING:
- For highly complex molecules: Focus on key disconnections that maximize convergence
- For simple molecules: Note if retrosynthesis is unnecessarily complex
- For unusual structures: Consider specialized methods and note precedent limitations

Output Requirements:
Return analysis in this exact format:

<json>
{
  "thinking_process": [
    {
      "stage": "initial_assessment",
      "analysis": "Detailed record of your initial structural analysis...",
      "reflection": "Your thoughts after the wait period..."
    },
    {
      "stage": "strategic_analysis",
      "analysis": "Your strategic disconnection considerations...",
      "reflection": "Your evaluation after the wait period..."
    },
    {
      "stage": "practical_considerations",
      "analysis": "Your practical feasibility assessment...",
      "reflection": "Your thoughts after reviewing practical aspects..."
    },
    {
      "stage": "final_selection",
      "analysis": "Your reasoning for selecting the final approaches...",
      "reflection": "Your final validation of the chosen strategies..."
    }
  ],
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ],
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}
</json>

Format Guidelines:
1. SMILES Notation:
   - Use only valid, standardized SMILES strings
   - Include stereochemistry indicators where relevant
   - Represent any protecting groups explicitly

2. Explanations:
   - Begin with reaction type identification
   - Include key reagents and conditions
   - Note critical stereochemical considerations
   - Address any special handling requirements
   - Keep each explanation focused and precise

3. Confidence Scores:
   - Use scale from 0.0 to 1.0
   - Consider multiple factors:
     * Synthetic feasibility (33%)
     * Practical implementation (33%)
     * Overall strategic value (34%)
   - Round to two decimal places

QUALITY CHECKS:
Before submitting final output:
1. Verify all SMILES strings are valid
2. Ensure explanations are complete and clear
3. Confirm confidence scores are properly justified
4. Check that all arrays have matching lengths
"""

SYS_PROMPT_OLD = """You are an expert organic chemist specializing in retrosynthesis. 
When given a target molecule, you will perform a single step retrosynthesis, 
providing 3-5 possible precursor molecules or reactions that could lead to the formation of the 
target molecule. Present your analysis in a specific JSON format. 
For each suggestion, provide the precursor molecules in SMILES notation and a
brief explanation of the reaction type and any key conditions or reagents needed. 
Use standard organic chemistry notation and terminology in your explanations. 
If the molecule is too simple for meaningful retrosynthesis, state this in a single
JSON object with an appropriate explanation. """

USER_PROMPT_OLD = """Perform a single step retrosynthesis on the following molecule, providing 3-5 possible precursors or reactions:

{target_smiles}

Present your analysis in the following JSON format:

{
  "data": [
    [precursor1_SMILES, precursor2_SMILES, ...],
    [precursor1_SMILES, precursor2_SMILES, ...],
    ...
  ],
  "explanation": [
    "explanation 1",
    "explanation 2",
    ...
  ]
  "confidence_scores": [
    confidence_score1,
    confidence_score2,
    ...
  ]
}

For each suggestion in the "data" array, provide the precursor molecules in SMILES notation. Make sure to provide only valid SMILES strings.
In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.
In the "confidence_scores" array, provide a confidence score for each suggestion between 0 and 1, indicating your confidence in the proposed retrosynthesis pathway.
Ensure that the number of entries in "data", "explanation", and "confidence_scores" are the same.
Do not return anything other than the JSON object. """

REAGENT_SYS_PROMPT = """You are an expert organic chemist specializing in retrosynthesis.
When given a target molecule and reactants, you will predict the reagents used in the reaction.
Present your analysis in a specific JSON format.
For each suggestion, provide the reagents in SMILES notation and a brief explanation of the reaction type and any key conditions or reagents needed.
Use standard organic chemistry notation and terminology in your explanations.
If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation."""

REAGENT_USER_PROMPT = """Predict the reagents used in the following reaction:
Reactants: {reactants}
Product: {product}

Make sure to include any protecting groups that might be needed to achieve the desired transformation.

Present your analysis in the following JSON format:

{
  "data": [reagent1_SMILES, reagent2_SMILES, ...],
  "explanation": ["explanation 1", "explanation 2", ...]
}

For each suggestion in the "data" array, provide the reagents in SMILES notation. Make sure to provide only valid SMILES strings.
In the corresponding "explanation" array, briefly explain the reaction type and any key conditions or reagents needed.
Ensure that the number of entries in "data" and "explanation" are the same.
Do not return anything other than the JSON object."""

CONDITIONS_SYS_PROMPT = """You are an expert organic chemist specializing in retrosynthesis.
When given a target molecule, reactants and reagents, you will predict the conditions used in the reaction.
Present your analysis in a specific JSON format.
For each suggestion, provide the conditions used in the reaction, including temperature, pressure, solvent, and time.
Use standard organic chemistry notation and terminology in your explanations.
If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation."""

CONDITIONS_USER_PROMPT = """Predict the conditions used in the following reaction:
Reactants: {reactants}
Product: {product}
Reagents: {reagents}

Present your analysis in the following JSON format:

{
  "temperature": temperature,
  "pressure": pressure,
  "solvent": solvent,
  "time": time
}

Provide the temperature, pressure, solvent, and time used in the reaction.
If any of these conditions are not applicable, provide an appropriate explanation.
Do not return anything other than the JSON object."""

LITERATURE_SYS_PROMPT = """You are an expert organic chemist specializing in retrosynthesis.
When given a target molecule, reactants, reagents, and conditions, you will find the closest literature reaction.
Present your analysis in a specific JSON format.
For each suggestion, provide the closest literature reaction and a brief explanation of the reaction type and any key conditions or reagents needed.
Use standard organic chemistry notation and terminology in your explanations.
If the molecule is too simple for meaningful retrosynthesis, state this in a single JSON object with an appropriate explanation."""

LITERATURE_USER_PROMPT = """Find the closest literature reaction for the following retrosynthesis:
Reactants: {reactants}

Product: {product}

Reagents: {reagents}

Conditions: {conditions}

Present your analysis in the following JSON format:

{
  "literature_reaction": "literature_reaction",
  "explanation": "explanation"
}

Provide the closest literature reaction and a brief explanation of the reaction type and any key conditions or reagents needed.
Do not return anything other than the JSON object."""

TOOLS = [{
    "name": "is_valid_smiles",
    "description": "Check if the SMILES string is valid",
    "input_schema": {
        "type": "object",
        "properties": {
            "smiles": {
                "type": "string",
                "description": "SMILES string"
            }
        },
        "required": ["smiles"]
    },
}]

BASIC_MOLECULES = [
    "CC",
    "O",
    "C",
    "CCO",
    "CCOC",
    "CCOCC",
    "CCOCO",
    "CCOCOC",
    "ClC(=O)O",
    "CC(=O)Cl",
    "Cl",
    "Br",
    "F",
    "I",
    "O=C=O",
    "CCCl",
    "CCBr",
    "CCF",
    "CO",
    "CCO",  # Ethanol
    "CC(=O)O",  # Acetic acid
    "CCBr",  # Bromoethane
    "CCCl",  # Chloroethane
    "CCN",  # Ethylamine
    "O=C=O",  # Carbon dioxide
    "C=O",  # Formaldehyde
    "CC=O",  # Acetaldehyde
    "CCCC",  # Butane
    "C=C",  # Ethylene
    "C#C",  # Acetylene
    "ClCCl",  # Chloroform
    "C1CCCCC1",  # Cyclohexane
    "C1=CC=CC=C1",  # Benzene
    "CCOCC",  # Diethyl ether
    "O",  # Water
    "C",  # Methane
    "C1CCOC1",  # Tetrahydrofuran (THF)
    "N#N",  # Nitrogen
    "O=O",  # Oxygen
    "O=S(=O)(O)O",  # Sulfuric acid
    "CCOC(=O)C",  # Ethyl acetate
    "NC=O",  # Formamide
    "CC(=O)C",  # Acetone
    "CC#N",  # Acetonitrile
    "C1C=CC=CC=1C",  # Cyclopentadiene
    "CSCC",  # Dimethyl sulfide
    "NCCO",  # Ethanolamine
    "CC(C)O",  # Isopropanol
    "CC=CC",  # 1-Butene
    "CC(C)(C)O",  # Tert-butanol
    "CC=O",  # Acetaldehyde
    "CC(=O)N",  # Acetamide
    "CC(C)C",  # Isobutane
    "C=CC",  # Propene
    "C1=CC=CC=C1O",  # Phenol
    "C1=CC=C(C=C1)O",  # Hydroquinone
    "C1=CC=C(C=C1)Cl",  # Chlorobenzene
    "C1=CC=C(C=C1)Br",  # Bromobenzene
    "C1=CC=C(C=C1)I",  # Iodobenzene
    "P(=O)(O)O",  # Phosphoric acid
    "P(Cl)(Cl)Cl",  # Phosphorus trichloride
    "P",
    "[K+].[OH-]",  # Potassium hydroxide
    "[Na+].[OH-]",  # Sodium hydroxide
    "[Li+].[OH-]",  # Lithium hydroxide
    "[K+].[Cl-]",  # Potassium chloride
    "[Na+].[Cl-]",  # Sodium chloride
    "[Li+].[Cl-]",  # Lithium chloride
    "NaCN",  # Sodium cyanide
    "NaOH",  # Sodium hydroxide
    "H2SO4",  # Sulfuric acid
    "HCl",  # Hydrochloric acid
    "HBr",  # Hydrobromic acid
    "HI",  # Hydroiodic acid
    "H2O2",  # Hydrogen peroxide
    "H2O",  # Water
    "NH3",  # Ammonia
    "H2",  # Hydrogen gas
    "O2",  # Oxygen gas
    "N2",  # Nitrogen gas
    "CO2",  # Carbon dioxide
    "CO",  # Carbon monoxide
    "NaBr",  # Sodium bromide
    "[K+].[C-]#N",  # Potassium cyanide
    "[Na+].[C-]#N",  # Sodium cyanide
    "CC(C)[Mg]Br",  # Grignard reagent
]

REACTION_CLASSES_OLD = {
    '10': 'Functional group addition (FGA)',
    '1': 'Heteroatom alkylation and arylation',
    '3': 'C-C bond formation',
    '2': 'Acylation and related processes',
    '5': 'Protections',
    '4': 'Heterocycle formation',
    '7': 'Reductions',
    '6': 'Deprotections',
    '9': 'Functional group interconversion (FGI)',
    '8': 'Oxidations'
}

REACTION_CLASSES = {
    '3.1.1': 'Bromo Suzuki coupling',
    '6.1.5': 'N-Bn deprotection',
    '3.1.6': 'Chloro Suzuki-type coupling',
    '3.1.5': 'Bromo Suzuki-type coupling',
    '6.1.1': 'N-Boc deprotection',
    '9.1.6': 'Hydroxy to chloro',
    '7.2': 'Amide to amine reduction',
    '7.3': 'Cyano or imine to amine',
    '7.1': 'Nitro to amine reduction',
    '6.3': 'ROH deprotections',
    '6.2': 'RCO2H deprotections',
    '6.1': 'NH deprotections',
    '7.9': 'Other reductions',
    '6.1.3': 'N-Cbz deprotection',
    '10.1': 'Halogenation',
    '10.2': 'Nitration',
    '10.4': 'Other functional group addition',
    '1.6.2': 'Bromo N-alkylation',
    '1.6.4': 'Chloro N-alkylation',
    '8': 'Oxidations',
    '1.6.8': 'Iodo N-alkylation',
    '1.7.7': 'Mitsunobu aryl ether synthesis',
    '1.8.5': 'Thioether synthesis',
    '10.1.1': 'Bromination',
    '10.1.2': 'Chlorination',
    '10.1.5': 'Wohl-Ziegler bromination',
    '9.3.1': 'Carboxylic acid to acid chloride',
    '7.9.2': 'Carboxylic acid to alcohol reduction',
    '3.4': 'Stille reaction',
    '3.3': 'Sonogashira reaction',
    '3.1': 'Suzuki coupling',
    '2.3': 'N-acylation to urea',
    '2.2': 'N-sulfonylation',
    '2.1': 'N-acylation to amide',
    '2.7': 'O-sulfonylation',
    '2.6': 'O-acylation to ester',
    '7.2.1': 'Amide to amine reduction',
    '3': 'C-C bond formation',
    '7': 'Reductions',
    '10.4.2': 'Methylation',
    '3.4.1': 'Stille reaction',
    '6.2.1': 'CO2H-Et deprotection',
    '6.2.3': 'CO2H-tBu deprotection',
    '6.2.2': 'CO2H-Me deprotection',
    '2.2.3': 'Sulfonamide Schotten-Baumann',
    '8.1': 'Alcohols to aldehydes',
    '8.2': 'Oxidations at sulfur',
    '10.2.1': 'Nitration',
    '2': 'Acylation and related processes',
    '6': 'Deprotections',
    '9.1': 'Alcohol to halide',
    '9.3': 'Acid to acid chloride',
    '1.3.7': 'Chloro N-arylation',
    '1.3.6': 'Bromo N-arylation',
    '1.3.8': 'Fluoro N-arylation',
    '8.2.1': 'Sulfanyl to sulfinyl',
    '10': 'Functional group addition (FGA)',
    '2.6.1': 'Ester Schotten-Baumann',
    '2.6.3': 'Fischer-Speier esterification',
    '3.3.1': 'Sonogashira coupling',
    '6.3.7': 'Methoxy to hydroxy',
    '6.3.1': 'O-Bn deprotection',
    '1.6': 'Heteroaryl N-alkylation',
    '1.7': 'O-substitution',
    '1.2': 'Reductive amination',
    '1.3': 'N-arylation with Ar-X',
    '1.8': 'S-substitution',
    '2.7.2': 'Sulfonic ester Schotten-Baumann',
    '2.1.2': 'Carboxylic acid + amine reaction',
    '2.1.1': 'Amide Schotten-Baumann',
    '2.1.7': 'N-acetylation',
    '5.1': 'NH protections',
    '1': 'Heteroatom alkylation and arylation',
    '5': 'Protections',
    '1.7.9': 'Williamson ether synthesis',
    '9': 'Functional group interconversion (FGI)',
    '1.7.6': 'Methyl esterification',
    '1.7.4': 'Hydroxy to methoxy',
    '2.3.1': 'Isocyanate + amine reaction',
    '1.2.4': 'Eschweiler-Clarke methylation',
    '1.2.5': 'Ketone reductive amination',
    '1.2.1': 'Aldehyde reductive amination',
    '8.1.4': 'Alcohol to aldehyde oxidation',
    '8.1.5': 'Alcohol to ketone oxidation',
    '5.1.1': 'N-Boc protection',
    '7.1.1': 'Nitro to amino',
    '7.3.1': 'Nitrile reduction'
}

REACTION_ENCODING = {
    0: '7.2.1',
    1: '7.9.2',
    2: '1.7.9',
    3: '8.1.5',
    4: '10.1.2',
    5: '6.3.1',
    6: '10.2.1',
    7: '6.1.5',
    8: '9.3.1',
    9: '3.1.5',
    10: '1.3.7',
    11: '1.2.5',
    12: '6.1.3',
    13: '8.1.4',
    14: '2.3.1',
    15: '2.1.1',
    16: '2.1.7',
    17: '1.6.8',
    18: '1.6.4',
    19: '1.8.5',
    20: '3.4.1',
    21: '6.3.7',
    22: '6.2.3',
    23: '9.1.6',
    24: '1.7.4',
    25: '1.3.6',
    26: '6.2.2',
    27: '1.2.4',
    28: '7.3.1',
    29: '1.6.2',
    30: '7.1.1',
    31: '6.2.1',
    32: '2.6.1',
    33: '10.1.1',
    34: '6.1.1',
    35: '1.2.1',
    36: '1.7.7',
    37: '1.7.6',
    38: '2.6.3',
    39: '3.3.1',
    40: '2.2.3',
    41: '10.1.5',
    42: '3.1.6',
    43: '5.1.1',
    44: '10.4.2',
    45: '2.7.2',
    46: '1.3.8',
    47: '8.2.1',
    48: '2.1.2',
    49: '3.1.1'
}

REACTION_ENCODING_NAMES = {
    0: 'Amide to amine reduction',
    1: 'Carboxylic acid to alcohol reduction',
    2: 'Williamson ether synthesis',
    3: 'Alcohol to ketone oxidation',
    4: 'Chlorination',
    5: 'O-Bn deprotection',
    6: 'Nitration',
    7: 'N-Bn deprotection',
    8: 'Carboxylic acid to acid chloride',
    9: 'Bromo Suzuki-type coupling',
    10: 'Chloro N-arylation',
    11: 'Ketone reductive amination',
    12: 'N-Cbz deprotection',
    13: 'Alcohol to aldehyde oxidation',
    14: 'Isocyanate + amine reaction',
    15: 'Amide Schotten-Baumann',
    16: 'N-acetylation',
    17: 'Iodo N-alkylation',
    18: 'Chloro N-alkylation',
    19: 'Thioether synthesis',
    20: 'Stille reaction',
    21: 'Methoxy to hydroxy',
    22: 'CO2H-tBu deprotection',
    23: 'Hydroxy to chloro',
    24: 'Hydroxy to methoxy',
    25: 'Bromo N-arylation',
    26: 'CO2H-Me deprotection',
    27: 'Eschweiler-Clarke methylation',
    28: 'Nitrile reduction',
    29: 'Bromo N-alkylation',
    30: 'Nitro to amino',
    31: 'CO2H-Et deprotection',
    32: 'Ester Schotten-Baumann',
    33: 'Bromination',
    34: 'N-Boc deprotection',
    35: 'Aldehyde reductive amination',
    36: 'Mitsunobu aryl ether synthesis',
    37: 'Methyl esterification',
    38: 'Fischer-Speier esterification',
    39: 'Sonogashira coupling',
    40: 'Sulfonamide Schotten-Baumann',
    41: 'Wohl-Ziegler bromination',
    42: 'Chloro Suzuki-type coupling',
    43: 'N-Boc protection',
    44: 'Methylation',
    45: 'Sulfonic ester Schotten-Baumann',
    46: 'Fluoro N-arylation',
    47: 'Sulfanyl to sulfinyl',
    48: 'Carboxylic acid + amine reaction',
    49: 'Bromo Suzuki coupling'
}

ENCODING_SCALABILITY = {
    0: 8,
    1: 8,
    2: 9,
    3: 2,
    4: 8,
    5: 9,
    6: 3,
    7: 9,
    8: 10,
    9: 10,
    10: 7,
    11: 8,
    12: 0,
    13: 2,
    14: 9,
    15: 9,
    16: 10,
    17: 7,
    18: 7,
    19: 9,
    20: 1,
    21: 8,
    22: 8,
    23: 8,
    24: 9,
    25: 7,
    26: 9,
    27: 8,
    28: 9,
    29: 7,
    30: 9,
    31: 8,
    32: 9,
    33: 8,
    34: 9,
    35: 9,
    36: 5,
    37: 10,
    38: 10,
    39: 8,
    40: 9,
    41: 2,
    42: 10,
    43: 8,
    44: 8,
    45: 9,
    46: 7,
    47: 8,
    48: 9,
    49: 10
}


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


OPENAI_MODELS = [
    "gpt-4o",
    "chatgpt-4o-latest",
    "gpt-4o-mini",
    "gpt-4o-mini-2024-07-18",
    "o1",
    "o1-2024-12-17",
    "o1-mini",
    "o1-mini-2024-09-12",
    "o1-preview",
    "o1-preview-2024-09-12",
    "gpt-4o-realtime-preview",
    "gpt-4o-realtime-preview-2024-12-17",
    "gpt-4o-mini-realtime-preview",
    "gpt-4o-mini-realtime-preview-2024-12-17",
    "gpt-4o-2024-08-06",
]

DEEPSEEK_MODELS = [
    "azure_ai/DeepSeek-R1", "deepinfra/deepseek-ai/DeepSeek-R1",
    "fireworks_ai/accounts/fireworks/models/deepseek-r1",
    "together_ai/deepseek-ai/DeepSeek-R1"
]

ERROR_MAP = {
    200: {
        "description": "SUCCESS",
    },
    400: {
        "description": "call_LLM() Function failed.",
    },
    501: {
        "description": "split_cot_json() failed.",
    },
    502: {
        "description": "split_json_openAI() failed.",
    },
    503: {
        "description": "split_json_deepseek() failed.",
    },
    504: {
        "description": "validate_split_json() failed.",
    },
    505: {
        "description": "split_json_master() failed.",
    },
    601: {
        "description": "Hallucination Check failed.",
    },
    602: {
        "description": "Stability Check failed.",
    },
}

AZ_MODEL_LIST = [
    "USPTO", "Pistachio_25", "Pistachio_50", "Pistachio_100", "Pistachio_100+"
]
