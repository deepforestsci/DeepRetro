VALID_SMILE_STRING = 'CC(=O)CC'

VALID_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot>\n<thinking>\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\n</thinking>\n\n<thinking>\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\n</thinking>\n\n<thinking>\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.  \n</thinking>\n\n<thinking>\nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n</thinking>\n\n</cot>\n\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

EMPTY_RESPONSE = ""

COT_TAG_MISSING_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<thinking>\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\n</thinking>\n\n<thinking>\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\n</thinking>\n\n<thinking>\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.  \n</thinking>\n\n<thinking>\nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n</thinking>\n\n\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

COT_BODY_MISSING_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot></cot>\n\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

THINKING_TAG_MISSING_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot>\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.\nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n\n</cot>\n\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

THINKING_BODY_MISSING = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot>\n<thinking></thinking>\n\n<thinking></thinking>\n\n<thinking></thinking>\n\n<thinking></thinking>\n\n</cot>\n\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

COT_AND_THINKING_VALID_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.  \nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n<json>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n</json>'

JSON_TAG_MISSING_CLAUDE_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot>\n<thinking>\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\n</thinking>\n\n<thinking>\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\n</thinking>\n\n<thinking>\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.  \n</thinking>\n\n<thinking>\nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n</thinking>\n\n</cot>\n{\n  "data": [\n    ["CC(=O)C", "CC=O"],\n    ["CCC(=O)Cl", "CMgBr"],\n    ["CC(O)CC"]\n  ],\n  "explanation": [\n    "Aldol condensation of acetone and acetaldehyde, proceeding via enolate formation, nucleophilic addition, and dehydration",\n    "Grignard addition of methylmagnesium bromide to propanoyl chloride, followed by acidic workup",\n    "Oxidation of butan-2-ol, e.g. using PCC or Swern conditions"\n  ],\n  "confidence_scores": [\n    0.9,\n    0.7,\n    0.8\n  ]\n}\n'

EMPTY_JSON_BODY_RESPONSE = 'Here is the single-step retrosynthesis analysis for the molecule CC(=O)CC:\n\n<cot>\n<thinking>\nThe target molecule CC(=O)CC contains a ketone functional group. Possible retrosynthetic disconnections to consider are:\n1) Disconnection of the C-C bond adjacent to the ketone, which could arise from an aldol condensation reaction.\n2) Disconnection of the C-C bond on the other side of the ketone, which could come from a Grignard addition to a carboxylic acid derivative like an ester.\n3) Reduction of the ketone to an alcohol, which could then be derived from an oxidation of the corresponding secondary alcohol.\n</thinking>\n\n<thinking>\nFor the aldol disconnection, the precursors would be acetone (CC(=O)C) and acetaldehyde (CC=O). The reaction would proceed via enolate formation of the acetone, followed by nucleophilic addition to the acetaldehyde. A subsequent dehydration step would give the α,β-unsaturated ketone product.\n</thinking>\n\n<thinking>\nFor the Grignard addition, the precursors would be propanoyl chloride (CCC(=O)Cl) and methylmagnesium bromide (CMgBr). The Grignard reagent would add to the carbonyl, followed by an acidic workup to give the final ketone product.  \n</thinking>\n\n<thinking>\nFor the alcohol reduction, the precursor would be butan-2-ol (CC(O)CC). Oxidation, potentially using a chromium reagent like pyridinium chlorochromate (PCC) or a Swern oxidation, would convert the secondary alcohol to the ketone.\n</thinking>\n\n</cot>\n\n<json></json>'

# Advance prompt vars
CLAUDE_ADV_RESPONSE = '''<cot>
<thinking type="initial_assessment">
The target molecule CC(=O)CC has the following structural features:
- Linear 4-carbon chain 
- Terminal methyl ketone functional group
- No stereogenic centers or double bond geometry to consider
- Relatively simple, small molecule

The molecule is a valid target for retrosynthetic analysis. The key strategic consideration will be how to construct the carbon chain and install the ketone.
</thinking>

<thinking type="strategic_analysis">
Potential strategic disconnections to consider:
1. Ketone synthesis via oxidation
   - Disconnect the C-C bond adjacent to ketone 
   - Reveals a secondary alcohol precursor
   - Oxidation of the alcohol would give the ketone
   - Simple, robust approach with many reagent options

2. Carbon chain homologation 
   - Disconnect at the central C-C bond
   - Suggests an acetyl anion equivalent and electrophile
   - Aldol-type reaction or Grignard addition
   - Provides convergent, modular approach

3. Ketone synthesis via organometallic addition
   - Disconnect C-C bond adjacent to ketone
   - Implies an acetyl anion equivalent added to acetaldehyde
   - Grignard, organolithium or Weinreb ketone synthesis
   - Very direct approach with good precedent

The carbon chain homologation and organometallic addition strategies look most promising. The oxidation approach, while robust, requires an extra reduction step to access the alcohol.
</thinking>

<thinking type="practical_considerations">
Practical considerations for the two main strategies:

Carbon chain homologation:
- Acetyl anion equivalents like Weinreb amide are non-trivial 
- Aldol reaction requires enolate formation, low temp
- Purification of aldol products can be challenging
- Grignard addition is straightforward but moisture-sensitive

Organometallic addition:
- Grignard and organolithium reagents require anhydrous conditions
- Acetaldehyde is low boiling point, needs low temp
- Weinreb ketone synthesis very robust and scalable
- Organometallic reagents can be expensive on scale

Both strategies are viable with proper execution. The organometallic addition using a Weinreb amide is particularly attractive for its directness and reliability.
</thinking>

<thinking type="final_selection">
Based on the strategic analysis and practical considerations, I propose the following routes in order of preference:

1. Organometallic addition to acetaldehyde 
   - Use ethylmagnesium bromide or ethyllithium
   - Cryogenic temperatures for acetaldehyde
   - Aqueous workup and distillation
   - Highest confidence based on precedent and directness

2. Weinreb ketone synthesis
   - Form Weinreb amide from acetic acid 
   - Add ethylmagnesium bromide at 0 °C
   - Aqueous workup and extraction
   - High confidence based on robustness and scalability

3. Aldol reaction
   - Form enolate of acetone using LDA at -78 °C
   - Add acetaldehyde, then aqueous workup
   - Purify by distillation or chromatography
   - Moderate confidence based on purification challenges

The proposed routes balance strategic elegance with practical considerations. The organometallic addition is the most direct while the Weinreb ketone synthesis offers the best overall profile. The aldol route is a viable backup if the others fail.
</thinking>
</cot>

<json>
{
  "thinking_process": [
    {
      "stage": "initial_assessment",
      "analysis": "The target molecule CC(=O)CC has a linear 4-carbon chain with a terminal methyl ketone. There are no stereogenic centers or double bond geometry considerations. It is a relatively simple, small molecule that is a valid target for retrosynthetic analysis. The key strategic consideration will be constructing the carbon chain and installing the ketone.",
      "reflection": "The molecule is suitable for retrosynthetic analysis, with the key challenge being carbon-carbon bond formation and ketone installation."
    },
    {
      "stage": "strategic_analysis", 
      "analysis": "Three potential strategic disconnections are: 1) Ketone synthesis via oxidation of a secondary alcohol, disconnecting the C-C bond adjacent to the ketone 2) Carbon chain homologation using an acetyl anion equivalent and an electrophile, disconnecting the central C-C bond 3) Ketone synthesis via organometallic addition of an acetyl anion equivalent to acetaldehyde, disconnecting the C-C bond adjacent to the ketone.",
      "reflection": "The carbon chain homologation and organometallic addition strategies are most promising. The oxidation approach requires an extra reduction step."
    },
    {
      "stage": "practical_considerations",
      "analysis": "For carbon chain homologation, acetyl anion equivalents like Weinreb amides are non-trivial, aldol reactions require low temperatures and enolate formation, and purification can be challenging. Grignard addition is moisture sensitive. For organometallic addition, anhydrous conditions are required, acetaldehyde is low boiling and needs low temperature, but Weinreb ketone synthesis is robust and scalable. Organometallics can be expensive on scale.",
      "reflection": "Both main strategies are viable with proper execution. The organometallic addition using a Weinreb amide is particularly attractive for its directness and reliability."
    },
    {
      "stage": "final_selection",
      "analysis": "The proposed routes in order of preference are: 1) Organometallic addition to acetaldehyde using ethylmagnesium bromide or ethyllithium 2) Weinreb ketone synthesis by forming the Weinreb amide from acetic acid and adding ethylmagnesium bromide 3) Aldol reaction by forming the enolate of acetone with LDA and adding acetaldehyde. The organometallic addition is most direct, Weinreb ketone synthesis is most robust and scalable, and the aldol route is a viable backup.",
      "reflection": "The proposed routes balance strategic elegance with practical considerations. The organometallic addition is the most direct while the Weinreb ketone synthesis offers the best overall profile. Confidence is high in these routes based on precedent, robustness and scalability."
    }
  ],
  "data": [
    ["CCC(=O)C"],
    ["CCNC(=O)C", "CCC(=O)Cl"],
    ["CC(O)CC"]
  ],
  "explanation": [
    "Direct organometallic addition of ethyl Grignard or organolithium reagent to acetaldehyde, followed by aqueous workup and distillation",
    "Weinreb ketone synthesis by forming Weinreb amide from acetic acid, adding ethylmagnesium bromide, then aqueous workup and extraction", 
    "Aldol reaction by forming enolate of acetone with LDA, adding acetaldehyde, then aqueous workup and purification by distillation or chromatography"
  ],
  "confidence_scores": [
    0.90,
    0.85,
    0.70
  ]
}
</json>
'''