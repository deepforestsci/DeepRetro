# Valid sample smile string
VALID_SMILE_STRING = 'CC(=O)CC'

# Claude model
CLAUDE_MODEL = "claude-3-opus-20240229"

CLAUDE_ADV_MODEL = "claude-3-opus-20240229:adv"

# OpenAI model
OPENAI_MODEL = "gpt-4o"

OPENAI_ADV_MODEL = "gpt-4o:adv"

# Deepseek model
DEEPSEEK_MODEL = "deepinfra/deepseek-ai/DeepSeek-R1"

DEEPSEEK_ADV_MODEL = "deepinfra/deepseek-ai/DeepSeek-R1:adv"

DEEPSEEK_FIREWORKS_MODEL = "fireworks_ai/accounts/fireworks/models/deepseek-r1:adv"


# Valid claude model response
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



CLAUDE_ADV_RESPONSE_COT_TAG_MISSING = '''
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



CLAUDE_ADV_RESPONSE_THINKING_TAG_MISSING = '''<cot>
The target molecule CC(=O)CC has the following structural features:
- Linear 4-carbon chain 
- Terminal methyl ketone functional group
- No stereogenic centers or double bond geometry to consider
- Relatively simple, small molecule

The molecule is a valid target for retrosynthetic analysis. The key strategic consideration will be how to construct the carbon chain and install the ketone.

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



CLAUDE_ADV_RESPONSE_COT_BODY_MISSING = '''<cot>
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




OPENAI_ADV_RESPONSE = '''<cot>
<thinking type="initial_assessment">
1. Core Framework
   - The molecule is a linear carbon chain with four carbon atoms.
   - It contains a ketone functional group at the second carbon.

2. Functional Group Analysis
   - The primary functional group is a ketone (C=O).
   - There are no protecting groups or additional functional groups.

3. Stereochemical Features
   - There are no stereogenic centers or double bond geometry to consider.
   - The molecule is achiral.

Reflection:
The molecule is relatively simple, with a straightforward linear structure and a single ketone functional group. This simplicity suggests that the retrosynthetic analysis will focus on the disconnection of the ketone group or adjacent carbon-carbon bonds.
</thinking>

<thinking type="strategic_analysis">
1. Strategic Bond Analysis
   - Disconnection of the C-C bond adjacent to the carbonyl group (C2-C3) could lead to an acylation reaction.
   - Disconnection of the C-C bond between the carbonyl carbon and the methyl group (C1-C2) could suggest a Grignard reaction.

2. Transform Consideration
   - Acylation: The molecule could be formed by acylating a propyl group with acetyl chloride.
   - Grignard Reaction: The ketone could be formed by reacting a Grignard reagent (ethylmagnesium bromide) with acetaldehyde.
   - Aldol Reaction: An aldol condensation could form the ketone from acetaldehyde and acetone.

3. Stereochemical Strategy
   - No stereochemical control is required due to the lack of stereocenters.

Reflection:
The disconnection strategies focus on well-known reactions such as acylation and Grignard reactions, which are practical and efficient for forming ketones. The aldol reaction offers an alternative approach, though it may be less direct.
</thinking>

<thinking type="practical_considerations">
1. Starting Material Evaluation
   - Acetyl chloride, propyl group precursors, and Grignard reagents are commercially available and relatively inexpensive.
   - Acetaldehyde and acetone are also readily available.

2. Reaction Conditions
   - Acylation typically requires an acid chloride and a base, with mild conditions.
   - Grignard reactions require anhydrous conditions and careful handling of reagents.
   - Aldol reactions may require basic or acidic conditions, depending on the specific method.

3. Process Considerations
   - Acylation and Grignard reactions are scalable and well-documented.
   - Aldol reactions may require careful control to avoid side reactions.

Reflection:
The practical aspects of these reactions are favorable, with readily available starting materials and well-established conditions. The Grignard reaction requires more careful handling but is a robust method for ketone synthesis.
</thinking>

<thinking type="final_selection">
1. Rank Solutions
   - Acylation: High practicality and simplicity (Rank 1)
   - Grignard Reaction: High synthetic value, moderate handling complexity (Rank 2)
   - Aldol Reaction: Moderate practicality, less direct (Rank 3)

2. Validate Selections
   - Acylation and Grignard reactions have strong precedent in literature.
   - Aldol reactions are less common for this specific transformation but feasible.

3. Final Assessment
   - Acylation offers the best balance of simplicity and practicality.
   - Grignard reaction provides a robust alternative with good synthetic value.
   - Aldol reaction is a viable, though less direct, option.

Reflection:
The selected strategies are both innovative and practical, with a strong emphasis on well-established methods. The confidence in these approaches is high due to their simplicity and precedent.
</thinking>
</cot>

<json>
{
  "thinking_process": [
    {
      "stage": "initial_assessment",
      "analysis": "The molecule is a simple linear ketone with no stereochemistry or additional functional groups.",
      "reflection": "The simplicity suggests focusing on disconnection of the ketone or adjacent C-C bonds."
    },
    {
      "stage": "strategic_analysis",
      "analysis": "Considered acylation, Grignard, and aldol reactions for ketone formation.",
      "reflection": "Focused on practical and well-known reactions for ketone synthesis."
    },
    {
      "stage": "practical_considerations",
      "analysis": "Evaluated starting materials, reaction conditions, and scalability.",
      "reflection": "Reactions are practical with readily available materials and established conditions."
    },
    {
      "stage": "final_selection",
      "analysis": "Ranked acylation highest for simplicity, followed by Grignard and aldol reactions.",
      "reflection": "Selected strategies are practical and supported by precedent."
    }
  ],
  "data": [
    ["CCC", "CC(=O)Cl"],
    ["CC[Mg]Br", "CC=O"],
    ["CC(=O)C", "CC=O"]
  ],
  "explanation": [
    "Acylation: Propyl group acylated with acetyl chloride to form the ketone. Requires a base and mild conditions.",
    "Grignard Reaction: Ethylmagnesium bromide reacts with acetaldehyde to form the ketone. Requires anhydrous conditions.",
    "Aldol Reaction: Acetone and acetaldehyde undergo aldol condensation to form the ketone. Requires basic or acidic conditions."
  ],
  "confidence_scores": [
    0.95,
    0.90,
    0.80
  ]
}
</json>
'''

DEEPSEEK_ADV_RESPONSE = '''<cot>
<thinking type="initial_assessment">
The target molecule (3-pentanone, CC(=O)CC) features:
1. Linear 5-carbon chain with central ketone
2. No stereochemical complexity
3. Symmetrical C3 position
4. Potential for multiple straightforward disconnections

Reflection: Initial analysis confirms simple aliphatic ketone structure. Symmetry at'''











MOLS_SMALL = {
    "Favipiravir": "NC(=O)C1=NC(F)=CN=C1O",
    "Chlorzoxazone": "ClC1=CC2=C(OC(=O)N2)C=C1",
    "Metronidazole": "CC1=NC=C(N1CCO)[N+]([O-])=O",
    "Theophylline": "CN1C2=C(NC=N2)C(=O)N(C)C1=O",
    "Propylparaben": "CCCOC(=O)C1=CC=C(O)C=C1",
    "Secnidazole": "CC(O)CN1C(C)=NC=C1N(=O)=O",
    "Tocainide": "CC(N)C(=O)NC1=C(C)C=CC=C1C",
    "Enprofylline": "CCCN1C2=C(NC=N2)C(=O)NC1=O",
    "Milrinone": "CC1=C(C=C(C#N)C(=O)N1)C1=CC=NC=C1",
    "Methyldopa": "C[C@](N)(CC1=CC=C(O)C(O)=C1)C(O)=O",
    "Methoxamine": "COC1=CC(C(O)C(C)N)=C(OC)C=C1",
    "Zonisamide": "NS(=O)(=O)CC1=NOC2=CC=CC=C12",
    "Flavone": "O=C1C=C(OC2=CC=CC=C12)C1=CC=CC=C1",
    "Diazoxide": "CC1=NS(=O)(=O)C2=C(N1)C=CC(Cl)=C2",
    "Naproxen": "COC1=CC2=C(C=C1)C=C(C=C2)[C@H](C)C(O)=O",
    "Nalidixic acid": "CCN1C=C(C(O)=O)C(=O)C2=C1N=C(C)C=C2",
    "Phenobarbital": "CCC1(C(=O)NC(=O)NC1=O)C1=CC=CC=C1",
    "Melatonin": "COC1=CC2=C(NC=C2CCNC(C)=O)C=C1",
    "Sulbactam": "[H][C@@]12CC(=O)N1[C@@]([H])(C(O)=O)C(C)(C)S2(=O)=O",
    "Riluzole": "NC1=NC2=C(S1)C=C(OC(F)(F)F)C=C2",
    "Didanosine": "OC[C@@H]1CC[C@@H](O1)N1C=NC2=C1NC=NC2=O",
    "Tilbroquinol": "CC1=C2C=CC=NC2=C(O)C(Br)=C1",
    "Sulfadiazine": "NC1=CC=C(C=C1)S(=O)(=O)NC1=NC=CC=N1",
    "Phenprocoumon": "CCC(C1=CC=CC=C1)C1=C(O)C2=C(OC1=O)C=CC=C2",
}