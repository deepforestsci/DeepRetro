# Ohuamine Hallucination Categories

Target: `CC(C)C[C@]12NC(=O)[C@@H]3[C@H](C)OC(=O)[C@](CC(C)C)(NC1=O)N32`

This note records the reviewed hallucination categories found in the
heuristic Ohuamine C pathways generated on 2026-06-11.

## Reviewed hallucinated steps

### Category A: Unsupported skeletal or ring-topology edit

Definition:
The proposed precursor-to-product jump changes the carbon skeleton,
side-chain topology, or ring connectivity without a chemically credible
source for the change.

Ohuamine example:
- Pathway `0`, step `1`
- Hallucination reason:
  one-carbon side-chain and topology change with no defensible leaving
  group or reagent source

Signals:
- carbon count or placement changes without a donor fragment
- ring closure/opening implied without a realistic bond-forming logic
- product connectivity cannot be explained by standard reduction,
  hydrolysis, or deprotection chemistry

### Category B: Incoherent protecting-group swap with fragment loss

Definition:
The step installs or swaps a protecting group while silently deleting an
existing protecting-group fragment, without including the deprotection
or fragment-removal chemistry required to make that possible.

Ohuamine example:
- Pathway `0`, step `4`
- Hallucination reason:
  Cbz-like benzyl carbonate fragment disappears while Boc appears, but
  only Boc-transfer reagents are present

Signals:
- benzyl, Boc, or other protecting-group fragments vanish without a
  hydrogenolysis, acid, or cleavage reagent
- net atom balance is only explainable if a missing deprotection step
  is inserted
- apparent one-step “protecting-group exchange” that actually requires
  two mechanistically distinct operations

### Category C: Missing external carbon or methylene source

Definition:
The step forms an exocyclic alkene or otherwise inserts a carbon unit
without any carbon- or methylene-donating reagent in the written
reactants.

Ohuamine example:
- Pathway `2`, step `3`
- Hallucination reason:
  ketone-to-exocyclic-alkene methylenation needs an external methylene
  source that is absent

Signals:
- ketone to terminal/exocyclic alkene conversion with no Wittig,
  Tebbe, Petasis, or equivalent reagent
- carbon count stays constant on paper while the bonding pattern implies
  new carbon content
- “olefination” proposed from only substrate + innocuous reagent

## Valid counter-patterns the model still overflags

These should be reinforced as valid examples during augmentation.

- Alkene hydrogenation on a preassembled Ohuamine-like cage
- Boc deprotection
- N-benzyl hydrogenolysis / debenzylation
- Intramolecular acyl-chloride aminolysis to close a lactam
- Acid-chloride formation from carboxylic acid
- Ester hydrolysis to free acid
- Reductive amination with benzaldehyde

## Current reviewed Ohuamine outcome

Reviewed heuristic sheet:
- `14` labelled steps
- `3` finalized hallucinations
- `11` finalized valid steps

Final hallucinated steps:
- Pathway `0`, step `1`
- Pathway `0`, step `4`
- Pathway `2`, step `3`
