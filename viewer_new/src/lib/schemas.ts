import { z } from "zod";

export const moleculeMetadataSchema = z
  .object({
    name: z.string().optional(),
    chemical_formula: z.string().optional(),
    mass: z.number().optional(),
    smiles: z.string().optional(),
    inchi: z.string().optional(),
  })
  .catchall(z.unknown());

export const moleculeRecordSchema = z
  .object({
    smiles: z.string(),
    reactant_metadata: moleculeMetadataSchema.optional(),
    reagent_metadata: moleculeMetadataSchema.optional(),
    product_metadata: moleculeMetadataSchema.optional(),
  })
  .catchall(z.unknown());

export const reactionMetricsSchema = z
  .object({
    scalabilityindex: z.union([z.string(), z.number()]).optional(),
    confidenceestimate: z.number().optional(),
    closestliterature: z.string().optional(),
  })
  .catchall(z.unknown());

export const pathwayStepSchema = z
  .object({
    step: z.union([z.string(), z.number()]).transform((value) => String(value)),
    reactants: z.array(moleculeRecordSchema).optional(),
    reagents: z.array(moleculeRecordSchema).optional(),
    products: z.array(moleculeRecordSchema).optional(),
    conditions: z.union([z.record(z.string(), z.unknown()), z.array(z.unknown())]).optional(),
    reactionmetrics: z.array(reactionMetricsSchema).optional(),
  })
  .catchall(z.unknown());

export const pathwayResultSchema = z
  .object({
    dependencies: z
      .record(z.string(), z.array(z.union([z.string(), z.number()])))
      .optional(),
    steps: z.array(pathwayStepSchema),
  })
  .catchall(z.unknown());

export const runtimeConfigSchema = z.object({
  instances: z.array(
    z.object({
      id: z.string(),
      label: z.string(),
      baseUrl: z.string(),
      defaults: z.object({
        model_type: z.string(),
        model_version: z.string(),
        advanced_prompt: z.boolean(),
        stability_flag: z.boolean(),
        hallucination_check: z.boolean(),
        use_protecting_group_feature: z.boolean().default(false),
      }),
    }),
  ),
  endpoints: z.object({
    retrosynthesis: z.string(),
    rerun: z.string(),
    partialRerun: z.string(),
    saveEdited: z.string(),
    health: z.string(),
  }),
});

export const advancedSettingsSchema = z.object({
  llm_models: z.record(
    z.string(),
    z.object({
      internal_name: z.string(),
      display_name: z.string(),
      supports_advanced_prompt: z.boolean(),
      supports_stability_check: z.boolean(),
      supports_hallucination_check: z.boolean(),
      supports_protecting_group_feature: z.boolean().optional(),
    }),
  ),
  az_models: z.record(z.string(), z.object({ display_name: z.string() })),
  defaults: z.object({
    model_type: z.string(),
    model_version: z.string(),
    advanced_prompt: z.boolean(),
    stability_flag: z.boolean(),
    hallucination_check: z.boolean(),
    use_protecting_group_feature: z.boolean().default(false),
  }),
});
