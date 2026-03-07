export type ToggleableModelFlags = {
  advanced_prompt: boolean;
  stability_flag: boolean;
  hallucination_check: boolean;
  use_protecting_group_feature: boolean;
};

export type InstanceDefaults = ToggleableModelFlags & {
  model_type: string;
  model_version: string;
};

export type ViewerRuntimeConfig = {
  instances: Array<{
    id: string;
    label: string;
    baseUrl: string;
    defaults: InstanceDefaults;
  }>;
  endpoints: {
    retrosynthesis: string;
    rerun: string;
    partialRerun: string;
    saveEdited: string;
    health: string;
  };
};

export type AdvancedSettingsConfig = {
  llm_models: Record<
    string,
    {
      internal_name: string;
      display_name: string;
      supports_advanced_prompt: boolean;
      supports_stability_check: boolean;
      supports_hallucination_check: boolean;
      supports_protecting_group_feature?: boolean;
    }
  >;
  az_models: Record<string, { display_name: string }>;
  defaults: InstanceDefaults;
};

export type MoleculeMetadata = {
  name?: string;
  chemical_formula?: string;
  mass?: number;
  smiles?: string;
  inchi?: string;
  [key: string]: unknown;
};

export type MoleculeRecord = {
  smiles: string;
  reactant_metadata?: MoleculeMetadata;
  reagent_metadata?: MoleculeMetadata;
  product_metadata?: MoleculeMetadata;
  [key: string]: unknown;
};

export type ReactionMetrics = {
  scalabilityindex?: string | number;
  confidenceestimate?: number;
  closestliterature?: string;
  [key: string]: unknown;
};

export type PathwayStep = {
  step: string;
  reactants?: MoleculeRecord[];
  reagents?: MoleculeRecord[];
  products?: MoleculeRecord[];
  conditions?: Record<string, unknown> | unknown[];
  reactionmetrics?: ReactionMetrics[];
  [key: string]: unknown;
};

export type PathwayResult = {
  dependencies?: Record<string, Array<string | number>>;
  steps: PathwayStep[];
  [key: string]: unknown;
};

export type InstanceRequestSettings = InstanceDefaults;

export type AnalysisRequest = InstanceRequestSettings & {
  smiles: string;
};

export type ViewerMode = "analyze" | "upload";
export type ViewerRunSource = "api" | "file";
export type ViewerRunStatus = "idle" | "pending" | "success" | "error";
export type HealthState = "unknown" | "healthy" | "unauthorized" | "error";

export type ValidationIssue = {
  path: string;
  message: string;
};

export type NormalizedMolecule = {
  key: string;
  role: "product" | "reactant" | "reagent";
  smiles: string;
  metadata?: MoleculeMetadata;
  label: string;
  formula?: string;
  mass?: number;
};

export type NormalizedStepNode = {
  id: string;
  stepId: string;
  title: string;
  isVirtualRoot: boolean;
  rawStep?: PathwayStep;
  products: NormalizedMolecule[];
  reactants: NormalizedMolecule[];
  reagents: NormalizedMolecule[];
  metrics?: ReactionMetrics;
  conditions?: Record<string, unknown> | unknown[];
  childIds: string[];
  parentIds: string[];
};

export type NormalizedPathwayGraph = {
  stepMap: Record<string, PathwayStep>;
  dependencyMap: Record<string, string[]>;
  nodes: NormalizedStepNode[];
  edges: Array<{
    id: string;
    source: string;
    target: string;
  }>;
  virtualRoot: NormalizedStepNode;
};

export type ViewerRun = {
  key: string;
  label: string;
  source: ViewerRunSource;
  status: ViewerRunStatus;
  instanceId?: string;
  request?: AnalysisRequest;
  rawResult?: PathwayResult;
  graph?: NormalizedPathwayGraph;
  error?: string;
  durationMs?: number;
  dirty: boolean;
  lastUpdatedAt: number;
};

export type HealthStatusMap = Record<
  string,
  {
    state: HealthState;
    message: string;
  }
>;

export type ViewerStoreState = {
  runtimeConfig?: ViewerRuntimeConfig;
  advancedSettings?: AdvancedSettingsConfig;
  mode: ViewerMode;
  apiKey: string;
  currentSmiles: string;
  activeRunKey?: string;
  selectedStepId?: string;
  runs: Record<string, ViewerRun>;
  instanceSettings: Record<string, InstanceRequestSettings>;
  health: HealthStatusMap;
  editorOpen: boolean;
};
