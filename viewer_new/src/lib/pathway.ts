import { ZodError } from "zod";

import { pathwayResultSchema, pathwayStepSchema } from "./schemas";
import type {
  MoleculeMetadata,
  MoleculeRecord,
  NormalizedMolecule,
  NormalizedPathwayGraph,
  NormalizedStepNode,
  PathwayResult,
  PathwayStep,
  ValidationIssue,
} from "../types/viewer";

const VIRTUAL_ROOT_STEP = "0";

function getMetadata(
  molecule: MoleculeRecord,
  role: "product" | "reactant" | "reagent",
): MoleculeMetadata | undefined {
  if (role === "product") {
    return molecule.product_metadata;
  }
  if (role === "reactant") {
    return molecule.reactant_metadata;
  }
  return molecule.reagent_metadata;
}

function normalizeMolecules(
  molecules: MoleculeRecord[] | undefined,
  role: "product" | "reactant" | "reagent",
  stepId: string,
): NormalizedMolecule[] {
  return (molecules ?? []).map((molecule, index) => {
    const metadata = getMetadata(molecule, role);
    const label = metadata?.name || metadata?.chemical_formula || molecule.smiles;

    return {
      key: `${stepId}-${role}-${index}`,
      role,
      smiles: molecule.smiles,
      metadata,
      label,
      formula: metadata?.chemical_formula,
      mass: metadata?.mass,
    };
  });
}

function normalizeDependencyMap(
  result: PathwayResult,
): Record<string, string[]> {
  const dependencyMap: Record<string, string[]> = {};
  for (const [parent, children] of Object.entries(result.dependencies ?? {})) {
    dependencyMap[String(parent)] = (children ?? []).map((child) => String(child));
  }

  if (result.steps.some((step) => String(step.step) === "1")) {
    dependencyMap[VIRTUAL_ROOT_STEP] = ["1"];
  } else {
    dependencyMap[VIRTUAL_ROOT_STEP] = [];
  }

  return dependencyMap;
}

function buildParentMap(
  dependencyMap: Record<string, string[]>,
): Record<string, string[]> {
  const parentMap: Record<string, string[]> = {};

  for (const [parent, children] of Object.entries(dependencyMap)) {
    for (const child of children) {
      if (!parentMap[child]) {
        parentMap[child] = [];
      }
      parentMap[child].push(parent);
    }
  }

  return parentMap;
}

function createVirtualRootNode(result: PathwayResult): NormalizedStepNode {
  const primaryProduct = result.steps[0]?.products?.[0];
  const products = primaryProduct
    ? normalizeMolecules([primaryProduct], "product", VIRTUAL_ROOT_STEP)
    : [];

  return {
    id: `step-${VIRTUAL_ROOT_STEP}`,
    stepId: VIRTUAL_ROOT_STEP,
    title: "Target Molecule",
    isVirtualRoot: true,
    products,
    reactants: [],
    reagents: [],
    childIds: [],
    parentIds: [],
    conditions: result.steps[0]?.conditions,
    metrics: result.steps[0]?.reactionmetrics?.[0],
  };
}

export function buildPathwayGraph(result: PathwayResult): NormalizedPathwayGraph {
  const parsed = pathwayResultSchema.parse(result);
  const dependencyMap = normalizeDependencyMap(parsed);
  const parentMap = buildParentMap(dependencyMap);
  const stepMap = Object.fromEntries(parsed.steps.map((step) => [String(step.step), step]));
  const virtualRoot = createVirtualRootNode(parsed);

  const nodes: NormalizedStepNode[] = [virtualRoot];
  const edges: NormalizedPathwayGraph["edges"] = [];

  const sortedSteps = [...parsed.steps].sort(
    (left, right) => Number(left.step) - Number(right.step),
  );

  for (const step of sortedSteps) {
    const stepId = String(step.step);
    nodes.push({
      id: `step-${stepId}`,
      stepId,
      title: `Step ${stepId}`,
      isVirtualRoot: false,
      rawStep: step,
      products: normalizeMolecules(step.products, "product", stepId),
      reactants: normalizeMolecules(step.reactants, "reactant", stepId),
      reagents: normalizeMolecules(step.reagents, "reagent", stepId),
      metrics: step.reactionmetrics?.[0],
      conditions: step.conditions,
      childIds: dependencyMap[stepId] ?? [],
      parentIds: parentMap[stepId] ?? [],
    });
  }

  nodes[0].childIds = dependencyMap[VIRTUAL_ROOT_STEP] ?? [];
  const validNodeIds = new Set(nodes.map((node) => node.id));

  for (const [parent, children] of Object.entries(dependencyMap)) {
    for (const child of children) {
      if (parent === child) {
        continue;
      }
      const source = `step-${parent}`;
      const target = `step-${child}`;
      if (!validNodeIds.has(source) || !validNodeIds.has(target)) {
        continue;
      }
      edges.push({
        id: `edge-${parent}-${child}`,
        source,
        target,
      });
    }
  }

  return {
    stepMap,
    dependencyMap,
    nodes,
    edges,
    virtualRoot,
  };
}

export function validatePathwayResult(input: unknown) {
  try {
    return {
      data: pathwayResultSchema.parse(input),
      issues: [] as ValidationIssue[],
    };
  } catch (error) {
    if (error instanceof ZodError) {
      return {
        data: undefined,
        issues: error.issues.map((issue) => ({
          path: issue.path.join("."),
          message: issue.message,
        })),
      };
    }
    return {
      data: undefined,
      issues: [{ path: "", message: "Unknown validation error" }],
    };
  }
}

export function validatePathwayStep(input: unknown) {
  try {
    return {
      data: pathwayStepSchema.parse(input),
      issues: [] as ValidationIssue[],
    };
  } catch (error) {
    if (error instanceof ZodError) {
      return {
        data: undefined,
        issues: error.issues.map((issue) => ({
          path: issue.path.join("."),
          message: issue.message,
        })),
      };
    }
    return {
      data: undefined,
      issues: [{ path: "", message: "Unknown validation error" }],
    };
  }
}

export function getStep(result: PathwayResult | undefined, stepId: string) {
  return result?.steps.find((step) => String(step.step) === stepId);
}

export function updateStepInResult(
  result: PathwayResult,
  stepId: string,
  nextStep: PathwayStep,
): PathwayResult {
  return {
    ...result,
    steps: result.steps.map((step) =>
      String(step.step) === stepId ? nextStep : step,
    ),
  };
}

export function canPartialRerunStep(step: PathwayStep | undefined) {
  if (!step) {
    return {
      enabled: false,
      reason: "Select a valid step first.",
    };
  }

  const products = step.products ?? [];
  if (products.length !== 1) {
    return {
      enabled: false,
      reason: "Partial rerun requires exactly one product SMILES on the selected step.",
    };
  }

  if (!products[0]?.smiles) {
    return {
      enabled: false,
      reason: "The selected step is missing a product SMILES string.",
    };
  }

  return {
    enabled: true,
    reason: "",
  };
}

export function prettyValidationIssues(issues: ValidationIssue[]) {
  return issues.map((issue) =>
    issue.path ? `${issue.path}: ${issue.message}` : issue.message,
  );
}

export function inferTopLevelSmiles(result: PathwayResult | undefined) {
  return result?.steps[0]?.products?.[0]?.smiles ?? "";
}

export function formatConfidence(value: number | undefined) {
  if (typeof value !== "number" || Number.isNaN(value)) {
    return "n/a";
  }
  return `${Math.round(value * 100)}%`;
}

export function safeJsonParse(text: string) {
  try {
    return {
      data: JSON.parse(text) as unknown,
      error: "",
    };
  } catch (error) {
    return {
      data: undefined,
      error: error instanceof Error ? error.message : "Invalid JSON",
    };
  }
}

export function stringifyResult(result: PathwayResult | undefined) {
  return JSON.stringify(result ?? { steps: [], dependencies: {} }, null, 2);
}
