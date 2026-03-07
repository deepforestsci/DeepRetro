import { describe, expect, test } from "vitest";

import {
  buildPathwayGraph,
  canPartialRerunStep,
  updateStepInResult,
  validatePathwayResult,
} from "./pathway";
import type { PathwayResult } from "../types/viewer";

const sampleResult: PathwayResult = {
  dependencies: {
    1: ["2"],
    2: [],
  },
  steps: [
    {
      step: "1",
      reactants: [{ smiles: "CCO", reactant_metadata: { chemical_formula: "C2H6O" } }],
      products: [{ smiles: "CC=O", product_metadata: { chemical_formula: "C2H4O" } }],
      reactionmetrics: [{ confidenceestimate: 0.88, scalabilityindex: "7" }],
      conditions: { solvent: "MeOH" },
    },
    {
      step: "2",
      reactants: [{ smiles: "O", reactant_metadata: { chemical_formula: "H2O" } }],
      products: [{ smiles: "CCO", product_metadata: { chemical_formula: "C2H6O" } }],
    },
  ],
};

describe("pathway utilities", () => {
  test("builds a graph with a virtual root", () => {
    const graph = buildPathwayGraph(sampleResult);

    expect(graph.virtualRoot.stepId).toBe("0");
    expect(graph.nodes).toHaveLength(3);
    expect(graph.edges.map((edge) => edge.id)).toContain("edge-0-1");
    expect(graph.edges.map((edge) => edge.id)).toContain("edge-1-2");
  });

  test("updates a step in place without mutating other steps", () => {
    const next = updateStepInResult(sampleResult, "2", {
      ...sampleResult.steps[1],
      products: [{ smiles: "CO", product_metadata: { chemical_formula: "CH4O" } }],
    });

    expect(next.steps[1].products?.[0]?.smiles).toBe("CO");
    expect(sampleResult.steps[1].products?.[0]?.smiles).toBe("CCO");
  });

  test("validates partial rerun eligibility from product count", () => {
    expect(canPartialRerunStep(sampleResult.steps[0]).enabled).toBe(true);
    expect(
      canPartialRerunStep({
        ...sampleResult.steps[0],
        products: [],
      }).enabled,
    ).toBe(false);
  });

  test("returns validation issues for malformed pathway payloads", () => {
    const validation = validatePathwayResult({ steps: {} });

    expect(validation.data).toBeUndefined();
    expect(validation.issues.length).toBeGreaterThan(0);
  });
});
