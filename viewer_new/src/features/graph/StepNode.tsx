import type { NodeProps } from "@xyflow/react";
import { Beaker, FlaskConical, Pencil, RefreshCcw } from "lucide-react";

import { formatConfidence } from "../../lib/pathway";
import { MoleculePreview } from "./MoleculePreview";
import type { StepFlowNode } from "../../lib/layout";

export function StepNode({
  data,
  selected,
}: NodeProps<StepFlowNode>) {
  const stepNode = data.stepNode;
  const product = stepNode.products[0];

  return (
    <div
      className={`step-node${selected ? " selected" : ""}${stepNode.isVirtualRoot ? " root" : ""}${data.highlighted ? "" : " muted"}`}
      onClick={() => data.onInspect(stepNode.stepId)}
      role="button"
      tabIndex={0}
    >
      <div className="step-node__header">
        <div>
          <p className="eyebrow">{stepNode.isVirtualRoot ? "Target" : "Synthesis step"}</p>
          <h3>{stepNode.title}</h3>
        </div>
        <div className="badge-row">
          {stepNode.metrics?.confidenceestimate !== undefined ? (
            <span className="badge">
              <Beaker size={12} />
              {formatConfidence(stepNode.metrics.confidenceestimate)}
            </span>
          ) : null}
          {stepNode.metrics?.scalabilityindex !== undefined ? (
            <span className="badge subtle">
              <FlaskConical size={12} />
              Scale {stepNode.metrics.scalabilityindex}
            </span>
          ) : null}
        </div>
      </div>

      {product ? (
        <div className="step-node__product">
          <MoleculePreview
            smiles={product.smiles}
            label={product.label}
            width={stepNode.isVirtualRoot ? 220 : 180}
            height={stepNode.isVirtualRoot ? 150 : 110}
            compact
          />
          <div>
            <strong>{product.label}</strong>
            <p>{product.smiles}</p>
          </div>
        </div>
      ) : (
        <div className="step-node__placeholder">No product data available</div>
      )}

      {!stepNode.isVirtualRoot ? (
        <>
          <div className="molecule-chip-group">
            {stepNode.reactants.slice(0, 4).map((reactant) => (
              <span className="molecule-chip" key={reactant.key}>
                {reactant.label}
              </span>
            ))}
            {stepNode.reactants.length > 4 ? (
              <span className="molecule-chip more">
                +{stepNode.reactants.length - 4} more reactants
              </span>
            ) : null}
          </div>

          <div className="step-node__actions">
            <button
              className="ghost-button"
              onClick={(event) => {
                event.stopPropagation();
                data.onInspect(stepNode.stepId);
              }}
              type="button"
            >
              Inspect
            </button>
            <button
              className="ghost-button"
              onClick={(event) => {
                event.stopPropagation();
                data.onEdit(stepNode.stepId);
              }}
              type="button"
            >
              <Pencil size={14} />
              Edit
            </button>
            <button
              className="ghost-button"
              onClick={(event) => {
                event.stopPropagation();
                data.onPartialRerun(stepNode.stepId);
              }}
              type="button"
            >
              <RefreshCcw size={14} />
              Partial rerun
            </button>
          </div>
        </>
      ) : null}
    </div>
  );
}
