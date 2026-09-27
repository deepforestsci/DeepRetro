import { Handle, Position, type NodeProps } from "@xyflow/react";
import { Pencil, RefreshCcw } from "lucide-react";

import { formatConfidence } from "../../lib/pathway";
import { MoleculePreview } from "./MoleculePreview";
import type { StepFlowNode } from "../../lib/layout";

function confidenceTier(value: number | undefined) {
  if (typeof value !== "number" || Number.isNaN(value)) {
    return "";
  }
  if (value >= 0.8) {
    return "ok";
  }
  if (value >= 0.5) {
    return "warn";
  }
  return "err";
}

export function StepNode({
  data,
  selected,
}: NodeProps<StepFlowNode>) {
  const stepNode = data.stepNode;
  const product = stepNode.products[0];
  const primaryFormula = product?.formula ?? product?.label;
  const confidence = stepNode.metrics?.confidenceestimate;
  const tier = confidenceTier(confidence);

  const counts = [
    `${stepNode.reactants.length} reactant${stepNode.reactants.length === 1 ? "" : "s"}`,
    stepNode.reagents.length
      ? `${stepNode.reagents.length} reagent${stepNode.reagents.length === 1 ? "" : "s"}`
      : "",
  ]
    .filter(Boolean)
    .join(" · ");

  return (
    <div
      className={`step-node${selected ? " selected" : ""}${stepNode.isVirtualRoot ? " root" : ""}${data.highlighted ? "" : " muted"}`}
      onClick={() => data.onInspect(stepNode.stepId)}
      role="button"
      tabIndex={0}
      style={{
        width: `${data.cardWidth}px`,
        minHeight: `${data.cardHeight}px`,
      }}
    >
      <Handle
        className="step-node__handle"
        isConnectable={false}
        position={Position.Left}
        type="target"
      />
      <Handle
        className="step-node__handle"
        isConnectable={false}
        position={Position.Right}
        type="source"
      />

      <div className="step-node__head">
        <span className="step-node__tag">
          {stepNode.isVirtualRoot ? "Target" : `Step ${stepNode.stepId}`}
        </span>
        {stepNode.isVirtualRoot ? (
          <span className="step-node__conf">input</span>
        ) : (
          <span className={`step-node__conf ${tier}`}>
            {formatConfidence(confidence)}
          </span>
        )}
      </div>

      {product ? (
        <div className="step-node__product solo">
          <MoleculePreview
            smiles={product.smiles}
            label={product.label}
            width={data.previewWidth}
            height={data.previewHeight}
            compact
          />
        </div>
      ) : (
        <div className="step-node__placeholder">No product data available</div>
      )}

      <div className="step-node__foot">
        <span className="step-node__formula">{primaryFormula ?? "—"}</span>
        <span className="step-node__count">
          {stepNode.isVirtualRoot ? product?.label ?? "" : counts}
        </span>
      </div>

      <div className="step-node__rail">
        <i
          className={tier}
          style={{
            width: stepNode.isVirtualRoot
              ? "100%"
              : typeof confidence === "number" && !Number.isNaN(confidence)
                ? `${Math.round(Math.min(Math.max(confidence, 0), 1) * 100)}%`
                : "0%",
          }}
        />
      </div>

      {!stepNode.isVirtualRoot ? (
        <div className="step-node__quick">
          <button
            aria-label={`Edit step ${stepNode.stepId}`}
            title="Edit step"
            onClick={(event) => {
              event.stopPropagation();
              data.onEdit(stepNode.stepId);
            }}
            type="button"
          >
            <Pencil size={13} />
          </button>
          <button
            aria-label={`Partial rerun from step ${stepNode.stepId}`}
            title="Partial rerun"
            onClick={(event) => {
              event.stopPropagation();
              data.onPartialRerun(stepNode.stepId);
            }}
            type="button"
          >
            <RefreshCcw size={13} />
          </button>
        </div>
      ) : null}
    </div>
  );
}
