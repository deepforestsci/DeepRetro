import { Handle, Position, type NodeProps } from "@xyflow/react";
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
  const primaryFormula = product?.formula ?? product?.label;

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

      <div className="step-node__overlay">
        <div className="step-node__overlay-header">
          <div>
            <p className="eyebrow">{stepNode.isVirtualRoot ? "Target" : "Synthesis step"}</p>
            <h3>{stepNode.isVirtualRoot ? "Target" : `Step ${stepNode.stepId}`}</h3>
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

        <div className="step-node__overlay-meta">
          {primaryFormula ? (
            <span className="badge subtle">{primaryFormula}</span>
          ) : null}
          <span className="molecule-chip">{stepNode.products.length} product</span>
          <span className="molecule-chip">{stepNode.reactants.length} reactants</span>
          {stepNode.reagents.length ? (
            <span className="molecule-chip">{stepNode.reagents.length} reagents</span>
          ) : null}
        </div>

        {!stepNode.isVirtualRoot ? (
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
        ) : null}
      </div>
    </div>
  );
}
