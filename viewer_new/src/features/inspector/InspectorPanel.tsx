import { ChevronDown, Pencil, RefreshCcw, Save } from "lucide-react";

import {
  canPartialRerunStep,
  formatConfidence,
  getStep,
  stringifyResult,
} from "../../lib/pathway";
import { MoleculePreview } from "../graph/MoleculePreview";
import type { ViewerRun } from "../../types/viewer";

type InspectorPanelProps = {
  run?: ViewerRun;
  selectedStepId?: string;
  onSelectStep: (stepId: string) => void;
  onEdit: () => void;
  onPartialRerun: () => void;
  onSaveEdits: () => void;
  saveDisabled: boolean;
};

function MoleculeSection({
  title,
  molecules,
}: {
  title: string;
  molecules: NonNullable<NonNullable<ViewerRun["graph"]>["nodes"][number]["products"]>;
}) {
  if (!molecules.length) {
    return null;
  }

  return (
    <section className="inspector-section">
      <div className="section-heading">
        <h3>{title}</h3>
        <span>{molecules.length}</span>
      </div>
      <div className="inspector-molecules">
        {molecules.map((molecule) => (
          <article className="inspector-molecule" key={molecule.key}>
            <MoleculePreview
              smiles={molecule.smiles}
              label={molecule.label}
              width={220}
              height={140}
            />
            <strong>{molecule.label}</strong>
            <p>{molecule.smiles}</p>
            <div className="badge-row">
              {molecule.formula ? <span className="badge">{molecule.formula}</span> : null}
              {typeof molecule.mass === "number" ? (
                <span className="badge subtle">{molecule.mass.toFixed(1)} g/mol</span>
              ) : null}
            </div>
          </article>
        ))}
      </div>
    </section>
  );
}

function ConditionCards({
  conditions,
}: {
  conditions: Record<string, unknown> | unknown[] | undefined;
}) {
  if (!conditions || Array.isArray(conditions)) {
    return (
      <div className="conditions-grid">
        <div className="condition-card">
          <span>Conditions</span>
          <strong>{conditions ? "Unstructured payload" : "Not provided"}</strong>
        </div>
      </div>
    );
  }

  const entries = Object.entries(conditions);
  if (!entries.length) {
    return (
      <div className="conditions-grid">
        <div className="condition-card">
          <span>Conditions</span>
          <strong>Not provided</strong>
        </div>
      </div>
    );
  }

  return (
    <div className="conditions-grid">
      {entries.map(([key, value]) => (
        <div className="condition-card" key={key}>
          <span>{key.replaceAll("_", " ")}</span>
          <strong>{typeof value === "string" ? value : JSON.stringify(value)}</strong>
        </div>
      ))}
    </div>
  );
}

export function InspectorPanel({
  run,
  selectedStepId,
  onSelectStep,
  onEdit,
  onPartialRerun,
  onSaveEdits,
  saveDisabled,
}: InspectorPanelProps) {
  if (!run?.graph) {
    return (
      <aside className="inspector">
        <section className="panel inspector-panel">
          <div className="panel__header">
            <div>
              <p className="eyebrow">Inspector</p>
              <h2>No selection</h2>
            </div>
          </div>
          <p>Select a pathway run to inspect its steps, molecules, and raw JSON.</p>
        </section>
      </aside>
    );
  }

  const activeStepId = selectedStepId ?? run.graph.virtualRoot.stepId;
  const stepNode = run.graph.nodes.find((node) => node.stepId === activeStepId);
  const rawStep = getStep(run.rawResult, activeStepId);
  const rerunState = canPartialRerunStep(rawStep);

  return (
    <aside className="inspector">
      <section className="panel inspector-panel">
        <div className="panel__header">
          <div>
            <p className="eyebrow">{run.label}</p>
            <h2>{stepNode?.title ?? "Selected step"}</h2>
          </div>
          <div className="badge-row">
            <span className={`status-pill ${run.status}`}>{run.status}</span>
            {run.dirty ? <span className="status-pill dirty">unsaved edits</span> : null}
          </div>
        </div>

        {stepNode ? (
          <>
            <div className="inspector-actions">
              {!stepNode.isVirtualRoot ? (
                <button className="ghost-button" onClick={onEdit} type="button">
                  <Pencil size={14} />
                  Edit step
                </button>
              ) : null}
              {!stepNode.isVirtualRoot ? (
                <button
                  className="ghost-button"
                  disabled={!rerunState.enabled}
                  onClick={onPartialRerun}
                  type="button"
                >
                  <RefreshCcw size={14} />
                  Partial rerun
                </button>
              ) : null}
              {run.source === "api" ? (
                <button
                  className="ghost-button"
                  disabled={saveDisabled}
                  onClick={onSaveEdits}
                  type="button"
                >
                  <Save size={14} />
                  Save edits
                </button>
              ) : null}
            </div>

            {!stepNode.isVirtualRoot && !rerunState.enabled ? (
              <p className="inline-warning">{rerunState.reason}</p>
            ) : null}

            <div className="inspector-grid">
              <div className="metric-card">
                <span>Confidence</span>
                <strong>{formatConfidence(stepNode.metrics?.confidenceestimate)}</strong>
              </div>
              <div className="metric-card">
                <span>Scalability</span>
                <strong>{stepNode.metrics?.scalabilityindex ?? "n/a"}</strong>
              </div>
              <div className="metric-card">
                <span>Parents</span>
                <strong>{stepNode.parentIds.length}</strong>
              </div>
              <div className="metric-card">
                <span>Children</span>
                <strong>{stepNode.childIds.length}</strong>
              </div>
            </div>

            <section className="inspector-section">
              <div className="section-heading">
                <h3>Steps</h3>
              </div>
              <div className="molecule-chip-group">
                {run.graph.nodes.map((node) => (
                  <button
                    className={`molecule-chip action-chip${node.stepId === activeStepId ? " active" : ""}`}
                    key={node.id}
                    onClick={() => onSelectStep(node.stepId)}
                    type="button"
                  >
                    {node.isVirtualRoot ? "Target" : `Step ${node.stepId}`}
                  </button>
                ))}
              </div>
            </section>

            <MoleculeSection title="Products" molecules={stepNode.products} />
            <MoleculeSection title="Reactants" molecules={stepNode.reactants} />
            <MoleculeSection title="Reagents" molecules={stepNode.reagents} />

            <section className="inspector-section">
              <div className="section-heading">
                <h3>Conditions</h3>
              </div>
              <ConditionCards conditions={stepNode.conditions} />
            </section>

            <details className="json-details">
              <summary>
                <span>Raw step JSON</span>
                <ChevronDown size={14} />
              </summary>
              <pre className="json-block">
                {JSON.stringify(rawStep ?? stepNode.rawStep ?? {}, null, 2)}
              </pre>
            </details>

            <details className="json-details">
              <summary>
                <span>Run JSON</span>
                <ChevronDown size={14} />
              </summary>
              <pre className="json-block tall">{stringifyResult(run.rawResult)}</pre>
            </details>
          </>
        ) : (
          <p>The selected step is not available in the active run.</p>
        )}
      </section>
    </aside>
  );
}
