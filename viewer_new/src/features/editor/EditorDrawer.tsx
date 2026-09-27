import { useEffect, useMemo, useState } from "react";
import CodeMirror from "@uiw/react-codemirror";
import { json as jsonLanguage } from "@codemirror/lang-json";
import { X } from "lucide-react";

import { getStep, safeJsonParse, stringifyResult } from "../../lib/pathway";
import type { MoleculeRecord, PathwayStep, ViewerRun } from "../../types/viewer";

type EditorDrawerProps = {
  open: boolean;
  run?: ViewerRun;
  stepId?: string;
  onClose: () => void;
  onApplyStructured: (stepId: string, nextStep: PathwayStep) => void;
  onApplyRaw: (input: unknown) => void;
};

function cloneStep(step: PathwayStep | undefined): PathwayStep {
  if (!step) {
    return {
      step: "",
      products: [],
      reactants: [],
      reagents: [],
      conditions: {},
      reactionmetrics: [{}],
    };
  }

  return JSON.parse(JSON.stringify(step)) as PathwayStep;
}

function makeBlankMolecule(): MoleculeRecord {
  return {
    smiles: "",
    product_metadata: {},
  };
}

function serializeConditionValue(value: unknown) {
  if (typeof value === "string") {
    return value;
  }
  if (value === undefined) {
    return "";
  }
  return JSON.stringify(value);
}

function parseConditionValue(value: string) {
  const trimmed = value.trim();
  if (!trimmed) {
    return "";
  }

  const shouldParseAsJson =
    trimmed.startsWith("{") ||
    trimmed.startsWith("[") ||
    trimmed.startsWith('"') ||
    trimmed === "true" ||
    trimmed === "false" ||
    trimmed === "null" ||
    /^-?\d+(\.\d+)?$/.test(trimmed);

  if (!shouldParseAsJson) {
    return value;
  }

  const parsed = safeJsonParse(trimmed);
  if (parsed.data !== undefined) {
    return parsed.data;
  }

  return value;
}

function ConditionEditor({
  conditions,
  onChange,
}: {
  conditions: PathwayStep["conditions"];
  onChange: (next: PathwayStep["conditions"]) => void;
}) {
  const entries = Array.isArray(conditions)
    ? []
    : Object.entries(conditions ?? {}).map(([key, value]) => ({
        key,
        value: serializeConditionValue(value),
      }));

  const updateEntries = (
    nextEntries: Array<{
      key: string;
      value: string;
    }>,
  ) => {
    const nextConditions: Record<string, unknown> = {};
    nextEntries.forEach((entry) => {
      if (!entry.key.trim()) {
        return;
      }
      nextConditions[entry.key] = parseConditionValue(entry.value);
    });
    onChange(nextConditions);
  };

  return (
    <section className="editor-section">
      <div className="section-heading">
        <h3>Conditions</h3>
        <button
          className="ghost-button"
          onClick={() =>
            updateEntries([
              ...entries,
              {
                key: "",
                value: "",
              },
            ])
          }
          type="button"
        >
          Add
        </button>
      </div>
      {Array.isArray(conditions) ? (
        <p className="inline-warning">
          This step stores conditions as an array. Switch to raw JSON to edit that payload safely.
        </p>
      ) : null}
      <div className="editor-condition-list">
        {entries.map((entry, index) => (
          <div className="editor-condition-row" key={`condition-${index}`}>
            <label>
              <span>Field</span>
              <input
                value={entry.key}
                onChange={(event) => {
                  const next = [...entries];
                  next[index] = {
                    ...next[index],
                    key: event.target.value,
                  };
                  updateEntries(next);
                }}
              />
            </label>
            <label>
              <span>Value</span>
              <input
                value={entry.value}
                onChange={(event) => {
                  const next = [...entries];
                  next[index] = {
                    ...next[index],
                    value: event.target.value,
                  };
                  updateEntries(next);
                }}
              />
            </label>
            <button
              className="ghost-button danger"
              onClick={() => updateEntries(entries.filter((_, candidate) => candidate !== index))}
              type="button"
            >
              Remove
            </button>
          </div>
        ))}
        {!entries.length && !Array.isArray(conditions) ? (
          <p className="empty-note">No structured conditions yet. Add fields above.</p>
        ) : null}
      </div>
    </section>
  );
}

function MoleculeEditor({
  label,
  molecules,
  metadataKey,
  onChange,
}: {
  label: string;
  molecules: MoleculeRecord[];
  metadataKey: "product_metadata" | "reactant_metadata" | "reagent_metadata";
  onChange: (next: MoleculeRecord[]) => void;
}) {
  return (
    <section className="editor-section">
      <div className="section-heading">
        <h3>{label}</h3>
        <button
          className="ghost-button"
          onClick={() =>
            onChange([
              ...molecules,
              {
                ...makeBlankMolecule(),
                [metadataKey]: {},
              },
            ])
          }
          type="button"
        >
          Add
        </button>
      </div>
      <div className="editor-molecule-list">
        {molecules.map((molecule, index) => {
          const metadata = molecule[metadataKey] ?? {};
          return (
            <div className="editor-molecule-row" key={`${label}-${index}`}>
              <label>
                <span>SMILES</span>
                <input
                  value={molecule.smiles}
                  onChange={(event) => {
                    const next = [...molecules];
                    next[index] = { ...next[index], smiles: event.target.value };
                    onChange(next);
                  }}
                />
              </label>
              <label>
                <span>Name</span>
                <input
                  value={String(metadata.name ?? "")}
                  onChange={(event) => {
                    const next = [...molecules];
                    next[index] = {
                      ...next[index],
                      [metadataKey]: {
                        ...metadata,
                        name: event.target.value,
                      },
                    };
                    onChange(next);
                  }}
                />
              </label>
              <label>
                <span>Formula</span>
                <input
                  value={String(metadata.chemical_formula ?? "")}
                  onChange={(event) => {
                    const next = [...molecules];
                    next[index] = {
                      ...next[index],
                      [metadataKey]: {
                        ...metadata,
                        chemical_formula: event.target.value,
                      },
                    };
                    onChange(next);
                  }}
                />
              </label>
              <label>
                <span>Mass</span>
                <input
                  type="number"
                  value={metadata.mass === undefined ? "" : String(metadata.mass)}
                  onChange={(event) => {
                    const next = [...molecules];
                    next[index] = {
                      ...next[index],
                      [metadataKey]: {
                        ...metadata,
                        mass: event.target.value ? Number(event.target.value) : undefined,
                      },
                    };
                    onChange(next);
                  }}
                />
              </label>
              <button
                className="ghost-button danger"
                onClick={() => onChange(molecules.filter((_, candidate) => candidate !== index))}
                type="button"
              >
                Remove
              </button>
            </div>
          );
        })}
      </div>
    </section>
  );
}

export function EditorDrawer({
  open,
  run,
  stepId,
  onClose,
  onApplyStructured,
  onApplyRaw,
}: EditorDrawerProps) {
  const selectedStep = useMemo(
    () => (run?.rawResult && stepId ? getStep(run.rawResult, stepId) : undefined),
    [run?.rawResult, stepId],
  );
  const [tab, setTab] = useState<"structured" | "raw">("structured");
  const [draft, setDraft] = useState<PathwayStep>(cloneStep(selectedStep));
  const [rawText, setRawText] = useState(stringifyResult(run?.rawResult));
  const [error, setError] = useState("");

  useEffect(() => {
    setDraft(cloneStep(selectedStep));
  }, [selectedStep]);

  useEffect(() => {
    setRawText(stringifyResult(run?.rawResult));
  }, [run?.rawResult]);

  if (!open || !run || !stepId) {
    return null;
  }

  return (
    <div className="editor-drawer__backdrop" onClick={onClose}>
      <aside
        className="editor-drawer"
        onClick={(event) => event.stopPropagation()}
      >
        <div className="panel__header">
          <div>
            <p className="eyebrow">{run.label}</p>
            <h2>Edit step {stepId}</h2>
          </div>
          <button className="ghost-button" onClick={onClose} type="button">
            <X size={16} />
          </button>
        </div>

        <div className="tab-row">
          <button
            className={tab === "structured" ? "active" : ""}
            onClick={() => setTab("structured")}
            type="button"
          >
            Structured
          </button>
          <button
            className={tab === "raw" ? "active" : ""}
            onClick={() => setTab("raw")}
            type="button"
          >
            Full pathway JSON
          </button>
        </div>

        {tab === "structured" ? (
          <div className="editor-drawer__content">
            <label>
              <span>Step number</span>
              <input
                value={draft.step}
                onChange={(event) =>
                  setDraft((current) => ({
                    ...current,
                    step: event.target.value,
                  }))
                }
              />
            </label>

            <MoleculeEditor
              label="Products"
              metadataKey="product_metadata"
              molecules={draft.products ?? []}
              onChange={(products) =>
                setDraft((current) => ({ ...current, products }))
              }
            />
            <MoleculeEditor
              label="Reactants"
              metadataKey="reactant_metadata"
              molecules={draft.reactants ?? []}
              onChange={(reactants) =>
                setDraft((current) => ({ ...current, reactants }))
              }
            />
            <MoleculeEditor
              label="Reagents"
              metadataKey="reagent_metadata"
              molecules={draft.reagents ?? []}
              onChange={(reagents) =>
                setDraft((current) => ({ ...current, reagents }))
              }
            />

            <ConditionEditor
              conditions={draft.conditions}
              onChange={(conditions) =>
                setDraft((current) => ({
                  ...current,
                  conditions,
                }))
              }
            />

            <section className="editor-section">
              <div className="section-heading">
                <h3>Reaction metrics</h3>
              </div>
              <div className="form-grid">
                <label>
                  <span>Scalability index</span>
                  <input
                    value={String(draft.reactionmetrics?.[0]?.scalabilityindex ?? "")}
                    onChange={(event) =>
                      setDraft((current) => ({
                        ...current,
                        reactionmetrics: [
                          {
                            ...(current.reactionmetrics?.[0] ?? {}),
                            scalabilityindex: event.target.value,
                          },
                        ],
                      }))
                    }
                  />
                </label>
                <label>
                  <span>Confidence estimate</span>
                  <input
                    type="number"
                    step="0.01"
                    value={String(draft.reactionmetrics?.[0]?.confidenceestimate ?? "")}
                    onChange={(event) =>
                      setDraft((current) => ({
                        ...current,
                        reactionmetrics: [
                          {
                            ...(current.reactionmetrics?.[0] ?? {}),
                            confidenceestimate: event.target.value
                              ? Number(event.target.value)
                              : undefined,
                          },
                        ],
                      }))
                    }
                  />
                </label>
                <label className="full-span">
                  <span>Closest literature</span>
                  <input
                    value={String(draft.reactionmetrics?.[0]?.closestliterature ?? "")}
                    onChange={(event) =>
                      setDraft((current) => ({
                        ...current,
                        reactionmetrics: [
                          {
                            ...(current.reactionmetrics?.[0] ?? {}),
                            closestliterature: event.target.value,
                          },
                        ],
                      }))
                    }
                  />
                </label>
              </div>
            </section>

            {error ? <p className="inline-error">{error}</p> : null}

            <div className="editor-drawer__footer">
              <button
                className="primary-button"
                onClick={() => {
                  try {
                    onApplyStructured(stepId, draft);
                    setError("");
                    onClose();
                  } catch (nextError) {
                    setError(
                      nextError instanceof Error ? nextError.message : "Failed to apply changes.",
                    );
                  }
                }}
                type="button"
              >
                Apply structured changes
              </button>
            </div>
          </div>
        ) : (
          <div className="editor-drawer__content">
            <CodeMirror
              value={rawText}
              extensions={[jsonLanguage()]}
              height="100%"
              onChange={setRawText}
              theme="dark"
            />
            {error ? <p className="inline-error">{error}</p> : null}
            <div className="editor-drawer__footer">
              <button
                className="primary-button"
                onClick={() => {
                  const parsed = safeJsonParse(rawText);
                  if (parsed.error || parsed.data === undefined) {
                    setError(parsed.error);
                    return;
                  }

                  try {
                    onApplyRaw(parsed.data);
                    setError("");
                    onClose();
                  } catch (nextError) {
                    setError(
                      nextError instanceof Error ? nextError.message : "Failed to apply JSON changes.",
                    );
                  }
                }}
                type="button"
              >
                Apply raw JSON
              </button>
            </div>
          </div>
        )}
      </aside>
    </div>
  );
}
