import { useMemo, useState } from "react";
import { Download, KeyRound, Play, RefreshCcw, Save, Sparkles, Upload } from "lucide-react";

import { GraphCanvas } from "../features/graph/GraphCanvas";
import { InspectorPanel } from "../features/inspector/InspectorPanel";
import { RunSidebar } from "../features/runs/RunSidebar";
import { UploadPanel } from "../features/upload/UploadPanel";
import { EditorDrawer } from "../features/editor/EditorDrawer";
import { useBootstrapViewer, useViewerActions } from "../features/runs/useViewerActions";
import { useViewerStore } from "../features/runs/store";
import { canPartialRerunStep, getStep } from "../lib/pathway";
import type { HealthStatusMap } from "../types/viewer";

function summarizeHealthStates(states: HealthStatusMap) {
  const values = Object.values(states);
  if (!values.length) {
    return "Health unknown";
  }
  if (values.every((value) => value.state === "healthy")) {
    return "All backends healthy";
  }
  if (values.some((value) => value.state === "healthy")) {
    return "Backend health degraded";
  }
  return "Backends unavailable";
}

export function App() {
  useBootstrapViewer();

  const runtimeConfig = useViewerStore((state) => state.runtimeConfig);
  const advancedSettings = useViewerStore((state) => state.advancedSettings);
  const health = useViewerStore((state) => state.health);
  const runs = useViewerStore((state) => state.runs);
  const instanceSettings = useViewerStore((state) => state.instanceSettings);
  const mode = useViewerStore((state) => state.mode);
  const apiKey = useViewerStore((state) => state.apiKey);
  const currentSmiles = useViewerStore((state) => state.currentSmiles);
  const activeRunKey = useViewerStore((state) => state.activeRunKey);
  const selectedStepId = useViewerStore((state) => state.selectedStepId);
  const editorOpen = useViewerStore((state) => state.editorOpen);
  const setMode = useViewerStore((state) => state.setMode);
  const setApiKey = useViewerStore((state) => state.setApiKey);
  const setCurrentSmiles = useViewerStore((state) => state.setCurrentSmiles);
  const setActiveRun = useViewerStore((state) => state.setActiveRun);
  const setSelectedStep = useViewerStore((state) => state.setSelectedStep);
  const updateInstanceSetting = useViewerStore((state) => state.updateInstanceSetting);
  const setEditorOpen = useViewerStore((state) => state.setEditorOpen);

  const [feedback, setFeedback] = useState("");
  const [error, setError] = useState("");

  const {
    activeRun,
    isAnalyzing,
    isSaving,
    isPartialRerunning,
    analyzeAll,
    loadUploadedResult,
    applyStructuredStepEdit,
    applyRawJsonEdit,
    saveRunEdits,
    partialRerun,
  } = useViewerActions();

  const selectedStep = useMemo(
    () =>
      activeRun?.rawResult && selectedStepId
        ? getStep(activeRun.rawResult, selectedStepId)
        : undefined,
    [activeRun?.rawResult, selectedStepId],
  );

  const partialRerunState = canPartialRerunStep(selectedStep);
  const saveDisabled = !activeRun || activeRun.source !== "api" || !activeRun.dirty || isSaving;
  const healthSummary = summarizeHealthStates(health);

  if (!runtimeConfig || !advancedSettings) {
    return (
      <main className="boot-screen">
        <div className="panel boot-card">
          <p className="eyebrow">DeepRetro</p>
          <h1>Loading `viewer_new`</h1>
          <p>Fetching runtime config, advanced settings, and backend capabilities.</p>
        </div>
      </main>
    );
  }

  return (
    <main className="app-shell">
      <header className="topbar">
        <div className="brand">
          <div className="brand-mark">
            <Sparkles size={18} />
          </div>
          <div>
            <p className="eyebrow">DeepRetro Viewer</p>
            <h1>Retrosynthesis Workbench</h1>
          </div>
        </div>

        <div className="topbar__controls">
          <div className="segmented-control" role="tablist" aria-label="viewer mode">
            <button
              className={mode === "analyze" ? "active" : ""}
              onClick={() => setMode("analyze")}
              type="button"
            >
              <Play size={14} />
              Analyze
            </button>
            <button
              className={mode === "upload" ? "active" : ""}
              onClick={() => setMode("upload")}
              type="button"
            >
              <Upload size={14} />
              Upload
            </button>
          </div>

          <label className="api-key-field">
            <KeyRound size={14} />
            <input
              placeholder="Backend API key"
              type="password"
              value={apiKey}
              onChange={(event) => setApiKey(event.target.value)}
            />
          </label>

          <label className="smiles-field">
            <span>SMILES</span>
            <input
              placeholder="Enter the target molecule SMILES"
              value={currentSmiles}
              onChange={(event) => setCurrentSmiles(event.target.value)}
            />
          </label>

          <button
            className="primary-button"
            disabled={mode !== "analyze" || isAnalyzing || !currentSmiles.trim()}
            onClick={async () => {
              try {
                setError("");
                setFeedback("");
                await analyzeAll(currentSmiles.trim(), false);
                setFeedback("Analysis completed. Switch between pathway runs from the sidebar.");
              } catch (nextError) {
                setError(
                  nextError instanceof Error ? nextError.message : "Analysis failed.",
                );
              }
            }}
            type="button"
          >
            <Play size={14} />
            {isAnalyzing ? "Running..." : "Run analysis"}
          </button>

          <button
            className="ghost-button"
            disabled={mode !== "analyze" || isAnalyzing || !currentSmiles.trim()}
            onClick={async () => {
              try {
                setError("");
                setFeedback("");
                await analyzeAll(currentSmiles.trim(), true);
                setFeedback("Rerun completed using the current instance settings.");
              } catch (nextError) {
                setError(nextError instanceof Error ? nextError.message : "Rerun failed.");
              }
            }}
            type="button"
          >
            <RefreshCcw size={14} />
            Run again
          </button>

          <button
            className="ghost-button"
            disabled={!activeRun?.rawResult}
            onClick={() => {
              if (!activeRun?.rawResult) {
                return;
              }
              const blob = new Blob([JSON.stringify(activeRun.rawResult, null, 2)], {
                type: "application/json",
              });
              const url = URL.createObjectURL(blob);
              const anchor = document.createElement("a");
              anchor.href = url;
              anchor.download = `${activeRun.label.replace(/\s+/g, "-").toLowerCase()}.json`;
              anchor.click();
              URL.revokeObjectURL(url);
            }}
            type="button"
          >
            <Download size={14} />
            Download JSON
          </button>

          <button
            className="ghost-button"
            disabled={saveDisabled}
            onClick={async () => {
              if (!activeRunKey) {
                return;
              }
              try {
                setError("");
                setFeedback("");
                await saveRunEdits(activeRunKey);
                setFeedback("Edits saved to the backend. Future partial reruns will use the updated pathway.");
              } catch (nextError) {
                setError(nextError instanceof Error ? nextError.message : "Save failed.");
              }
            }}
            type="button"
          >
            <Save size={14} />
            Save edits
          </button>
        </div>

        <div className="topbar__status">
          <span className="status-pill healthy">{healthSummary}</span>
          {feedback ? <span className="status-pill success">{feedback}</span> : null}
          {error ? <span className="status-pill error">{error}</span> : null}
        </div>
      </header>

      <div className="workspace">
        <RunSidebar
          runtimeConfig={runtimeConfig}
          advancedSettings={advancedSettings}
          health={health}
          instanceSettings={instanceSettings}
          runs={runs}
          activeRunKey={activeRunKey}
          onRunSelect={setActiveRun}
          onUpdateSetting={updateInstanceSetting}
        />

        <section className="main-column">
          {mode === "upload" ? (
            <UploadPanel
              onLoad={(fileName, input) => {
                try {
                  setError("");
                  setFeedback("");
                  const runKey = loadUploadedResult(fileName, input);
                  setMode("analyze");
                  setActiveRun(runKey);
                  setFeedback(`Loaded ${fileName} into the new viewer.`);
                } catch (nextError) {
                  setError(nextError instanceof Error ? nextError.message : "Upload failed.");
                }
              }}
            />
          ) : null}

          <section className="panel graph-panel">
            <div className="panel__header">
              <div>
                <p className="eyebrow">Pathway graph</p>
                <h2>{activeRun?.label ?? "No run selected"}</h2>
              </div>
              {activeRun?.dirty ? <span className="status-pill dirty">local edits pending</span> : null}
            </div>

            <GraphCanvas
              run={activeRun}
              selectedStepId={selectedStepId}
              onSelectStep={(stepId) => setSelectedStep(stepId)}
              onEditStep={(stepId) => {
                setSelectedStep(stepId);
                setEditorOpen(true);
              }}
              onPartialRerunStep={async (stepId) => {
                if (!activeRunKey) {
                  return;
                }

                try {
                  setError("");
                  setFeedback("");
                  await partialRerun(activeRunKey, stepId);
                  setFeedback(`Partial rerun completed from step ${stepId}.`);
                } catch (nextError) {
                  setError(
                    nextError instanceof Error ? nextError.message : "Partial rerun failed.",
                  );
                }
              }}
            />
          </section>
        </section>

        <InspectorPanel
          run={activeRun}
          selectedStepId={selectedStepId}
          onEdit={() => setEditorOpen(true)}
          onPartialRerun={async () => {
            if (!activeRunKey || !selectedStepId || !partialRerunState.enabled) {
              setError(partialRerunState.reason || "Partial rerun is not available.");
              return;
            }

            try {
              setError("");
              setFeedback("");
              await partialRerun(activeRunKey, selectedStepId);
              setFeedback(`Partial rerun completed from step ${selectedStepId}.`);
            } catch (nextError) {
              setError(
                nextError instanceof Error ? nextError.message : "Partial rerun failed.",
              );
            }
          }}
          onSaveEdits={async () => {
            if (!activeRunKey) {
              return;
            }
            try {
              setError("");
              setFeedback("");
              await saveRunEdits(activeRunKey);
              setFeedback("Edits saved to the backend.");
            } catch (nextError) {
              setError(nextError instanceof Error ? nextError.message : "Save failed.");
            }
          }}
          saveDisabled={saveDisabled}
        />
      </div>

      <EditorDrawer
        open={editorOpen}
        run={activeRun}
        stepId={selectedStepId && selectedStepId !== "0" ? selectedStepId : undefined}
        onClose={() => setEditorOpen(false)}
        onApplyStructured={(stepId, nextStep) => {
          if (!activeRunKey) {
            return;
          }
          applyStructuredStepEdit(activeRunKey, stepId, nextStep);
          setFeedback(`Applied structured changes to step ${stepId}.`);
        }}
        onApplyRaw={(input) => {
          if (!activeRunKey) {
            return;
          }
          applyRawJsonEdit(activeRunKey, input);
          setFeedback("Applied raw pathway JSON changes.");
        }}
      />

      {(isPartialRerunning || isSaving) && (
        <div className="busy-overlay">
          <div className="panel busy-card">
            <p className="eyebrow">Working</p>
            <h2>{isPartialRerunning ? "Running partial synthesis" : "Saving pathway edits"}</h2>
            <p>The new viewer is waiting for the backend to finish.</p>
          </div>
        </div>
      )}
    </main>
  );
}
