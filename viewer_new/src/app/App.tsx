import { useMemo, useState } from "react";
import { Download, KeyRound, Moon, Play, RefreshCcw, Save, Sun, Upload } from "lucide-react";

import { GraphCanvas } from "../features/graph/GraphCanvas";
import { InspectorPanel } from "../features/inspector/InspectorPanel";
import { RunSidebar } from "../features/runs/RunSidebar";
import { UploadPanel } from "../features/upload/UploadPanel";
import { EditorDrawer } from "../features/editor/EditorDrawer";
import { useBootstrapViewer, useViewerActions } from "../features/runs/useViewerActions";
import { useTheme } from "../features/theme/useTheme";
import { useViewerStore } from "../features/runs/store";
import { canPartialRerunStep, getStep } from "../lib/pathway";
import type { HealthStatusMap } from "../types/viewer";

function summarizeHealthStates(states: HealthStatusMap) {
  const values = Object.values(states);
  if (!values.length) {
    return { label: "Health unknown", tone: "muted" } as const;
  }
  if (values.every((value) => value.state === "healthy")) {
    return { label: "All backends healthy", tone: "ok" } as const;
  }
  if (values.some((value) => value.state === "healthy")) {
    return { label: "Backend health degraded", tone: "warn" } as const;
  }
  return { label: "Backends unavailable", tone: "err" } as const;
}

function BrandMark() {
  return (
    <svg
      aria-hidden="true"
      fill="none"
      height="26"
      stroke="currentColor"
      strokeWidth="1.8"
      viewBox="0 0 24 24"
      width="26"
    >
      <path d="M12 2.5l8.2 4.75v9.5L12 21.5l-8.2-4.75v-9.5z" />
      <path
        d="M12 7.2v9.6M7.8 9.6l8.4 4.8M16.2 9.6l-8.4 4.8"
        opacity="0.55"
        strokeWidth="1.2"
      />
    </svg>
  );
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

  const handleUploadedFiles = (files: Array<{ fileName: string; input: unknown }>) => {
    const loadedRunKeys = files.map(({ fileName, input }) =>
      loadUploadedResult(fileName, input),
    );
    const latestRunKey = loadedRunKeys[loadedRunKeys.length - 1];
    if (latestRunKey) {
      setActiveRun(latestRunKey);
    }
    setFeedback(
      `Loaded ${files.length} file${files.length === 1 ? "" : "s"} into the new viewer.`,
    );
  };

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
  const { theme, toggleTheme } = useTheme();

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
            <BrandMark />
          </div>
          <div>
            <h1 className="brand-name">DeepRetro</h1>
            <p className="brand-sub">Retrosynthesis</p>
          </div>
        </div>

        <div className="segmented-control" role="tablist" aria-label="viewer mode">
          <button
            className={mode === "analyze" ? "active" : ""}
            onClick={() => setMode("analyze")}
            type="button"
          >
            <Play size={13} />
            Analyze
          </button>
          <button
            className={mode === "upload" ? "active" : ""}
            onClick={() => setMode("upload")}
            type="button"
          >
            <Upload size={13} />
            Upload
          </button>
        </div>

        <label className="smiles-field">
          <span>SMILES</span>
          <input
            placeholder="Target molecule SMILES"
            spellCheck={false}
            value={currentSmiles}
            onChange={(event) => setCurrentSmiles(event.target.value)}
          />
        </label>

        <label className="api-key-field">
          <KeyRound size={13} />
          <input
            placeholder="API key"
            type="password"
            value={apiKey}
            onChange={(event) => setApiKey(event.target.value)}
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
          <Play size={13} />
          {isAnalyzing ? "Running..." : "Run analysis"}
        </button>

        <button
          aria-label="Run again"
          className="ghost-button icon-button"
          disabled={mode !== "analyze" || isAnalyzing || !currentSmiles.trim()}
          title="Run again"
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
          <RefreshCcw size={15} />
        </button>

        <button
          aria-label="Download JSON"
          className="ghost-button icon-button"
          disabled={!activeRun?.rawResult}
          title="Download JSON"
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
          <Download size={15} />
        </button>

        <button
          aria-label="Save edits"
          className="ghost-button icon-button"
          disabled={saveDisabled}
          title="Save edits"
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
          <Save size={15} />
        </button>

        <div className={`health-chip ${healthSummary.tone} topbar-spacer`}>
          <span className="dot" />
          {healthSummary.label}
        </div>

        <button
          aria-label={theme === "light" ? "Switch to dark theme" : "Switch to light theme"}
          className="ghost-button icon-button"
          title={theme === "light" ? "Switch to dark theme" : "Switch to light theme"}
          onClick={toggleTheme}
          type="button"
        >
          {theme === "light" ? <Moon size={15} /> : <Sun size={15} />}
        </button>
      </header>

      {feedback || error ? (
        <div className="status-strip">
          {feedback ? <span className="status-pill success">{feedback}</span> : null}
          {error ? <span className="status-pill error">{error}</span> : null}
        </div>
      ) : (
        <div />
      )}

      <div className="workspace">
        <RunSidebar
          runtimeConfig={runtimeConfig}
          advancedSettings={advancedSettings}
          health={health}
          instanceSettings={instanceSettings}
          runs={runs}
          activeRunKey={activeRunKey}
          onRunSelect={setActiveRun}
          onUploadFiles={(files) => {
            try {
              setError("");
              setFeedback("");
              handleUploadedFiles(files);
            } catch (nextError) {
              setError(nextError instanceof Error ? nextError.message : "Upload failed.");
            }
          }}
          onUpdateSetting={updateInstanceSetting}
        />

        <section className="main-column">
          {mode === "upload" ? (
            <UploadPanel
              onLoadFiles={(files) => {
                try {
                  setError("");
                  setFeedback("");
                  handleUploadedFiles(files);
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
          onSelectStep={setSelectedStep}
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
