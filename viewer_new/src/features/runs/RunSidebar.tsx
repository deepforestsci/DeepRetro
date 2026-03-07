import { useId } from "react";
import { Activity, CheckCircle2, CircleAlert, LoaderCircle, Upload } from "lucide-react";

import type {
  AdvancedSettingsConfig,
  HealthStatusMap,
  InstanceRequestSettings,
  ViewerRun,
  ViewerRuntimeConfig,
} from "../../types/viewer";

type RunSidebarProps = {
  runtimeConfig?: ViewerRuntimeConfig;
  advancedSettings?: AdvancedSettingsConfig;
  health: HealthStatusMap;
  instanceSettings: Record<string, InstanceRequestSettings>;
  runs: Record<string, ViewerRun>;
  activeRunKey?: string;
  onRunSelect: (runKey: string) => void;
  onUploadFiles: (files: Array<{ fileName: string; input: unknown }>) => void;
  onUpdateSetting: (
    instanceId: string,
    patch: Partial<InstanceRequestSettings>,
  ) => void;
};

function StatusIcon({ run }: { run?: ViewerRun }) {
  if (!run) {
    return <Activity size={16} />;
  }
  if (run.status === "pending") {
    return <LoaderCircle className="spin" size={16} />;
  }
  if (run.status === "success") {
    return <CheckCircle2 size={16} />;
  }
  if (run.status === "error") {
    return <CircleAlert size={16} />;
  }
  return <Activity size={16} />;
}

function ToggleRow({
  label,
  checked,
  disabled,
  onChange,
}: {
  label: string;
  checked: boolean;
  disabled?: boolean;
  onChange: (checked: boolean) => void;
}) {
  return (
    <label className={`toggle-row${disabled ? " disabled" : ""}`}>
      <span className="toggle-row__label">{label}</span>
      <input
        className="toggle-switch__input"
        checked={checked}
        disabled={disabled}
        type="checkbox"
        onChange={(event) => onChange(event.target.checked)}
      />
      <span className="toggle-switch" aria-hidden="true">
        <span className="toggle-switch__thumb" />
      </span>
    </label>
  );
}

export function RunSidebar({
  runtimeConfig,
  advancedSettings,
  health,
  instanceSettings,
  runs,
  activeRunKey,
  onRunSelect,
  onUploadFiles,
  onUpdateSetting,
}: RunSidebarProps) {
  const uploadInputId = useId();
  const uploadedRuns = Object.values(runs).filter((run) => run.source === "file");

  return (
    <aside className="sidebar">
      <section className="panel sidebar-panel">
        <div className="panel__header">
          <div>
            <p className="eyebrow">Compare mode</p>
            <h2>Backend pathways</h2>
          </div>
        </div>
        <div className="stack">
          {runtimeConfig?.instances.map((instance) => {
            const runKey = `api:${instance.id}`;
            const run = runs[runKey];
            const settings = instanceSettings[instance.id] ?? instance.defaults;
            const modelCapabilities =
              advancedSettings?.llm_models[settings.model_type] ?? null;

            return (
              <article
                className={`instance-card${activeRunKey === runKey ? " active" : ""}`}
                key={instance.id}
              >
                <button
                  className="instance-card__header"
                  onClick={() => {
                    if (run) {
                      onRunSelect(runKey);
                    }
                  }}
                  type="button"
                >
                  <div>
                    <h3>{instance.label}</h3>
                    <p>{instance.baseUrl}</p>
                  </div>
                  <div className="status-pill">
                    <StatusIcon run={run} />
                    <span>{run?.status ?? "idle"}</span>
                  </div>
                </button>

                <div className="instance-card__meta">
                  <span
                    className={`health-pill ${health[instance.id]?.state ?? "unknown"}`}
                  >
                    {health[instance.id]?.message ?? "Health unknown"}
                  </span>
                  {run?.durationMs ? <span>{run.durationMs} ms</span> : null}
                </div>

                <div className="form-grid">
                  <label>
                    <span>LLM model</span>
                    <select
                      value={settings.model_type}
                      onChange={(event) =>
                        onUpdateSetting(instance.id, {
                          model_type: event.target.value,
                        })
                      }
                    >
                      {advancedSettings
                        ? Object.entries(advancedSettings.llm_models).map(([key, model]) => (
                            <option key={key} value={key}>
                              {model.display_name}
                            </option>
                          ))
                        : (
                          <option value={settings.model_type}>{settings.model_type}</option>
                        )}
                    </select>
                  </label>
                  <label>
                    <span>AZ model</span>
                    <select
                      value={settings.model_version}
                      onChange={(event) =>
                        onUpdateSetting(instance.id, {
                          model_version: event.target.value,
                        })
                      }
                    >
                      {advancedSettings
                        ? Object.entries(advancedSettings.az_models).map(([key, model]) => (
                            <option key={key} value={key}>
                              {model.display_name}
                            </option>
                          ))
                        : (
                          <option value={settings.model_version}>{settings.model_version}</option>
                        )}
                    </select>
                  </label>
                </div>

                <div className="toggle-group">
                  <ToggleRow
                    label="Advanced prompt"
                    checked={settings.advanced_prompt}
                    disabled={!modelCapabilities?.supports_advanced_prompt}
                    onChange={(checked) =>
                      onUpdateSetting(instance.id, { advanced_prompt: checked })
                    }
                  />
                  <ToggleRow
                    label="Stability check"
                    checked={settings.stability_flag}
                    disabled={!modelCapabilities?.supports_stability_check}
                    onChange={(checked) =>
                      onUpdateSetting(instance.id, { stability_flag: checked })
                    }
                  />
                  <ToggleRow
                    label="Hallucination check"
                    checked={settings.hallucination_check}
                    disabled={!modelCapabilities?.supports_hallucination_check}
                    onChange={(checked) =>
                      onUpdateSetting(instance.id, { hallucination_check: checked })
                    }
                  />
                  <ToggleRow
                    label="Protecting groups"
                    checked={settings.use_protecting_group_feature}
                    onChange={(checked) =>
                      onUpdateSetting(instance.id, {
                        use_protecting_group_feature: checked,
                      })
                    }
                  />
                </div>

                {run?.error ? <p className="inline-error">{run.error}</p> : null}
              </article>
            );
          })}
        </div>
      </section>

      {uploadedRuns.length ? (
        <section className="panel sidebar-panel">
          <div className="panel__header">
            <div>
              <p className="eyebrow">Local runs</p>
              <h2>Uploaded pathways</h2>
            </div>
            <label className="ghost-button upload-trigger" htmlFor={uploadInputId}>
              <Upload size={16} />
              Add files
            </label>
            <input
              id={uploadInputId}
              accept=".json,application/json"
              className="sr-only"
              type="file"
              multiple
              onChange={async (event) => {
                const files = Array.from(event.target.files ?? []);
                if (!files.length) {
                  return;
                }

                try {
                  const parsedFiles = await Promise.all(
                    files.map(async (file) => ({
                      fileName: file.name,
                      input: JSON.parse(await file.text()),
                    })),
                  );
                  onUploadFiles(parsedFiles);
                } catch (error) {
                  console.error("Failed to load uploaded files from sidebar.", error);
                } finally {
                  event.target.value = "";
                }
              }}
            />
          </div>
          <div className="stack">
            {uploadedRuns.map((run) => (
              <button
                className={`instance-card__header uploaded${activeRunKey === run.key ? " active" : ""}`}
                key={run.key}
                onClick={() => onRunSelect(run.key)}
                type="button"
              >
                <div>
                  <h3>{run.label}</h3>
                  <p>Uploaded JSON</p>
                </div>
                <div className="status-pill success">
                  <CheckCircle2 size={16} />
                  <span>loaded</span>
                </div>
              </button>
            ))}
          </div>
        </section>
      ) : null}
    </aside>
  );
}
