import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, test, vi } from "vitest";

import { RunSidebar } from "./RunSidebar";
import type {
  AdvancedSettingsConfig,
  ViewerRun,
  ViewerRuntimeConfig,
} from "../../types/viewer";

const runtimeConfig: ViewerRuntimeConfig = {
  instances: [
    {
      id: "backend-a",
      label: "Backend A · Claude 4 Opus",
      baseUrl: "http://backend-a.example:5000",
      defaults: {
        model_type: "claude4opus",
        model_version: "Pistachio_100+",
        advanced_prompt: true,
        stability_flag: true,
        hallucination_check: true,
        use_protecting_group_feature: false,
      },
    },
  ],
  endpoints: {
    retrosynthesis: "/api/retrosynthesis",
    rerun: "/api/rerun_retrosynthesis",
    partialRerun: "/api/partial_rerun",
    saveEdited: "/api/save_edited_result",
    health: "/api/health",
  },
};

const advancedSettings: AdvancedSettingsConfig = {
  llm_models: {
    claude4opus: {
      internal_name: "claude-opus",
      display_name: "Claude 4 Opus",
      supports_advanced_prompt: true,
      supports_stability_check: true,
      supports_hallucination_check: true,
    },
  },
  az_models: {
    "Pistachio_100+": { display_name: "Pistachio (100+)" },
  },
  defaults: runtimeConfig.instances[0].defaults,
};

const uploadedRun: ViewerRun = {
  key: "file:aspirin-pathway.json",
  label: "aspirin-pathway.json",
  source: "file",
  status: "success",
  dirty: false,
  lastUpdatedAt: 1,
};

describe("RunSidebar", () => {
  test("renders a compact backend rail and can request expansion", async () => {
    const onCollapsedChange = vi.fn();
    const user = userEvent.setup();

    render(
      <RunSidebar
        activeRunKey={undefined}
        advancedSettings={advancedSettings}
        collapsed
        health={{}}
        instanceSettings={{ "backend-a": runtimeConfig.instances[0].defaults }}
        onCollapsedChange={onCollapsedChange}
        onRunSelect={vi.fn()}
        onUpdateSetting={vi.fn()}
        onUploadFiles={vi.fn()}
        runs={{ [uploadedRun.key]: uploadedRun }}
        runtimeConfig={runtimeConfig}
      />,
    );

    expect(screen.getByRole("button", { name: "Select Backend A · Claude 4 Opus" })).toBeVisible();
    expect(screen.getByRole("button", { name: "Select uploaded aspirin-pathway.json" })).toBeVisible();
    expect(screen.getByLabelText("Add pathway files")).toBeVisible();
    expect(screen.queryByText("LLM model")).not.toBeInTheDocument();

    await user.click(screen.getByRole("button", { name: "Expand backend panel" }));
    expect(onCollapsedChange).toHaveBeenCalledWith(false);
  });
});
