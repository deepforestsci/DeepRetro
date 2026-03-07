import { create } from "zustand";

import { storage } from "../../config/storage";
import type {
  AdvancedSettingsConfig,
  HealthStatusMap,
  InstanceRequestSettings,
  ViewerMode,
  ViewerRun,
  ViewerRuntimeConfig,
  ViewerStoreState,
} from "../../types/viewer";

type ViewerActions = {
  bootstrap: (payload: {
    runtimeConfig: ViewerRuntimeConfig;
    advancedSettings: AdvancedSettingsConfig;
  }) => void;
  setMode: (mode: ViewerMode) => void;
  setApiKey: (apiKey: string) => void;
  setCurrentSmiles: (smiles: string) => void;
  setHealth: (health: HealthStatusMap) => void;
  setActiveRun: (runKey?: string) => void;
  setSelectedStep: (stepId?: string) => void;
  setRuns: (
    updater: (runs: Record<string, ViewerRun>) => Record<string, ViewerRun>,
  ) => void;
  updateInstanceSetting: (
    instanceId: string,
    patch: Partial<InstanceRequestSettings>,
  ) => void;
  setEditorOpen: (open: boolean) => void;
};

type ViewerStore = ViewerStoreState & ViewerActions;

export const useViewerStore = create<ViewerStore>((set, get) => ({
  mode: "analyze",
  apiKey: storage.getApiKey(),
  currentSmiles: storage.getLastSmiles(),
  activeRunKey: storage.getActiveRun() || undefined,
  selectedStepId: undefined,
  runs: {},
  instanceSettings: {},
  health: {},
  editorOpen: false,
  bootstrap: ({ runtimeConfig, advancedSettings }) => {
    const existing = get().instanceSettings;
    const nextSettings = Object.fromEntries(
      runtimeConfig.instances.map((instance) => [
        instance.id,
        existing[instance.id] ?? instance.defaults,
      ]),
    );

    set({
      runtimeConfig,
      advancedSettings,
      instanceSettings: nextSettings,
    });
  },
  setMode: (mode) => set({ mode }),
  setApiKey: (apiKey) => {
    storage.setApiKey(apiKey);
    set({ apiKey });
  },
  setCurrentSmiles: (smiles) => {
    storage.setLastSmiles(smiles);
    set({ currentSmiles: smiles });
  },
  setHealth: (health) => set({ health }),
  setActiveRun: (runKey) => {
    storage.setActiveRun(runKey ?? "");
    const nextRun = runKey ? get().runs[runKey] : undefined;
    const nextStepId =
      nextRun?.graph?.nodes.find((node) => !node.isVirtualRoot)?.stepId ??
      nextRun?.graph?.virtualRoot.stepId;
    set({
      activeRunKey: runKey,
      selectedStepId: nextStepId,
    });
  },
  setSelectedStep: (stepId) => set({ selectedStepId: stepId }),
  setRuns: (updater) => set((state) => ({ runs: updater(state.runs) })),
  updateInstanceSetting: (instanceId, patch) =>
    set((state) => ({
      instanceSettings: {
        ...state.instanceSettings,
        [instanceId]: {
          ...state.instanceSettings[instanceId],
          ...patch,
        },
      },
    })),
  setEditorOpen: (editorOpen) => set({ editorOpen }),
}));
