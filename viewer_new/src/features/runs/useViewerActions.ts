import { useEffect } from "react";
import { useMutation } from "@tanstack/react-query";

import { loadAdvancedSettings, loadRuntimeConfig } from "../../config/loaders";
import { createFlaskClient } from "../../api/flaskClient";
import { storage } from "../../config/storage";
import {
  buildPathwayGraph,
  canPartialRerunStep,
  getStep,
  inferTopLevelSmiles,
  updateStepInResult,
  validatePathwayResult,
  validatePathwayStep,
} from "../../lib/pathway";
import { useViewerStore } from "./store";
import type {
  AnalysisRequest,
  HealthStatusMap,
  ViewerRun,
} from "../../types/viewer";

function createRun(
  key: string,
  label: string,
  status: ViewerRun["status"],
  source: ViewerRun["source"],
  payload: Partial<ViewerRun>,
): ViewerRun {
  return {
    key,
    label,
    source,
    status,
    dirty: false,
    lastUpdatedAt: Date.now(),
    ...payload,
  };
}

export function useBootstrapViewer() {
  const bootstrap = useViewerStore((state) => state.bootstrap);

  useEffect(() => {
    let cancelled = false;

    Promise.all([loadRuntimeConfig(), loadAdvancedSettings()])
      .then(([runtimeConfig, advancedSettings]) => {
        if (!cancelled) {
          bootstrap({ runtimeConfig, advancedSettings });
        }
      })
      .catch((error) => {
        console.error("Failed to bootstrap viewer:", error);
      });

    return () => {
      cancelled = true;
    };
  }, [bootstrap]);
}

export function useViewerActions() {
  const runtimeConfig = useViewerStore((state) => state.runtimeConfig);
  const instanceSettings = useViewerStore((state) => state.instanceSettings);
  const apiKey = useViewerStore((state) => state.apiKey);
  const runs = useViewerStore((state) => state.runs);
  const activeRunKey = useViewerStore((state) => state.activeRunKey);
  const setRuns = useViewerStore((state) => state.setRuns);
  const setHealth = useViewerStore((state) => state.setHealth);
  const setActiveRun = useViewerStore((state) => state.setActiveRun);
  const setSelectedStep = useViewerStore((state) => state.setSelectedStep);

  const healthMutation = useMutation({
    mutationFn: async () => {
      if (!runtimeConfig) {
        return {};
      }

      const statusMap: HealthStatusMap = {};

      if (!apiKey) {
        for (const instance of runtimeConfig.instances) {
          statusMap[instance.id] = {
            state: "unauthorized",
            message: "API key required",
          };
        }
        return statusMap;
      }

      await Promise.all(
        runtimeConfig.instances.map(async (instance) => {
          const client = createFlaskClient({
            baseUrl: instance.baseUrl,
            endpoints: runtimeConfig.endpoints,
            apiKey,
          });

          try {
            await client.health();
            statusMap[instance.id] = {
              state: "healthy",
              message: "Healthy",
            };
          } catch (error) {
            const message = error instanceof Error ? error.message : "Health check failed";
            statusMap[instance.id] = {
              state: message.toLowerCase().includes("unauthorized")
                ? "unauthorized"
                : "error",
              message,
            };
          }
        }),
      );

      return statusMap;
    },
    onSuccess: (health) => {
      setHealth(health);
    },
  });

  const analyzeMutation = useMutation({
    mutationFn: async ({
      smiles,
      rerun,
    }: {
      smiles: string;
      rerun: boolean;
    }) => {
      if (!runtimeConfig) {
        throw new Error("Runtime config has not loaded yet.");
      }

      if (!apiKey) {
        throw new Error("An API key is required before running the viewer.");
      }

      const pendingRuns = runtimeConfig.instances.reduce<Record<string, ViewerRun>>(
        (accumulator, instance) => {
          const key = `api:${instance.id}`;
          accumulator[key] = createRun(key, instance.label, "pending", "api", {
            instanceId: instance.id,
            request: {
              smiles,
              ...instanceSettings[instance.id],
            },
          });
          return accumulator;
        },
        {},
      );

      setRuns((previous) => ({
        ...previous,
        ...pendingRuns,
      }));

      const results = await Promise.all(
        runtimeConfig.instances.map(async (instance) => {
          const request: AnalysisRequest = {
            smiles,
            ...instanceSettings[instance.id],
          };

          const startedAt = performance.now();
          const client = createFlaskClient({
            baseUrl: instance.baseUrl,
            endpoints: runtimeConfig.endpoints,
            apiKey,
          });

          try {
            const rawResult = rerun
              ? await client.rerun(request)
              : await client.retrosynthesis(request);
            return {
              key: `api:${instance.id}`,
              run: createRun(`api:${instance.id}`, instance.label, "success", "api", {
                instanceId: instance.id,
                request,
                rawResult,
                graph: buildPathwayGraph(rawResult),
                durationMs: Math.round(performance.now() - startedAt),
              }),
            };
          } catch (error) {
            return {
              key: `api:${instance.id}`,
              run: createRun(`api:${instance.id}`, instance.label, "error", "api", {
                instanceId: instance.id,
                request,
                error: error instanceof Error ? error.message : "Request failed",
                durationMs: Math.round(performance.now() - startedAt),
              }),
            };
          }
        }),
      );

      return { smiles, results };
    },
    onSuccess: ({ smiles, results }) => {
      storage.setLastSmiles(smiles);
      setRuns((previous) => {
        const nextRuns = { ...previous };
        for (const { key, run } of results) {
          nextRuns[key] = run;
        }
        return nextRuns;
      });

      const firstSuccessfulRun = results.find((item) => item.run.status === "success");
      setActiveRun(firstSuccessfulRun?.key ?? results[0]?.key);
      setSelectedStep("0");
      void healthMutation.mutateAsync();
    },
  });

  const saveMutation = useMutation({
    mutationFn: async (runKey: string) => {
      if (!runtimeConfig) {
        throw new Error("Runtime config has not loaded yet.");
      }

      const run = runs[runKey];
      if (!run?.instanceId || !run.rawResult || !run.request) {
        throw new Error("Only API runs with loaded results can be saved.");
      }

      const instance = runtimeConfig.instances.find(
        (candidate) => candidate.id === run.instanceId,
      );
      if (!instance) {
        throw new Error("Missing backend instance configuration.");
      }

      const client = createFlaskClient({
        baseUrl: instance.baseUrl,
        endpoints: runtimeConfig.endpoints,
        apiKey,
      });

      const smiles = run.request.smiles || inferTopLevelSmiles(run.rawResult);
      await client.saveEditedResult(smiles, run.rawResult);
      return runKey;
    },
    onSuccess: (runKey) => {
      setRuns((previous) => ({
        ...previous,
        [runKey]: {
          ...previous[runKey],
          dirty: false,
          lastUpdatedAt: Date.now(),
        },
      }));
    },
  });

  const partialRerunMutation = useMutation({
    mutationFn: async ({
      runKey,
      stepId,
    }: {
      runKey: string;
      stepId: string;
    }) => {
      if (!runtimeConfig) {
        throw new Error("Runtime config has not loaded yet.");
      }

      const run = runs[runKey];
      if (!run?.instanceId || !run.request || !run.rawResult) {
        throw new Error("Select an API run before triggering partial rerun.");
      }

      const step = getStep(run.rawResult, stepId);
      const eligibility = canPartialRerunStep(step);
      if (!eligibility.enabled) {
        throw new Error(eligibility.reason);
      }

      const instance = runtimeConfig.instances.find(
        (candidate) => candidate.id === run.instanceId,
      );
      if (!instance) {
        throw new Error("Missing backend instance configuration.");
      }

      const client = createFlaskClient({
        baseUrl: instance.baseUrl,
        endpoints: runtimeConfig.endpoints,
        apiKey,
      });

      const request = {
        smiles: run.request.smiles,
        steps: stepId,
        ...instanceSettings[run.instanceId],
      };

      const rawResult = await client.partialRerun(request);
      return { runKey, rawResult };
    },
    onSuccess: ({ runKey, rawResult }) => {
      const currentRun = runs[runKey];
      setRuns((previous) => ({
        ...previous,
        [runKey]: {
          ...previous[runKey],
          rawResult,
          graph: buildPathwayGraph(rawResult),
          dirty: false,
          error: undefined,
          status: "success",
          lastUpdatedAt: Date.now(),
          request: currentRun?.request
            ? {
                ...currentRun.request,
                smiles: currentRun.request.smiles,
              }
            : currentRun?.request,
        },
      }));
      setActiveRun(runKey);
      setSelectedStep("0");
    },
  });

  useEffect(() => {
    if (!runtimeConfig) {
      return;
    }
    void healthMutation.mutateAsync();
  }, [apiKey, runtimeConfig, healthMutation]);

  return {
    runtimeConfig,
    activeRun: activeRunKey ? runs[activeRunKey] : undefined,
    activeRunKey,
    isAnalyzing: analyzeMutation.isPending,
    isSaving: saveMutation.isPending,
    isPartialRerunning: partialRerunMutation.isPending,
    refreshHealth: () => healthMutation.mutateAsync(),
    analyzeAll: (smiles: string, rerun = false) =>
      analyzeMutation.mutateAsync({ smiles, rerun }),
    loadUploadedResult: (fileName: string, input: unknown) => {
      const validation = validatePathwayResult(input);
      if (!validation.data) {
        throw new Error(validation.issues.map((issue) => issue.message).join("\n"));
      }

      const key = `file:${fileName}`;
      const run = createRun(key, fileName, "success", "file", {
        rawResult: validation.data,
        graph: buildPathwayGraph(validation.data),
      });

      setRuns((previous) => ({
        ...previous,
        [key]: run,
      }));
      setActiveRun(key);
      setSelectedStep("0");
      return key;
    },
    applyStructuredStepEdit: (runKey: string, stepId: string, nextStep: unknown) => {
      const validation = validatePathwayStep(nextStep);
      if (!validation.data) {
        throw new Error(validation.issues.map((issue) => issue.message).join("\n"));
      }

      const run = runs[runKey];
      if (!run?.rawResult) {
        throw new Error("No run data is loaded.");
      }

      const rawResult = updateStepInResult(run.rawResult, stepId, validation.data);
      const graph = buildPathwayGraph(rawResult);

      setRuns((previous) => ({
        ...previous,
        [runKey]: {
          ...previous[runKey],
          rawResult,
          graph,
          dirty: true,
          lastUpdatedAt: Date.now(),
        },
      }));
    },
    applyRawJsonEdit: (runKey: string, input: unknown) => {
      const validation = validatePathwayResult(input);
      if (!validation.data) {
        throw new Error(validation.issues.map((issue) => issue.message).join("\n"));
      }

      setRuns((previous) => ({
        ...previous,
        [runKey]: {
          ...previous[runKey],
          rawResult: validation.data,
          graph: buildPathwayGraph(validation.data),
          dirty: true,
          lastUpdatedAt: Date.now(),
        },
      }));
      setSelectedStep("0");
    },
    saveRunEdits: (runKey: string) => saveMutation.mutateAsync(runKey),
    partialRerun: (runKey: string, stepId: string) =>
      partialRerunMutation.mutateAsync({ runKey, stepId }),
  };
}
