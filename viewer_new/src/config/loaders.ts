import type { AdvancedSettingsConfig, ViewerRuntimeConfig } from "../types/viewer";
import { advancedSettingsSchema, runtimeConfigSchema } from "../lib/schemas";

async function fetchJson<T>(path: string): Promise<T> {
  const response = await fetch(path);
  if (!response.ok) {
    throw new Error(`Failed to load ${path}: ${response.status}`);
  }
  return response.json() as Promise<T>;
}

export async function loadRuntimeConfig() {
  const data = await fetchJson<ViewerRuntimeConfig>("/runtime-config.json");
  return runtimeConfigSchema.parse(data);
}

export async function loadAdvancedSettings() {
  const data = await fetchJson<AdvancedSettingsConfig>("/advanced-settings.json");
  return advancedSettingsSchema.parse(data);
}
