const keys = {
  apiKey: "deepretro.viewer.apiKey",
  lastSmiles: "deepretro.viewer.lastSmiles",
  activeRun: "deepretro.viewer.activeRun",
  theme: "deepretro.viewer.theme",
} as const;

function read(key: string) {
  if (typeof window === "undefined" || !window.localStorage) {
    return "";
  }
  try {
    return window.localStorage.getItem(key) ?? "";
  } catch {
    return "";
  }
}

function write(key: string, value: string) {
  if (typeof window === "undefined" || !window.localStorage) {
    return;
  }

  try {
    if (value) {
      window.localStorage.setItem(key, value);
    } else {
      window.localStorage.removeItem(key);
    }
  } catch {
    // Storage can be unavailable (private mode, blocked cookies); persistence
    // is best-effort and the viewer must keep working without it.
  }
}

export const storage = {
  keys,
  getApiKey: () => read(keys.apiKey),
  setApiKey: (value: string) => write(keys.apiKey, value),
  getLastSmiles: () => read(keys.lastSmiles),
  setLastSmiles: (value: string) => write(keys.lastSmiles, value),
  getActiveRun: () => read(keys.activeRun),
  setActiveRun: (value: string) => write(keys.activeRun, value),
  getTheme: () => read(keys.theme),
  setTheme: (value: string) => write(keys.theme, value),
};
