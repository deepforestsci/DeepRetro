const keys = {
  apiKey: "deepretro.viewer.apiKey",
  lastSmiles: "deepretro.viewer.lastSmiles",
  activeRun: "deepretro.viewer.activeRun",
} as const;

function read(key: string) {
  if (typeof window === "undefined") {
    return "";
  }
  return window.localStorage.getItem(key) ?? "";
}

function write(key: string, value: string) {
  if (typeof window === "undefined") {
    return;
  }

  if (value) {
    window.localStorage.setItem(key, value);
  } else {
    window.localStorage.removeItem(key);
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
};
