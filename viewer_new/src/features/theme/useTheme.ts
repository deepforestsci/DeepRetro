import { useCallback, useEffect, useState } from "react";

import { storage } from "../../config/storage";

export type ViewerTheme = "light" | "dark";

function resolveInitialTheme(): ViewerTheme {
  if (typeof window !== "undefined") {
    const fromUrl = new URLSearchParams(window.location.search).get("theme");
    if (fromUrl === "light" || fromUrl === "dark") {
      return fromUrl;
    }
  }
  const stored = storage.getTheme();
  if (stored === "light" || stored === "dark") {
    return stored;
  }
  if (
    typeof window !== "undefined" &&
    window.matchMedia?.("(prefers-color-scheme: dark)").matches
  ) {
    return "dark";
  }
  return "light";
}

export function useTheme() {
  const [theme, setTheme] = useState<ViewerTheme>(resolveInitialTheme);

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
    storage.setTheme(theme);
  }, [theme]);

  const toggleTheme = useCallback(() => {
    setTheme((current) => (current === "light" ? "dark" : "light"));
  }, []);

  return { theme, toggleTheme };
}
