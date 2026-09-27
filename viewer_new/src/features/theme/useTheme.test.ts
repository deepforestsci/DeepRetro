import { act, renderHook } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it } from "vitest";

import { storage } from "../../config/storage";
import { useTheme } from "./useTheme";

beforeAll(() => {
  // The vitest jsdom sandbox does not expose localStorage; install an
  // in-memory stand-in so persistence behaviour can be asserted.
  const store = new Map<string, string>();
  Object.defineProperty(window, "localStorage", {
    configurable: true,
    value: {
      getItem: (key: string) => store.get(key) ?? null,
      setItem: (key: string, value: string) => {
        store.set(key, String(value));
      },
      removeItem: (key: string) => {
        store.delete(key);
      },
      clear: () => {
        store.clear();
      },
    },
  });
});

afterEach(() => {
  window.localStorage.clear();
  delete document.documentElement.dataset.theme;
  window.history.replaceState(null, "", window.location.pathname);
});

describe("useTheme", () => {
  it("defaults to light when nothing is stored", () => {
    const { result } = renderHook(() => useTheme());

    expect(result.current.theme).toBe("light");
    expect(document.documentElement.dataset.theme).toBe("light");
  });

  it("restores a stored theme", () => {
    storage.setTheme("dark");

    const { result } = renderHook(() => useTheme());

    expect(result.current.theme).toBe("dark");
    expect(document.documentElement.dataset.theme).toBe("dark");
  });

  it("ignores invalid stored values", () => {
    storage.setTheme("solarized");

    const { result } = renderHook(() => useTheme());

    expect(result.current.theme).toBe("light");
  });

  it("prefers a theme passed via URL query", () => {
    storage.setTheme("light");
    window.history.replaceState(null, "", "?theme=dark");

    const { result } = renderHook(() => useTheme());

    expect(result.current.theme).toBe("dark");
  });

  it("toggles and persists the theme", () => {
    const { result } = renderHook(() => useTheme());

    act(() => {
      result.current.toggleTheme();
    });

    expect(result.current.theme).toBe("dark");
    expect(document.documentElement.dataset.theme).toBe("dark");
    expect(storage.getTheme()).toBe("dark");

    act(() => {
      result.current.toggleTheme();
    });

    expect(result.current.theme).toBe("light");
    expect(storage.getTheme()).toBe("light");
  });
});
