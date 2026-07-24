# Viewer redesign implementation — detailed instructions

Companion to the approved spec:
`docs/superpowers/specs/2026-07-25-viewer-theme-redesign-design.md`.
Reference mockups (open in any browser): `viewer_new/mockups/mockup-a-instrument.html`
(light) and `viewer_new/mockups/mockup-b-console.html` (dark). Reference screenshots:
`/tmp/mockup-a2.png`, `/tmp/mockup-b2.png`.

Goal: reskin `viewer_new` so mockup A is the light theme and mockup B is the dark
theme, switched by a Sun/Moon toggle in the topbar, persisted to localStorage.

## Status checklist

- [x] 1. Theme storage key — `src/config/storage.ts`: added `theme` key +
      `getTheme`/`setTheme` (localStorage key `deepretro.viewer.theme`).
- [x] 2. Theme hook — `src/features/theme/useTheme.ts`: returns
      `{ theme, toggleTheme }`; initial value = stored theme, else
      `prefers-color-scheme`, else light; effect sets
      `document.documentElement.dataset.theme` and persists. Unit tests in
      `useTheme.test.ts` (default / restore / invalid value / toggle+persist).
- [x] 3. Fonts — `viewer_new/index.html`: Google Fonts links for Inter,
      Space Grotesk, JetBrains Mono.
- [x] 4. Stylesheet — `src/styles/global.css`: full rewrite. `:root` holds the
      light "Instrument" tokens; `html[data-theme="dark"]` overrides with the
      "Night console" tokens (amber accent, graphite panels, hex-grid canvas via
      data-URI on `.canvas-flow`, RDKit SVG inversion via
      `filter: invert(0.93) hue-rotate(180deg)` on `.molecule-svg svg`).
      All existing component class names were kept.
- [x] 5. `App.tsx` — single-row 52px topbar:
      `brand (hexagon SVG mark + "DeepRetro" / "RETROSYNTHESIS") · segmented
      Analyze/Upload · SMILES field (mono) · API-key field · "Run analysis"
      primary button · icon-only ghost buttons (Run again / Download JSON /
      Save edits) · spacer · health chip (colored dot + text) · theme toggle
      (Moon in light, Sun in dark)`. Feedback/error pills render in a
      `.status-strip` row below the topbar only when non-empty.
      `summarizeHealthStates` returns `{ label, tone: ok|warn|err|muted }` and
      the chip uses `health-chip <tone>`. Uses `useTheme()`.
- [x] 6. `StepNode.tsx` — replaced the hover overlay (which hid the molecule)
      with: header row (`step-node__tag` "Target"/"Step N" + `step-node__conf`
      confidence, tinted ok ≥0.8 / warn ≥0.5 / err below), molecule preview,
      footer (`step-node__formula` + `step-node__count` "N reactants"),
      `step-node__rail` (3px confidence bar, width = confidence·100%,
      accent-colored for the target node), and `step-node__quick` — hover/
      selected-revealed icon buttons for Edit + Partial rerun (stopPropagation,
      aria-labels).
- [x] 7. Graph chrome theming — `layout.ts`: removed hardcoded white
      stroke/glow from `toFlowEdges` and the marker `color`; `GraphCanvas.tsx`:
      removed `defaultEdgeOptions` inline stroke, `MiniMap` `nodeColor`, and
      `Background` `color`/`gap` props. Colors now come from CSS
      (`.react-flow__edge-path`, `.react-flow__marker path`,
      `.react-flow__background circle`, minimap/controls rules in global.css).
- [x] 8. Verification — `cd viewer_new && bun install && bun run typecheck &&
      bun run test && bun run build` all pass, plus visual check (step 9 below).

Verified visually: `/tmp/viewer-light.png`, `/tmp/viewer-dark.png` (1440x900,
empty state — no backend was running). Not yet visually verified: graph nodes,
inspector, and editor drawer with a real pathway loaded; upload a pathway JSON
or run an analysis and check both themes, especially molecule legibility in
dark mode (inverted RDKit SVGs) and the step-node confidence rails.

Note: the vitest jsdom sandbox exposes no `window.localStorage`; theme tests
stub it in-memory (see `useTheme.test.ts`), and `storage.ts` now degrades
gracefully when storage is missing or blocked.

## If picking this up fresh

1. Read the spec doc above, then diff `git log --oneline` for what landed.
2. Any unchecked item: the exact target markup/CSS is in the mockup HTML files —
   copy class structure and tokens from there; global.css already contains all
   needed classes, so component work is only JSX restructuring.
3. Do not rename existing CSS classes used by RunSidebar / InspectorPanel /
   EditorDrawer / UploadPanel — the stylesheet reskins them as-is.

## Visual verification recipe (step 9)

```bash
cd viewer_new && bun run dev --port 5199 &
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
# light
"$CHROME" --headless=new --disable-gpu --hide-scrollbars \
  --window-size=1440,900 --screenshot=/tmp/viewer-light.png http://localhost:5199
# dark (theme persists via localStorage; simplest is toggling in the UI, or:)
"$CHROME" --headless=new --disable-gpu --hide-scrollbars \
  --window-size=1440,900 --virtual-time-budget=8000 \
  --screenshot=/tmp/viewer-dark.png \
  "http://localhost:5199/#dark"   # then click the moon icon manually if needed
```

Note: the theme init reads localStorage before render, so for a true dark
screenshot either click the toggle in a real browser or pre-seed
`localStorage["deepretro.viewer.theme"]="dark"` via devtools/Playwright.
Check: topbar is one row; molecule drawings legible in both themes (dark mode
inverts RDKit SVGs); step nodes show tag/confidence/rail; edges + arrowheads
follow theme; no hover overlay hiding molecules.

## Manual theme override

`localStorage.setItem("deepretro.viewer.theme", "dark" | "light")` then reload;
delete the key to fall back to the OS preference.
