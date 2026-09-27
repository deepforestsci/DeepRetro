# DeepRetro viewer redesign — dual-theme workbench

Date: 2026-07-25
Status: approved (user chose to keep both mockup directions as light/dark themes)

## Goal

Replace `viewer_new`'s soft consumer look (pastel gradients, glassy panels, 32px radii,
3-row topbar) with a professional, data-dense workbench. Two visual directions were
mocked up (`viewer_new/mockups/`); the user liked both, so they ship as one token
system with a light/dark toggle:

- **Light — "Instrument"** (mockup A): white panels, 1px hairline borders, 7–10px
  radii, restrained blue accent `#2557d6`, dot-grid canvas.
- **Dark — "Night console"** (mockup B): graphite panels (`#0e1217`/`#151b22`),
  amber accent `#f0a72e`, benzene-hexagon canvas texture, Space Grotesk headings.

Shared: Inter for UI, JetBrains Mono for SMILES/formulas/metrics/eyebrows, a 52px
single-row toolbar, and a confidence rail on every step node (green ≥ 0.8, amber
≥ 0.5, red below).

## Architecture

- **Tokens**: all colors/fonts as CSS variables on `:root` (light) overridden by
  `html[data-theme="dark"]`. Components never hardcode colors.
- **Theme state**: `src/features/theme/useTheme.ts` — reads
  `deepretro.viewer.theme` from localStorage, falls back to
  `prefers-color-scheme`, writes `data-theme` on `<html>`, persists on toggle.
  Toggle button (Sun/Moon) lives at the right end of the topbar.
- **Molecules in dark mode**: RDKit SVGs are black-on-white; dark theme applies
  `filter: invert(0.92) hue-rotate(180deg)` so structures render as light line
  work on a dark plate, keeping heteroatom hues roughly correct.
- **Graph chrome**: edge/marker/background colors move out of inline styles
  (`layout.ts`, `GraphCanvas.tsx`) into CSS so they follow the theme. Light canvas
  uses React Flow's dot Background; dark hides it and uses a hex-grid data-URI.
- **Step nodes**: hover overlay removed (it hid the molecule). New card: mono
  step tag + confidence value header, molecule, formula + component-count footer,
  confidence rail, and hover-revealed Edit / Partial-rerun icon buttons.

## Components touched

`global.css` (rewrite), `index.html` (fonts), `App.tsx` (single-row topbar,
status strip, theme toggle), `StepNode.tsx` (new card), `GraphCanvas.tsx` +
`layout.ts` (theme-aware graph chrome), `storage.ts` (theme key),
`useTheme.ts` (+ test). All other components keep their markup; the stylesheet
reskins them via existing class names.

## Error handling

- No stored theme + no matchMedia (jsdom/SSR) → default light.
- Invalid stored value → treated as unset.

## Testing

- Unit test for `useTheme`: default, persistence, toggle, `data-theme` attribute.
- Existing tests must stay green; `tsc -b`, eslint, and `vite build` must pass.
- Visual verification: headless Chrome screenshots of both themes.
