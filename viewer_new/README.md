# `viewer_new`

Standalone Bun + Vite + React viewer for DeepRetro.

## Commands

```bash
bun install
bun run dev
bun run build
bun run typecheck
bun run test
```

## Notes

- The legacy `viewer/` app is intentionally untouched.
- `bun run sync:config` copies the root `config/advanced_settings.json` into `public/advanced-settings.json`.
- Runtime backend configuration lives in `public/runtime-config.json`.
- Multiple backends are displayed in a vertically scrollable comparison panel. Use the panel button beside "Backend pathways" to collapse it to a compact 72px run rail and reclaim graph space. Uploaded pathways and the add-files action remain available in the compact rail.
