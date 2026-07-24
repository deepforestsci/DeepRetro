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
