import { cp, mkdir } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const scriptDir = dirname(fileURLToPath(import.meta.url));
const source = resolve(scriptDir, "../../config/advanced_settings.json");
const target = resolve(scriptDir, "../public/advanced-settings.json");

await mkdir(dirname(target), { recursive: true });
await cp(source, target);

console.log(`Synced advanced settings to ${target}`);
