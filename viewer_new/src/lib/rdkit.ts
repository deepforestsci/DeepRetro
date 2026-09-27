import rdkitLoaderModule, { type RDKitLoader, type RDKitModule } from "@rdkit/rdkit";
import rdkitWasm from "@rdkit/rdkit/dist/RDKit_minimal.wasm?url";

const svgCache = new Map<string, string>();

let rdkitPromise: Promise<RDKitModule> | null = null;

const initRDKitModule = (
  (rdkitLoaderModule as unknown as { default?: RDKitLoader }).default ??
  (rdkitLoaderModule as unknown as RDKitLoader)
) as RDKitLoader;

export async function getRDKit() {
  if (!rdkitPromise) {
    rdkitPromise = initRDKitModule({
      locateFile: () => rdkitWasm,
    });
  }
  return rdkitPromise;
}

export async function renderMoleculeSvg(
  smiles: string,
  width = 220,
  height = 140,
) {
  const cacheKey = `${smiles}:${width}:${height}`;
  const cached = svgCache.get(cacheKey);
  if (cached) {
    return cached;
  }

  const RDKit = await getRDKit();
  const molecule = RDKit.get_mol(smiles);

  if (!molecule || !molecule.is_valid()) {
    molecule?.delete();
    throw new Error(`Unable to render SMILES: ${smiles}`);
  }

  try {
    if (!molecule.has_coords()) {
      molecule.set_new_coords(true);
    }
    const svg = molecule.get_svg(width, height);
    svgCache.set(cacheKey, svg);
    return svg;
  } finally {
    molecule.delete();
  }
}
