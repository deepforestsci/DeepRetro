import { useEffect, useState } from "react";

import { renderMoleculeSvg } from "../../lib/rdkit";

type MoleculePreviewProps = {
  smiles: string;
  label: string;
  width?: number;
  height?: number;
  compact?: boolean;
};

export function MoleculePreview({
  smiles,
  label,
  width = 180,
  height = 120,
  compact = false,
}: MoleculePreviewProps) {
  const [svg, setSvg] = useState("");
  const [error, setError] = useState("");

  useEffect(() => {
    let cancelled = false;
    setSvg("");
    setError("");

    if (!smiles) {
      setError("Missing SMILES");
      return () => {
        cancelled = true;
      };
    }

    renderMoleculeSvg(smiles, width, height)
      .then((nextSvg) => {
        if (!cancelled) {
          setSvg(nextSvg);
        }
      })
      .catch((nextError) => {
        if (!cancelled) {
          setError(nextError instanceof Error ? nextError.message : "Render failed");
        }
      });

    return () => {
      cancelled = true;
    };
  }, [height, smiles, width]);

  if (error) {
    return (
      <div className={`molecule-fallback${compact ? " compact" : ""}`}>
        <strong>{label}</strong>
        <span>{smiles}</span>
        <small>{error}</small>
      </div>
    );
  }

  if (!svg) {
    return (
      <div className={`molecule-loading${compact ? " compact" : ""}`}>
        Rendering...
      </div>
    );
  }

  return (
    <div
      className={`molecule-svg${compact ? " compact" : ""}`}
      aria-label={`Molecule ${label}`}
      dangerouslySetInnerHTML={{ __html: svg }}
    />
  );
}
