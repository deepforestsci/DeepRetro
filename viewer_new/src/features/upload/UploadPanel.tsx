import { useId, useState } from "react";
import { FileUp } from "lucide-react";

type UploadPanelProps = {
  onLoadFiles: (files: Array<{ fileName: string; input: unknown }>) => void;
};

export function UploadPanel({ onLoadFiles }: UploadPanelProps) {
  const inputId = useId();
  const [error, setError] = useState("");

  return (
    <section className="panel upload-panel">
      <div className="panel__header">
        <div>
          <p className="eyebrow">Upload mode</p>
          <h2>Inspect existing pathway JSON</h2>
        </div>
      </div>
      <label className="upload-dropzone" htmlFor={inputId}>
        <FileUp size={28} />
        <strong>Choose one or more `.json` pathway exports</strong>
        <span>
          Uploaded pathways stay local in the viewer until you explicitly save edited API runs.
        </span>
      </label>
      <input
        id={inputId}
        accept=".json,application/json"
        className="sr-only"
        type="file"
        multiple
        onChange={async (event) => {
          const files = Array.from(event.target.files ?? []);
          if (!files.length) {
            return;
          }

          try {
            setError("");
            const parsedFiles = await Promise.all(
              files.map(async (file) => ({
                fileName: file.name,
                input: JSON.parse(await file.text()),
              })),
            );
            onLoadFiles(parsedFiles);
          } catch (nextError) {
            setError(
              nextError instanceof Error ? nextError.message : "Failed to load the file.",
            );
          } finally {
            event.target.value = "";
          }
        }}
      />
      {error ? <p className="inline-error">{error}</p> : null}
    </section>
  );
}
