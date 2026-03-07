import { useId, useState } from "react";
import { FileUp } from "lucide-react";

type UploadPanelProps = {
  onLoad: (fileName: string, input: unknown) => void;
};

export function UploadPanel({ onLoad }: UploadPanelProps) {
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
        <strong>Choose a `.json` pathway export</strong>
        <span>
          Uploaded pathways stay local in the viewer until you explicitly save edited API runs.
        </span>
      </label>
      <input
        id={inputId}
        accept=".json,application/json"
        className="sr-only"
        type="file"
        onChange={async (event) => {
          const file = event.target.files?.[0];
          if (!file) {
            return;
          }

          try {
            setError("");
            const text = await file.text();
            onLoad(file.name, JSON.parse(text));
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
