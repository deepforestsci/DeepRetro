import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, test, vi } from "vitest";

import { UploadPanel } from "./UploadPanel";

describe("UploadPanel", () => {
  test("parses a valid JSON upload and forwards it to the callback", async () => {
    const onLoad = vi.fn();
    const user = userEvent.setup();

    render(<UploadPanel onLoad={onLoad} />);

    const input = screen.getByLabelText(/choose a `.json` pathway export/i);
    const file = new File([JSON.stringify({ steps: [] })], "pathway.json", {
      type: "application/json",
    });

    await user.upload(input, file);

    expect(onLoad).toHaveBeenCalledWith("pathway.json", { steps: [] });
  });

  test("surfaces a parsing error for malformed JSON", async () => {
    const onLoad = vi.fn();
    const user = userEvent.setup();

    render(<UploadPanel onLoad={onLoad} />);

    const input = screen.getByLabelText(/choose a `.json` pathway export/i);
    const file = new File(["not-json"], "broken.json", {
      type: "application/json",
    });

    await user.upload(input, file);

    expect(onLoad).not.toHaveBeenCalled();
    expect(screen.getByText(/unexpected token/i)).toBeInTheDocument();
  });
});
