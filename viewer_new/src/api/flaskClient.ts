import type {
  AnalysisRequest,
  PathwayResult,
  ViewerRuntimeConfig,
} from "../types/viewer";

type InstanceEndpointConfig = {
  baseUrl: string;
  endpoints: ViewerRuntimeConfig["endpoints"];
  apiKey: string;
};

async function requestJson<T>(
  url: string,
  init: RequestInit,
): Promise<T> {
  const response = await fetch(url, init);
  const data = (await response.json().catch(() => ({}))) as T & {
    error?: string;
  };

  if (!response.ok) {
    throw new Error(data.error || `Request failed with ${response.status}`);
  }

  return data;
}

function withHeaders(apiKey: string) {
  return {
    "Content-Type": "application/json",
    "X-API-KEY": apiKey,
  };
}

export function createFlaskClient(config: InstanceEndpointConfig) {
  const buildUrl = (path: string) => `${config.baseUrl}${path}`;

  return {
    async retrosynthesis(body: AnalysisRequest) {
      return requestJson<PathwayResult>(buildUrl(config.endpoints.retrosynthesis), {
        method: "POST",
        headers: withHeaders(config.apiKey),
        body: JSON.stringify(body),
      });
    },
    async rerun(body: AnalysisRequest) {
      return requestJson<PathwayResult>(buildUrl(config.endpoints.rerun), {
        method: "POST",
        headers: withHeaders(config.apiKey),
        body: JSON.stringify(body),
      });
    },
    async partialRerun(body: AnalysisRequest & { steps: string }) {
      return requestJson<PathwayResult>(buildUrl(config.endpoints.partialRerun), {
        method: "POST",
        headers: withHeaders(config.apiKey),
        body: JSON.stringify(body),
      });
    },
    async saveEditedResult(smiles: string, result: PathwayResult) {
      return requestJson<{ message: string }>(buildUrl(config.endpoints.saveEdited), {
        method: "POST",
        headers: withHeaders(config.apiKey),
        body: JSON.stringify({ smiles, result }),
      });
    },
    async health() {
      return requestJson<{ status: string }>(buildUrl(config.endpoints.health), {
        method: "GET",
        headers: {
          "X-API-KEY": config.apiKey,
        },
      });
    },
  };
}
