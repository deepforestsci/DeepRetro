import ELK from "elkjs/lib/elk.bundled.js";
import type { Edge, Node } from "@xyflow/react";
import { MarkerType, Position } from "@xyflow/react";

import type {
  NormalizedPathwayGraph,
  NormalizedStepNode,
} from "../types/viewer";

export type StepFlowNode = Node<{
  stepNode: NormalizedStepNode;
  highlighted: boolean;
  onInspect: (stepId: string) => void;
  onEdit: (stepId: string) => void;
  onPartialRerun: (stepId: string) => void;
}>;

const elk = new ELK();

function getDimensions(stepNode: NormalizedStepNode) {
  return {
    width: stepNode.isVirtualRoot ? 250 : 230,
    height: stepNode.isVirtualRoot ? 250 : 230,
  };
}

function toFlowEdges(graph: NormalizedPathwayGraph) {
  return graph.edges.map((edge) => ({
    ...edge,
    type: "smoothstep" as const,
    animated: false,
    markerEnd: {
      type: MarkerType.ArrowClosed,
      width: 24,
      height: 24,
      color: "#9bd8ff",
    },
    style: {
      stroke: "#9bd8ff",
      strokeOpacity: 0.85,
      strokeWidth: 3,
    },
  }));
}

function buildFallbackFlowGraph(
  graph: NormalizedPathwayGraph,
  options: {
    onInspect: (stepId: string) => void;
    onEdit: (stepId: string) => void;
    onPartialRerun: (stepId: string) => void;
    matchesSearch: (stepNode: NormalizedStepNode) => boolean;
  },
) {
  const depthMap = new Map<string, number>([["step-0", 0]]);
  const queue = ["step-0"];

  while (queue.length) {
    const currentId = queue.shift()!;
    const currentDepth = depthMap.get(currentId) ?? 0;
    for (const edge of graph.edges) {
      if (edge.source !== currentId) {
        continue;
      }
      if (!depthMap.has(edge.target)) {
        depthMap.set(edge.target, currentDepth + 1);
        queue.push(edge.target);
      }
    }
  }

  const rowsByDepth = new Map<number, number>();

  return {
    nodes: graph.nodes.map((stepNode) => {
      const depth = depthMap.get(stepNode.id) ?? 0;
      const row = rowsByDepth.get(depth) ?? 0;
      rowsByDepth.set(depth, row + 1);

      return {
        id: stepNode.id,
        type: "step",
        position: {
          x: depth * 420,
          y: row * 280,
        },
        sourcePosition: Position.Right,
        targetPosition: Position.Left,
        data: {
          stepNode,
          highlighted: options.matchesSearch(stepNode),
          onInspect: options.onInspect,
          onEdit: options.onEdit,
          onPartialRerun: options.onPartialRerun,
        },
        draggable: false,
      };
    }),
    edges: toFlowEdges(graph),
  };
}

export async function buildFlowGraph(
  graph: NormalizedPathwayGraph,
  options: {
    onInspect: (stepId: string) => void;
    onEdit: (stepId: string) => void;
    onPartialRerun: (stepId: string) => void;
    matchesSearch: (stepNode: NormalizedStepNode) => boolean;
  },
): Promise<{
  nodes: StepFlowNode[];
  edges: Edge[];
}> {
  try {
    const elkGraph = await elk.layout({
      id: "deepretro-pathway",
      layoutOptions: {
        "elk.algorithm": "layered",
        "elk.direction": "RIGHT",
        "elk.layered.spacing.nodeNodeBetweenLayers": "180",
        "elk.spacing.nodeNode": "120",
        "elk.padding": "[top=56,left=56,bottom=56,right=120]",
      },
      children: graph.nodes.map((stepNode) => {
        const { width, height } = getDimensions(stepNode);
        return {
          id: stepNode.id,
          width,
          height,
        };
      }),
      edges: graph.edges.map((edge) => ({
        id: edge.id,
        sources: [edge.source],
        targets: [edge.target],
      })),
    });

    const layoutNodes = new Map(
      (elkGraph.children ?? []).map((child) => [child.id, child]),
    );

    if (!layoutNodes.size) {
      return buildFallbackFlowGraph(graph, options);
    }

    return {
      nodes: graph.nodes.map((stepNode) => {
        const layoutNode = layoutNodes.get(stepNode.id);
        return {
          id: stepNode.id,
          type: "step",
          position: {
            x: layoutNode?.x ?? 0,
            y: layoutNode?.y ?? 0,
          },
          sourcePosition: Position.Right,
          targetPosition: Position.Left,
          data: {
            stepNode,
            highlighted: options.matchesSearch(stepNode),
            onInspect: options.onInspect,
            onEdit: options.onEdit,
            onPartialRerun: options.onPartialRerun,
          },
          draggable: false,
        };
      }),
      edges: toFlowEdges(graph),
    };
  } catch (error) {
    console.error("ELK layout failed, using fallback placement.", error);
    return buildFallbackFlowGraph(graph, options);
  }
}
