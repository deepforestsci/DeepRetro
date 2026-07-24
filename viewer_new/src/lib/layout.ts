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
  cardWidth: number;
  cardHeight: number;
  previewWidth: number;
  previewHeight: number;
  onInspect: (stepId: string) => void;
  onEdit: (stepId: string) => void;
  onPartialRerun: (stepId: string) => void;
}>;

const elk = new ELK();

function getMoleculeSizing(stepNode: NormalizedStepNode) {
  const smiles = stepNode.products[0]?.smiles ?? "";
  const atomCount = (smiles.match(/[A-Z][a-z]?/g) ?? []).length;
  const complexity = Math.max(smiles.length, atomCount * 4);

  if (stepNode.isVirtualRoot) {
    if (complexity > 110) {
      return {
        cardWidth: 500,
        cardHeight: 460,
        previewWidth: 420,
        previewHeight: 380,
      };
    }

    return {
      cardWidth: 430,
      cardHeight: 390,
      previewWidth: 360,
      previewHeight: 320,
    };
  }

  if (complexity > 110) {
    return {
      cardWidth: 430,
      cardHeight: 390,
      previewWidth: 360,
      previewHeight: 320,
    };
  }

  if (complexity > 75) {
    return {
      cardWidth: 370,
      cardHeight: 340,
      previewWidth: 310,
      previewHeight: 275,
    };
  }

  return {
    cardWidth: stepNode.isVirtualRoot ? 430 : 290,
    cardHeight: stepNode.isVirtualRoot ? 390 : 290,
    previewWidth: stepNode.isVirtualRoot ? 360 : 235,
    previewHeight: stepNode.isVirtualRoot ? 320 : 210,
  };
}

function toFlowEdges(graph: NormalizedPathwayGraph) {
  // Edge and marker colors are theme-driven from global.css.
  return graph.edges.map((edge) => ({
    ...edge,
    type: "smoothstep" as const,
    animated: false,
    markerEnd: {
      type: MarkerType.ArrowClosed,
      width: 18,
      height: 18,
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
      const sizing = getMoleculeSizing(stepNode);

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
          ...sizing,
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
        "elk.layered.spacing.nodeNodeBetweenLayers": "280",
        "elk.spacing.nodeNode": "210",
        "elk.padding": "[top=72,left=72,bottom=72,right=140]",
      },
      children: graph.nodes.map((stepNode) => {
        const sizing = getMoleculeSizing(stepNode);
        return {
          id: stepNode.id,
          width: sizing.cardWidth,
          height: sizing.cardHeight,
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
        const sizing = getMoleculeSizing(stepNode);
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
            ...sizing,
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
