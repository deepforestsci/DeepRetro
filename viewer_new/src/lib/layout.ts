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
  const moleculeCount =
    stepNode.products.length + stepNode.reactants.length + stepNode.reagents.length;
  const height = Math.min(380, 180 + moleculeCount * 28);

  return {
    width: stepNode.isVirtualRoot ? 300 : 340,
    height,
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
  const elkGraph = await elk.layout({
    id: "deepretro-pathway",
    layoutOptions: {
      "elk.algorithm": "layered",
      "elk.direction": "RIGHT",
      "elk.layered.spacing.nodeNodeBetweenLayers": "120",
      "elk.spacing.nodeNode": "64",
      "elk.padding": "[top=40,left=40,bottom=40,right=80]",
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
    edges: graph.edges.map((edge) => ({
      ...edge,
      animated: false,
      markerEnd: {
        type: MarkerType.ArrowClosed,
      },
      style: {
        stroke: "#48556a",
        strokeWidth: 2,
      },
    })),
  };
}
