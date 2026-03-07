import { useEffect, useMemo, useRef, useState } from "react";
import {
  Background,
  Controls,
  type Edge,
  MiniMap,
  ReactFlow,
  ReactFlowProvider,
  useReactFlow,
  type NodeTypes,
} from "@xyflow/react";
import { Expand, Search } from "lucide-react";

import "@xyflow/react/dist/style.css";

import { buildFlowGraph, type StepFlowNode } from "../../lib/layout";
import { StepNode } from "./StepNode";
import type { ViewerRun } from "../../types/viewer";

const nodeTypes = {
  step: StepNode,
} satisfies NodeTypes;

type GraphCanvasProps = {
  run?: ViewerRun;
  selectedStepId?: string;
  onSelectStep: (stepId: string) => void;
  onEditStep: (stepId: string) => void;
  onPartialRerunStep: (stepId: string) => void;
};

function GraphCanvasInner({
  run,
  selectedStepId,
  onSelectStep,
  onEditStep,
  onPartialRerunStep,
}: GraphCanvasProps) {
  const [nodes, setNodes] = useState<StepFlowNode[]>([]);
  const [edges, setEdges] = useState<Edge[]>([]);
  const [search, setSearch] = useState("");
  const shellRef = useRef<HTMLDivElement | null>(null);
  const reactFlow = useReactFlow();

  const searchQuery = search.trim().toLowerCase();

  const matchesSearch = useMemo(() => {
    return (stepNode: NonNullable<ViewerRun["graph"]>["nodes"][number]) => {
      if (!searchQuery) {
        return true;
      }

      const haystacks = [
        stepNode.title,
        ...stepNode.products.map((molecule) => `${molecule.label} ${molecule.smiles}`),
        ...stepNode.reactants.map((molecule) => `${molecule.label} ${molecule.smiles}`),
        ...stepNode.reagents.map((molecule) => `${molecule.label} ${molecule.smiles}`),
      ];

      return haystacks.some((value) => value.toLowerCase().includes(searchQuery));
    };
  }, [searchQuery]);

  useEffect(() => {
    let cancelled = false;

    if (!run?.graph) {
      setNodes([]);
      setEdges([]);
      return () => {
        cancelled = true;
      };
    }

    buildFlowGraph(run.graph, {
      onInspect: onSelectStep,
      onEdit: onEditStep,
      onPartialRerun: onPartialRerunStep,
      matchesSearch,
    }).then((flow) => {
      if (!cancelled) {
        setNodes(flow.nodes);
        setEdges(flow.edges);
      }
    });

    return () => {
      cancelled = true;
    };
  }, [matchesSearch, onEditStep, onPartialRerunStep, onSelectStep, run]);

  useEffect(() => {
    if (!nodes.length) {
      return;
    }
    const timeout = window.setTimeout(() => {
      reactFlow.fitView({
        duration: 350,
        padding: 0.16,
        includeHiddenNodes: true,
      });
    }, 40);

    return () => {
      window.clearTimeout(timeout);
    };
  }, [nodes, reactFlow]);

  if (!run) {
    return (
      <div className="canvas-empty">
        <h3>No pathway selected</h3>
        <p>Run an analysis or upload a pathway JSON to populate the canvas.</p>
      </div>
    );
  }

  if (run.status === "error") {
    return (
      <div className="canvas-empty error">
        <h3>Pathway request failed</h3>
        <p>{run.error}</p>
      </div>
    );
  }

  if (!run.graph) {
    return (
      <div className="canvas-empty">
        <h3>Waiting for graph data</h3>
        <p>The active run has not produced a pathway graph yet.</p>
      </div>
    );
  }

  return (
    <div className="canvas-shell" ref={shellRef}>
      <div className="canvas-toolbar">
        <label className="search-field">
          <Search size={16} />
          <input
            type="search"
            placeholder="Search by step, molecule, or SMILES"
            value={search}
            onChange={(event) => setSearch(event.target.value)}
          />
        </label>
        <button
          className="ghost-button"
          onClick={() => {
            void shellRef.current?.requestFullscreen?.();
          }}
          type="button"
        >
          <Expand size={14} />
          Fullscreen
        </button>
      </div>
      <ReactFlow
        fitView
        nodes={nodes.map((node) => ({
          ...node,
          selected: node.data.stepNode.stepId === selectedStepId,
        }))}
        edges={edges}
        nodeTypes={nodeTypes}
        minZoom={0.2}
        maxZoom={1.8}
        onNodeClick={(_, node) => onSelectStep(node.data.stepNode.stepId)}
        proOptions={{ hideAttribution: true }}
      >
        <MiniMap
          pannable
          zoomable
          nodeColor={(node) =>
            node.id === "step-0" ? "#7ad2ac" : node.selected ? "#ff8663" : "#8fa7ca"
          }
        />
        <Controls position="bottom-left" showInteractive={false} />
        <Background gap={24} size={1.2} color="#213043" />
      </ReactFlow>
    </div>
  );
}

export function GraphCanvas(props: GraphCanvasProps) {
  return (
    <ReactFlowProvider>
      <GraphCanvasInner {...props} />
    </ReactFlowProvider>
  );
}
