"""Tests for the agent tool registry, executors, and dispatcher."""

from __future__ import annotations

from deepretro.agents.sandbox import SubprocessSandbox
from deepretro.agents.tools import build_tool_registry


def keep_all(product: str, pathways: list) -> tuple[int, list]:
    """Hallucination checker double that keeps every pathway."""
    del product
    return 200, pathways


def reject_all(product: str, pathways: list) -> tuple[int, list]:
    """Hallucination checker double that rejects every pathway."""
    del product, pathways
    return 400, []


class TestValidateSmiles:
    def test_valid_smiles_returns_canonical(self) -> None:
        registry = build_tool_registry()
        result = registry.execute("validate_smiles", {"smiles": "C(O)C"})
        assert result["valid"] is True
        assert result["canonical_smiles"] == "CCO"

    def test_invalid_smiles_reports_error(self) -> None:
        registry = build_tool_registry()
        result = registry.execute("validate_smiles", {"smiles": "not_a_smiles!!!"})
        assert result["valid"] is False
        assert result["error"]


class TestCheckStability:
    def test_returns_assessment(self) -> None:
        registry = build_tool_registry()
        result = registry.execute("check_stability", {"smiles": "CCO"})
        assert "assessment" in result
        assert "stability_score" in result


class TestCheckHallucination:
    def test_registered_when_checker_provided(self) -> None:
        registry = build_tool_registry(hallucination_checker=keep_all)
        assert "check_hallucination" in registry.names

    def test_omitted_when_no_checker(self) -> None:
        registry = build_tool_registry(hallucination_checker=None)
        assert "check_hallucination" not in registry.names

    def test_keep_all_reports_not_hallucination(self) -> None:
        registry = build_tool_registry(hallucination_checker=keep_all)
        result = registry.execute(
            "check_hallucination", {"product": "CCO", "reactants": ["CC", "O"]}
        )
        assert result["is_hallucination"] is False

    def test_reject_all_reports_hallucination(self) -> None:
        registry = build_tool_registry(hallucination_checker=reject_all)
        result = registry.execute(
            "check_hallucination", {"product": "CCO", "reactants": ["CC", "O"]}
        )
        assert result["is_hallucination"] is True


class TestRunPython:
    def test_absent_in_structured_backend(self) -> None:
        registry = build_tool_registry(tool_backend="structured")
        assert "run_python" not in registry.names

    def test_present_and_executes_in_sandbox_backend(self) -> None:
        registry = build_tool_registry(
            tool_backend="sandbox", sandbox=SubprocessSandbox()
        )
        assert "run_python" in registry.names
        result = registry.execute("run_python", {"code": "print(6 * 7)"})
        assert result["ok"] is True
        assert "42" in result["stdout"]


class TestDispatcher:
    def test_unknown_tool_returns_error(self) -> None:
        registry = build_tool_registry()
        result = registry.execute("does_not_exist", {})
        assert "error" in result

    def test_executor_exception_is_captured(self) -> None:
        # Missing required argument raises inside the executor; the dispatcher
        # must turn that into an error dict rather than propagating.
        registry = build_tool_registry()
        result = registry.execute("validate_smiles", {})
        assert "error" in result


class TestSchemas:
    def test_schemas_use_openai_function_format(self) -> None:
        registry = build_tool_registry(
            hallucination_checker=keep_all,
            tool_backend="sandbox",
            sandbox=SubprocessSandbox(),
        )
        assert registry.schemas
        for schema in registry.schemas:
            assert schema["type"] == "function"
            assert "name" in schema["function"]
            assert "parameters" in schema["function"]
        schema_names = {s["function"]["name"] for s in registry.schemas}
        assert schema_names == set(registry.names)


class TestAzTools:
    def test_absent_by_default(self) -> None:
        registry = build_tool_registry()
        assert "az_suggest" not in registry.names
        assert "az_route" not in registry.names

    def test_present_when_enabled(self) -> None:
        registry = build_tool_registry(az_tools=True)
        assert "az_suggest" in registry.names
        assert "az_route" in registry.names

    def test_az_suggest_executes(self, monkeypatch) -> None:
        import deepretro.utils.az as az

        monkeypatch.setattr(
            az,
            "az_single_step",
            lambda smiles, az_model="USPTO", top_k=5: [
                {"precursors": ["CC(=O)OC(C)=O", "O=C(O)c1ccccc1O"], "probability": 0.5}
            ],
        )
        registry = build_tool_registry(az_tools=True)
        result = registry.execute("az_suggest", {"smiles": "CC(=O)Oc1ccccc1C(=O)O"})
        assert result["suggestions"][0]["probability"] == 0.5
        assert result["suggestions"][0]["precursors"][0] == "CC(=O)OC(C)=O"

    def test_az_route_executes(self, monkeypatch) -> None:
        import deepretro.utils.az as az

        monkeypatch.setattr(
            az,
            "run_az",
            lambda smiles, az_model="USPTO": (
                True,
                [{"type": "mol", "smiles": smiles}],
            ),
        )
        registry = build_tool_registry(az_tools=True)
        result = registry.execute("az_route", {"smiles": "CCO"})
        assert result["solved"] is True
        assert result["route"]["smiles"] == "CCO"

    def test_az_tool_error_is_captured(self, monkeypatch) -> None:
        import deepretro.utils.az as az

        def boom(*args, **kwargs):
            raise RuntimeError("az down")

        monkeypatch.setattr(az, "az_single_step", boom)
        registry = build_tool_registry(az_tools=True)
        result = registry.execute("az_suggest", {"smiles": "CCO"})
        assert "error" in result
