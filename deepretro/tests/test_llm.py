"""Tests for DeepRetro LLM utilities."""

from __future__ import annotations

import inspect
import os
from collections.abc import Callable

import pytest
from dotenv import load_dotenv

from deepretro.utils.llm import (
    build_addon_prompt,
    build_messages,
    call_LLM,
    filter_hallucination_pathways,
    filter_metadata_by_pathways,
    filter_stable_pathways,
    llm_pipeline,
    obtain_prompt,
    parse_cot_response,
    parse_deepseek_response,
    parse_openai_response,
    parse_response,
    split_cot_json,
    split_json_deepseek,
    split_json_master,
    split_json_openAI,
    validate_json_response,
    validate_split_json,
)
from deepretro.utils.llm_helpers import (
    ChatMessage,
    build_completion_params,
    coerce_response_text,
    extract_json_payload,
    extract_tag_content,
    infer_provider,
    is_enabled,
    looks_like_openai_model,
    looks_like_openai_reasoning_model,
    normalize_completion_model,
    resolve_model_selection,
    resolve_output_token_limit,
    split_prompt_mode,
    strip_code_fences,
    strip_provider_prefix,
)

PUBLIC_LLM_FUNCTIONS: tuple[Callable[..., object], ...] = (
    call_LLM,
    llm_pipeline,
    parse_response,
    validate_json_response,
)


def test_openai_model_selection_uses_chat_completion_params() -> None:
    """OpenAI chat models should use OpenAI-compatible completion settings."""
    selection = resolve_model_selection("openai/gpt-4o-mini")

    assert selection.provider == "openai"
    assert selection.family == "openai"
    assert selection.completion_model == "openai/gpt-4o-mini"
    assert selection.output_token_param == "max_completion_tokens"
    assert selection.supports_seed is True
    assert selection.supports_reasoning_effort is False


def test_model_helper_functions_classify_and_normalize_models() -> None:
    """Helper functions should classify providers and normalize model aliases."""
    base_model, prompt_mode = split_prompt_mode("openai/gpt-4o-mini:adv")
    deepseek_model = normalize_completion_model("fireworks/deepseek-v3p2", "deepseek")

    assert (base_model, prompt_mode) == ("openai/gpt-4o-mini", "advanced")
    assert split_prompt_mode("claude-opus-4-6") == ("claude-opus-4-6", "standard")
    assert strip_provider_prefix("openai/gpt-4o-mini") == "gpt-4o-mini"
    assert strip_provider_prefix("deepseek-ai/DeepSeek-R1") == "deepseek-ai/DeepSeek-R1"
    assert looks_like_openai_model("gpt-4o-mini") is True
    assert looks_like_openai_model("claude-opus-4-6") is False
    assert looks_like_openai_reasoning_model("openai/gpt-5") is True
    assert looks_like_openai_reasoning_model("openai/gpt-4o-mini") is False
    assert infer_provider("gpt-unlisted-test-model") == "openai"
    assert infer_provider("chatgpt-4o-latest") == "openai"
    assert infer_provider("openai/gpt-4o-mini") == "openai"
    assert infer_provider("claude-opus-4-6") == "anthropic"
    assert infer_provider("fireworks_ai/accounts/fireworks/models/deepseek-r1") == (
        "deepseek"
    )
    assert infer_provider("custom/model") == "generic"
    assert deepseek_model == "fireworks_ai/accounts/fireworks/models/deepseek-r1"
    assert (
        normalize_completion_model(
            "fireworks_ai/accounts/fireworks/models/deepseek-r1",
            "deepseek",
        )
        == "fireworks_ai/accounts/fireworks/models/deepseek-r1"
    )


def test_openai_reasoning_model_drops_seed_and_sets_reasoning_effort() -> None:
    """OpenAI reasoning models should receive reasoning controls, not seeds."""
    messages: list[ChatMessage] = [{"role": "user", "content": "Return OK."}]

    params = build_completion_params(
        model="openai/gpt-5",
        messages=messages,
        max_completion_tokens=32,
        temperature=0.2,
        enable_thinking=True,
        thinking_effort="medium",
    )

    assert params["model"] == "openai/gpt-5"
    assert params["max_completion_tokens"] == 8192
    assert params["temperature"] == 1
    assert params["reasoning_effort"] == "medium"
    assert "seed" not in params


def test_completion_params_cover_seed_metadata_and_token_branches() -> None:
    """Completion parameter construction should cover provider-specific branches."""
    messages: list[ChatMessage] = [{"role": "user", "content": "Return OK."}]
    openai_selection = resolve_model_selection("openai/gpt-4o-mini")
    anthropic_params = build_completion_params(
        model="claude-opus-4-6",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.2,
        enable_thinking=False,
        metadata={"task": "demo"},
    )

    assert resolve_output_token_limit(openai_selection, 32, enable_thinking=True) == 32
    assert anthropic_params["model"] == "claude-opus-4-6"
    assert anthropic_params["max_tokens"] == 64
    assert anthropic_params["temperature"] == 0.2
    assert anthropic_params["metadata"] == {"task": "demo"}
    assert "seed" not in anthropic_params
    assert "reasoning_effort" not in anthropic_params


def test_response_text_and_json_helpers_cover_supported_payload_shapes() -> None:
    """Text and JSON helpers should handle strings, lists, fences, tags, and raw JSON."""
    assert coerce_response_text("OK") == "OK"
    assert coerce_response_text(None) == ""
    assert coerce_response_text([{"text": "O"}, {"text": "K"}, {"ignored": "x"}]) == (
        "OK"
    )
    assert coerce_response_text({"text": "OK"}) == "{'text': 'OK'}"
    assert strip_code_fences('```json\n{"data": []}\n```') == '{"data": []}'
    assert strip_code_fences('{"data": []}') == '{"data": []}'
    assert extract_tag_content("<json>{}</json>", "json") == "{}"
    assert extract_tag_content("missing", "json") is None
    assert extract_json_payload('<json>{"data": []}</json>') == '{"data": []}'
    assert extract_json_payload('```json\n{"data": []}\n```') == '{"data": []}'
    assert extract_json_payload('prefix {"data": []} suffix') == '{"data": []}'
    assert extract_json_payload("prefix [1, 2] suffix") == "[1, 2]"
    assert extract_json_payload("No JSON here") is None
    assert is_enabled(True) is True
    assert is_enabled(False) is False
    assert is_enabled("True") is True
    assert is_enabled("false") is False


def test_prompt_helpers_build_default_custom_and_addon_messages() -> None:
    """Prompt helpers should select prompts and preserve custom messages."""
    custom_messages: list[ChatMessage] = [{"role": "user", "content": "Reply OK"}]
    default_messages = build_messages("CCO", "openai/gpt-4o-mini", addon="ADDON")
    openai_prompt = obtain_prompt("openai/gpt-4o-mini")
    advanced_prompt = obtain_prompt("claude-opus-4-6", prompt_mode="advanced")

    assert openai_prompt[2] == 16384
    assert advanced_prompt[2] == 16384
    assert build_addon_prompt("CCO") == ""
    assert build_addon_prompt("C1CCCCCC1")
    assert default_messages[0]["role"] == "system"
    assert default_messages[1]["role"] == "user"
    assert default_messages[1]["content"].endswith("ADDON")
    assert (
        build_messages("CCO", "openai/gpt-4o-mini", messages=custom_messages)
        is custom_messages
    )


def test_parse_response_accepts_openai_tagged_json() -> None:
    """OpenAI responses should parse JSON payloads without chain-of-thought."""
    response = """
    Here is the answer:
    <json>
    {"data": [["CCO"]], "explanation": ["ethanol"], "confidence_scores": [0.8]}
    </json>
    """

    status, thinking_steps, json_content = parse_response(
        response, "openai/gpt-4o-mini"
    )

    assert status == 200
    assert thinking_steps == []
    assert '"data"' in json_content


def test_parser_wrappers_cover_success_and_error_paths() -> None:
    """Provider parser wrappers should cover tagged, raw, and invalid responses."""
    cot_response = "<cot><thinking>step</thinking></cot><json>{}</json>"
    deepseek_response = '<think>reason</think><json>{"data": []}</json>'

    assert parse_cot_response(cot_response) == (200, ["step"], "{}")
    assert parse_cot_response("missing") == (501, [], "")
    assert parse_cot_response("<cot>no thinking</cot><json>{}</json>") == (
        501,
        [],
        "",
    )
    assert parse_openai_response('{"data": []}') == (200, [], '{"data": []}')
    assert parse_openai_response("missing") == (502, [], "")
    assert parse_deepseek_response(deepseek_response) == (
        200,
        ["reason"],
        '{"data": []}',
    )
    assert parse_deepseek_response('{"data": []}') == (200, [], '{"data": []}')
    assert parse_deepseek_response("missing") == (503, [], "")
    assert split_cot_json(cot_response) == (200, ["step"], "{}")
    assert split_json_openAI('{"data": []}') == (200, '{"data": []}')
    assert split_json_deepseek(deepseek_response) == (
        200,
        ["reason"],
        '{"data": []}',
    )
    assert split_json_master('{"data": []}', "openai/gpt-4o-mini") == (
        200,
        [],
        '{"data": []}',
    )
    assert parse_response(deepseek_response, "deepseek-ai/DeepSeek-R1")[0] == 200
    assert parse_response(cot_response, "claude-opus-4-6") == (200, ["step"], "{}")


def test_validate_json_response_returns_typed_pathways() -> None:
    """JSON payload validation should coerce confidence values to floats."""
    status, pathways, explanations, confidence = validate_json_response(
        '{"data": [["CCO", "O"]], "explanation": ["demo"], "confidence_scores": [1]}'
    )

    assert status == 200
    assert pathways == [["CCO", "O"]]
    assert explanations == ["demo"]
    assert confidence == [1.0]


def test_json_validation_and_filter_helpers_cover_valid_and_invalid_payloads() -> None:
    """Validation and pathway filters should keep metadata aligned."""
    json_content = (
        '{"data": [["CCO"], ["CCN"]], "explanation": ["a", "b"], '
        '"confidence_scores": [0.1, 0.9]}'
    )

    assert validate_split_json(json_content) == (
        200,
        [["CCO"], ["CCN"]],
        ["a", "b"],
        [0.1, 0.9],
    )
    assert validate_json_response("{bad json") == (504, [], [], [])
    assert filter_stable_pathways([["CCO"]]) == [["CCO"]]
    assert filter_hallucination_pathways("CCO", [["CCO"]]) == [["CCO"]]
    assert filter_metadata_by_pathways(
        retained_pathways=[["CCN"], ["CCN"]],
        pathways=[["CCO"], ["CCN"], ["CCN"]],
        explanations=["a", "b", "c"],
        confidence=[0.1, 0.8, 0.9],
    ) == ([["CCN"], ["CCN"]], ["b", "c"], [0.8, 0.9])


def test_llm_pipeline_returns_empty_result_when_real_litellm_call_fails() -> None:
    """Pipeline should return empty outputs when a real LiteLLM call fails."""
    assert llm_pipeline(
        molecule="CCO",
        model="not-a-real-provider/not-a-real-model",
        messages=[{"role": "user", "content": "Reply OK"}],
        max_output_tokens=1,
        enable_thinking=False,
    ) == ([], [], [])


@pytest.mark.slow
def test_llm_pipeline_anthropic_endpoint_returns_pathway() -> None:
    """Call a real Anthropic endpoint through llm_pipeline."""
    load_dotenv()
    if not os.getenv("ANTHROPIC_API_KEY"):
        pytest.skip("ANTHROPIC_API_KEY is not configured")

    messages: list[ChatMessage] = [
        {
            "role": "user",
            "content": (
                "Return exactly this payload and no extra text:\n"
                "<cot><thinking>forced test response</thinking></cot>"
                '<json>{"data": [["N#N"]], "explanation": ["test route"], '
                '"confidence_scores": [0.99]}</json>'
            ),
        }
    ]

    pathways, explanations, confidence = llm_pipeline(
        molecule="CCO",
        model=os.getenv(
            "DEEPRETRO_TEST_ANTHROPIC_MODEL",
            "claude-haiku-4-5-20251001",
        ),
        messages=messages,
        max_output_tokens=512,
        enable_thinking=False,
    )

    if not pathways:
        pytest.skip(
            "Anthropic endpoint did not return a valid pipeline output with the "
            "configured key/model"
        )

    assert pathways == [["N#N"]]
    assert explanations == ["test route"]
    assert confidence == [0.99]


def test_llm_module_has_no_agentic_or_protecting_group_public_switches() -> None:
    """The PR should not expose agentic or protecting-group LLM switches."""
    parameters = inspect.signature(llm_pipeline).parameters

    assert "agentic_mode" not in parameters
    assert "use_protecting_group_feature" not in parameters


def test_llm_module_does_not_import_legacy_src_or_protecting_group_code() -> None:
    """Packaged LLM utilities should not depend on legacy ``src`` modules."""
    import deepretro.utils.llm as llm_module

    source = inspect.getsource(llm_module)

    assert "from src." not in source
    assert "protecting_group" not in source


def test_public_llm_functions_have_annotations_and_usage_examples() -> None:
    """Public LLM functions should meet DeepRetro docstring standards."""
    import deepretro.utils.llm as llm_module
    import deepretro.utils.llm_helpers as helper_module

    public_functions = [
        obj
        for module in (llm_module, helper_module)
        for name, obj in inspect.getmembers(module, inspect.isfunction)
        if not name.startswith("_") and obj.__module__ == module.__name__
    ]

    assert public_functions
    for function in public_functions:
        signature = inspect.signature(function)
        docstring = inspect.getdoc(function)

        assert docstring, f"{function.__name__} is missing a docstring"
        assert "Examples\n--------" in docstring, (
            f"{function.__name__} is missing usage examples"
        )
        assert "Returns\n-------" in docstring, (
            f"{function.__name__} is missing a Returns section"
        )
        assert signature.return_annotation is not inspect.Signature.empty, (
            f"{function.__name__} is missing a return type annotation"
        )
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, (
                f"{function.__name__}.{parameter.name} is missing a type annotation"
            )


@pytest.mark.slow
def test_call_llm_openai_endpoint_returns_text() -> None:
    """Call a real OpenAI endpoint through LiteLLM; no LLM mocks are used."""
    load_dotenv()
    if not os.getenv("OPENAI_API_KEY"):
        pytest.skip("OPENAI_API_KEY is not configured")

    status, response_text = call_LLM(
        molecule="CCO",
        model=os.getenv("DEEPRETRO_TEST_OPENAI_MODEL", "openai/gpt-4o-mini"),
        temperature=0.0,
        messages=[{"role": "user", "content": "Reply with exactly: OK"}],
        max_output_tokens=16,
        enable_thinking=False,
    )

    if status != 200 and "AuthenticationError" in response_text:
        pytest.skip("OpenAI endpoint rejected the configured OPENAI_API_KEY")

    assert status == 200
    assert response_text.strip()
