"""Focused tests for DeepRetro public LLM utilities."""

from __future__ import annotations

import json
import os
import re
from typing import Any

import pytest

from deepretro.utils import llm
from deepretro.utils.llm import (
    build_messages,
    call_LLM,
    llm_pipeline,
    parse_cot_response,
    parse_deepseek_response,
    parse_openai_response,
    parse_response,
    split_cot_json,
    split_json_deepseek,
    split_json_master,
    split_json_openAI,
    validate_split_json,
)
from deepretro.utils.llm_helpers import ChatMessage, PromptMode
from deepretro.utils.llm_interface import LLMRequest, LLMResponse

LIVE_TEST_MOLECULE = "CC(=O)Oc1ccccc1C(=O)O"


def test_build_messages_preserves_custom_messages_and_fills_default_prompt() -> None:
    custom_messages: list[ChatMessage] = [{"role": "user", "content": "Use this"}]

    assert (
        build_messages("CCO", "openai/gpt-4o-mini", messages=custom_messages)
        is custom_messages
    )

    messages = build_messages("CCO", "openai/gpt-4o-mini", addon="Extra constraint")

    assert [message["role"] for message in messages] == ["system", "user"]
    assert "{target_smiles}" not in messages[1]["content"]
    assert "CCO" in messages[1]["content"]
    assert messages[1]["content"].endswith("Extra constraint")


def test_public_parsers_extract_payloads_and_preserve_error_codes() -> None:
    cot_response = "<cot><thinking>step</thinking></cot><json>{}</json>"
    cot_with_fenced_json = (
        '<cot><thinking>step</thinking></cot>\n```json\n{"data": []}\n```'
    )
    deepseek_response = '<think>reason</think><json>{"data": []}</json>'

    assert parse_cot_response(cot_response) == (200, ["step"], "{}")
    assert parse_cot_response(cot_with_fenced_json) == (
        200,
        ["step"],
        '{"data": []}',
    )
    assert parse_cot_response("missing tags") == (501, [], "")
    assert parse_openai_response('prefix {"data": []} suffix') == (
        200,
        [],
        '{"data": []}',
    )
    assert parse_openai_response("missing json") == (502, [], "")
    assert parse_deepseek_response(deepseek_response) == (
        200,
        ["reason"],
        '{"data": []}',
    )
    assert parse_deepseek_response('{"data": []}') == (200, [], '{"data": []}')
    assert parse_deepseek_response("missing json") == (503, [], "")


def test_backward_compatible_parser_wrappers_delegate_to_provider_parsers() -> None:
    cot_response = "<cot><thinking>step</thinking></cot><json>{}</json>"
    deepseek_response = '<think>reason</think><json>{"data": []}</json>'

    assert split_cot_json(cot_response) == (200, ["step"], "{}")
    assert split_json_openAI('```json\n{"data": []}\n```') == (200, '{"data": []}')
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
    assert parse_response(cot_response, "claude-opus-4-6") == (200, ["step"], "{}")


def test_validate_split_json_returns_aligned_typed_pipeline_data() -> None:
    payload = (
        '{"data": [["CCO"], ["CCN"]], "explanation": ["ethanol", "ethylamine"], '
        '"confidence_scores": [1, "0.5"]}'
    )

    assert validate_split_json(payload) == (
        200,
        [["CCO"], ["CCN"]],
        ["ethanol", "ethylamine"],
        [1.0, 0.5],
    )
    assert validate_split_json("{bad json") == (504, [], [], [])
    assert validate_split_json('{"data": [["CCO"]]}') == (504, [], [], [])


def test_call_llm_builds_request_for_selected_interface(monkeypatch: Any) -> None:
    captured: dict[str, LLMRequest] = {}

    class FakeInterface:
        def call(self, request: LLMRequest) -> LLMResponse:
            captured["request"] = request
            return LLMResponse(status_code=200, text="OK")

    monkeypatch.setattr(
        llm, "create_llm_interface", lambda *args, **kwargs: FakeInterface()
    )

    status, text = call_LLM(
        molecule="CCO",
        model="openai/gpt-4o-mini",
        messages=[{"role": "user", "content": "Return OK"}],
        prompt_mode="advanced",
        enable_thinking=False,
        max_output_tokens=32,
    )

    assert (status, text) == (200, "OK")
    assert captured["request"] == LLMRequest(
        molecule="CCO",
        model="openai/gpt-4o-mini",
        messages=[{"role": "user", "content": "Return OK"}],
        temperature=0.0,
        prompt_mode="advanced",
        enable_thinking=False,
        max_output_tokens=32,
    )


def assert_retro_response_is_tagged_and_valid(
    response_text: str,
    model: str,
) -> None:
    assert "<cot>" in response_text
    assert "</cot>" in response_text
    assert re.search(r"<thinking\b[^>]*>", response_text)
    assert "</thinking>" in response_text

    parse_status, thinking_steps, json_content = parse_response(response_text, model)
    assert parse_status == 200, response_text[:500]
    assert thinking_steps

    parsed_json = json.loads(json_content)
    assert {"data", "explanation", "confidence_scores"} <= set(parsed_json)

    validation_status, pathways, explanations, confidence = validate_split_json(
        json_content
    )
    assert validation_status == 200
    assert pathways
    assert len(pathways) == len(explanations) == len(confidence)
    assert all(isinstance(pathway, list) and pathway for pathway in pathways)
    assert all(
        isinstance(smiles, str) and smiles for pathway in pathways for smiles in pathway
    )
    assert all(0 <= score <= 1 for score in confidence)


@pytest.mark.slow
@pytest.mark.parametrize("prompt_mode", ["standard", "advanced"])
def test_live_anthropic_retrosynthesis_prompt_returns_parseable_tagged_json(
    prompt_mode: PromptMode,
) -> None:
    if not os.getenv("ANTHROPIC_API_KEY"):
        pytest.skip("ANTHROPIC_API_KEY is not configured")

    model = os.getenv(
        "DEEPRETRO_TEST_ANTHROPIC_MODEL",
        "anthropic/claude-sonnet-4-6:adv",
    )

    status, response_text = call_LLM(
        molecule=LIVE_TEST_MOLECULE,
        model=model,
        prompt_mode=prompt_mode,
        max_output_tokens=8192,
        enable_thinking=False,
    )

    assert status == 200, response_text[:500]
    assert_retro_response_is_tagged_and_valid(response_text, model)


def test_llm_pipeline_falls_back_from_deepseek_after_failed_attempt(
    monkeypatch: Any,
) -> None:
    created_models: list[str] = []

    class FakeInterface:
        def __init__(self, model: str) -> None:
            self.model = model

        def call(self, request: LLMRequest) -> LLMResponse:
            if request.model == "deepseek-ai/DeepSeek-R1":
                return LLMResponse(status_code=400, text="unavailable")
            return LLMResponse(
                status_code=200,
                text="<cot><thinking>step</thinking></cot>"
                '<json>{"data": [["CCO"]], "explanation": ["fallback"], '
                '"confidence_scores": [0.7]}</json>',
            )

        def parse_response(self, response_text: str) -> tuple[int, list[str], str]:
            assert response_text.startswith("<cot>")
            return (
                200,
                ["step"],
                '{"data": [["CCO"]], "explanation": ["fallback"], '
                '"confidence_scores": [0.7]}',
            )

    def fake_create_interface(model: str, **kwargs: object) -> FakeInterface:
        del kwargs
        created_models.append(model)
        return FakeInterface(model)

    monkeypatch.setattr(llm, "create_llm_interface", fake_create_interface)
    monkeypatch.setattr(
        llm,
        "validity_check",
        lambda molecule, pathways, explanations, confidence: (
            pathways,
            explanations,
            confidence,
        ),
    )

    assert llm_pipeline(
        molecule="CCO",
        model="deepseek-ai/DeepSeek-R1",
        max_output_tokens=128,
        enable_thinking=False,
    ) == ([["CCO"]], ["fallback"], [0.7])
    assert created_models[:2] == ["deepseek-ai/DeepSeek-R1", llm.FALLBACK_MODEL]
