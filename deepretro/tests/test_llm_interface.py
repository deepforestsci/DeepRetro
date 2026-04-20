"""Focused tests for DeepRetro LLM provider utilities."""

from __future__ import annotations

import re
from typing import Any

from deepretro.utils.llm_helpers import (
    ChatMessage,
    build_completion_params,
    extract_json_payload,
    infer_provider,
    split_prompt_mode,
)
from deepretro.utils.llm_interface import (
    LLMRequest,
    LLMResponse,
    create_llm_interface,
    parse_cot_response,
)


def test_split_prompt_mode_only_treats_adv_suffix_as_advanced() -> None:
    assert split_prompt_mode("openai/gpt-5:adv") == ("openai/gpt-5", "advanced")
    assert split_prompt_mode("openai/gpt-5") == ("openai/gpt-5", "standard")
    assert split_prompt_mode("fireworks/deepseek-v3p2:adv") == (
        "fireworks/deepseek-v3p2",
        "advanced",
    )
    assert split_prompt_mode("anthropic/claude-opus-4-6:adv") == (
        "anthropic/claude-opus-4-6",
        "advanced",
    )
    assert split_prompt_mode("openai/gpt-5:preview") == (
        "openai/gpt-5:preview",
        "standard",
    )


def test_infer_provider_classifies_supported_model_identifier_shapes() -> None:
    assert infer_provider("openai/gpt-4o-mini") == "openai"
    assert infer_provider("gpt-unlisted-test-model") == "openai"
    assert infer_provider("chatgpt-4o-latest") == "openai"
    assert infer_provider("o3-mini") == "openai"
    assert infer_provider("anthropic/claude-opus-4-6") == "anthropic"
    assert infer_provider("claude-opus-4-6") == "anthropic"
    assert infer_provider("fireworks/deepseek-v3p2") == "deepseek"
    assert infer_provider("deepseek-ai/DeepSeek-R1") == "deepseek"
    assert infer_provider("custom/model") == "generic"


def test_create_llm_interface_selects_provider_parser() -> None:
    assert create_llm_interface("openai/gpt-4o-mini").parse_response(
        '{"data": []}'
    ) == (200, [], '{"data": []}')
    assert create_llm_interface("gpt-unlisted-test-model").parse_response(
        '{"data": []}'
    ) == (200, [], '{"data": []}')
    assert create_llm_interface("claude-opus-4-6").parse_response(
        "<cot><thinking>step</thinking></cot><json>{}</json>"
    ) == (200, ["step"], "{}")
    assert create_llm_interface("fireworks/deepseek-v3p2").parse_response(
        '<think>reason</think><json>{"data": []}</json>'
    ) == (200, ["reason"], '{"data": []}')
    assert create_llm_interface("custom/model").parse_response(
        "<cot><thinking>generic step</thinking></cot><json>{}</json>"
    ) == (200, ["generic step"], "{}")


def test_completion_params_follow_provider_contracts() -> None:
    messages: list[ChatMessage] = [{"role": "user", "content": "Return OK."}]

    openai_params = build_completion_params(
        model="openai/gpt-4o-mini",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.2,
        metadata={"task": "retrosynthesis"},
    )
    assert openai_params == {
        "model": "openai/gpt-4o-mini",
        "messages": messages,
        "max_completion_tokens": 64,
        "temperature": 0.2,
        "seed": 42,
        "metadata": {"task": "retrosynthesis"},
    }

    reasoning = build_completion_params(
        model="openai/gpt-5",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.2,
        thinking_effort="high",
    )
    assert reasoning["max_completion_tokens"] == 8192
    assert reasoning["temperature"] == 1
    assert reasoning["reasoning_effort"] == "high"
    assert "seed" not in reasoning

    anthropic = build_completion_params(
        model="claude-opus-4-6",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.2,
        enable_thinking=False,
    )
    assert anthropic == {
        "model": "claude-opus-4-6",
        "messages": messages,
        "max_tokens": 64,
        "temperature": 0.2,
    }

    anthropic_reasoning = build_completion_params(
        model="anthropic/claude-opus-4-6",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.2,
        enable_thinking=True,
    )
    assert anthropic_reasoning == {
        "model": "anthropic/claude-opus-4-6",
        "messages": messages,
        "max_tokens": 8192,
        "temperature": 1,
        "reasoning_effort": "medium",
    }

    deepseek = build_completion_params("fireworks/deepseek-v3p2", messages, 64, 0.2)
    assert deepseek["model"] == "fireworks_ai/accounts/fireworks/models/deepseek-r1"
    assert deepseek["max_tokens"] == 64
    assert "max_completion_tokens" not in deepseek


def test_interface_build_messages_preserves_custom_messages_and_fills_prompt() -> None:
    custom_messages: list[ChatMessage] = [{"role": "user", "content": "Use this"}]
    request = LLMRequest(
        molecule="CCO",
        model="openai/gpt-4o-mini",
        messages=custom_messages,
    )

    assert (
        create_llm_interface(request.model).build_messages(request) is custom_messages
    )

    messages = create_llm_interface("openai/gpt-4o-mini").build_messages(
        LLMRequest("CCO", "openai/gpt-4o-mini")
    )

    assert [message["role"] for message in messages] == ["system", "user"]
    assert "{target_smiles}" not in messages[1]["content"]
    assert "CCO" in messages[1]["content"]


def test_provider_parsers_extract_payloads_and_preserve_error_codes() -> None:
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
    assert create_llm_interface("openai/gpt-4o-mini").parse_response(
        'prefix {"data": []} suffix'
    ) == (200, [], '{"data": []}')
    assert create_llm_interface("openai/gpt-4o-mini").parse_response(
        "missing json"
    ) == (502, [], "")
    assert create_llm_interface("deepseek-ai/DeepSeek-R1").parse_response(
        deepseek_response
    ) == (200, ["reason"], '{"data": []}')
    assert create_llm_interface("deepseek-ai/DeepSeek-R1").parse_response(
        '{"data": []}'
    ) == (200, [], '{"data": []}')
    assert create_llm_interface("deepseek-ai/DeepSeek-R1").parse_response(
        "missing json"
    ) == (503, [], "")


def test_json_payload_helpers() -> None:
    assert extract_json_payload('<json>{"data": []}</json>') == '{"data": []}'
    assert extract_json_payload('```json\n{"data": []}\n```') == '{"data": []}'
    assert extract_json_payload("prefix [1, 2] suffix") == "[1, 2]"
    assert extract_json_payload("No JSON here") is None


def test_interface_call_wraps_litellm_completion(monkeypatch: Any) -> None:
    captured: dict[str, Any] = {}

    class FakeMessage:
        content = [{"text": "O"}, {"text": "K"}]

    class FakeChoice:
        message = FakeMessage()

    class FakeResponse:
        choices = [FakeChoice()]

    def fake_completion(**params: object) -> FakeResponse:
        captured["params"] = params
        return FakeResponse()

    monkeypatch.setattr(
        "deepretro.utils.llm_interface.completion",
        fake_completion,
    )

    response = create_llm_interface("openai/gpt-4o-mini").call(
        LLMRequest(
            molecule="CCO",
            model="openai/gpt-4o-mini",
            messages=[{"role": "user", "content": "Return OK"}],
            max_output_tokens=32,
        )
    )

    assert response == LLMResponse(status_code=200, text="OK")
    assert captured["params"]["model"] == "openai/gpt-4o-mini"
    assert captured["params"]["messages"] == [{"role": "user", "content": "Return OK"}]
    assert captured["params"]["max_completion_tokens"] == 32


def test_interface_call_returns_error_after_retries(monkeypatch: Any) -> None:
    attempts = 0

    def fake_completion(**params: object) -> object:
        nonlocal attempts
        del params
        attempts += 1
        raise RuntimeError("network unavailable")

    monkeypatch.setattr("deepretro.utils.llm_interface.completion", fake_completion)

    response = create_llm_interface("openai/gpt-4o-mini").call(
        LLMRequest(
            molecule="CCO",
            model="openai/gpt-4o-mini",
            messages=[{"role": "user", "content": "Return OK"}],
        )
    )

    assert attempts == 2
    assert response.status_code == 400
    assert re.search("network unavailable", response.text)
