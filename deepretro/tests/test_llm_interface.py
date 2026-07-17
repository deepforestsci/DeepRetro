"""Focused tests for DeepRetro LLM provider utilities."""

from __future__ import annotations

import re
from typing import Any

from deepretro.utils.llm_helpers import (
    ChatMessage,
    apply_prompt_caching,
    build_completion_params,
    extract_json_payload,
    infer_provider,
    resolve_model_selection,
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
    assert infer_provider("custom/model") == "anthropic"


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

    # added seperately to test reasoning effort, reasoning requires temperature=1
    # max_tokens should also be set to something higher for reasoning models
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

    # anthropic models use max_tokens instead of max_completion_tokens
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

    # added seperately to test reasoning effort, reasoning requires temperature=1,
    # max_tokens should also be set to something higher for reasoning models
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



def test_prompt_caching_marks_system_and_final_user_for_anthropic(
    monkeypatch: Any,
) -> None:
    monkeypatch.delenv("DEEPRETRO_PROMPT_CACHING", raising=False)
    selection = resolve_model_selection("claude-opus-4-6")
    messages: list[ChatMessage] = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "target molecule"},
    ]

    cached = apply_prompt_caching(messages, selection)

    assert cached[0]["content"] == [
        {
            "type": "text",
            "text": "system prompt",
            "cache_control": {"type": "ephemeral"},
        }
    ]
    assert cached[1]["content"] == [
        {
            "type": "text",
            "text": "target molecule",
            "cache_control": {"type": "ephemeral"},
        }
    ]
    # The original messages are not mutated.
    assert messages[0]["content"] == "system prompt"
    assert messages[1]["content"] == "target molecule"


def test_prompt_caching_marks_only_system_for_multi_turn_tool_conversation(
    monkeypatch: Any,
) -> None:
    monkeypatch.delenv("DEEPRETRO_PROMPT_CACHING", raising=False)
    selection = resolve_model_selection("claude-opus-4-6")
    messages: list[Any] = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "target molecule"},
        {"role": "assistant", "content": "", "tool_calls": [{"id": "1"}]},
        {"role": "tool", "tool_call_id": "1", "content": "tool result"},
    ]

    cached = apply_prompt_caching(messages, selection)

    assert cached[0]["content"][0]["cache_control"] == {"type": "ephemeral"}
    # The trailing tool message keeps its plain-string content (no breakpoint).
    assert cached[3]["content"] == "tool result"
    assert isinstance(cached[1]["content"], str)


def test_prompt_caching_is_noop_for_non_anthropic_providers() -> None:
    messages: list[ChatMessage] = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "target molecule"},
    ]
    for model in ("openai/gpt-4o-mini", "deepseek-ai/DeepSeek-R1"):
        selection = resolve_model_selection(model)
        assert apply_prompt_caching(messages, selection) is messages


def test_prompt_caching_is_noop_without_system_message() -> None:
    selection = resolve_model_selection("claude-opus-4-6")
    messages: list[ChatMessage] = [{"role": "user", "content": "target molecule"}]
    assert apply_prompt_caching(messages, selection) is messages


def test_prompt_caching_can_be_disabled_by_env(monkeypatch: Any) -> None:
    monkeypatch.setenv("DEEPRETRO_PROMPT_CACHING", "off")
    selection = resolve_model_selection("claude-opus-4-6")
    messages: list[ChatMessage] = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "target molecule"},
    ]
    assert apply_prompt_caching(messages, selection) is messages


def test_build_completion_params_adds_cache_control_for_anthropic(
    monkeypatch: Any,
) -> None:
    monkeypatch.delenv("DEEPRETRO_PROMPT_CACHING", raising=False)
    messages: list[ChatMessage] = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "target molecule"},
    ]

    anthropic_params = build_completion_params(
        model="claude-opus-4-6",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.0,
        enable_thinking=False,
    )
    system_block = anthropic_params["messages"][0]["content"][0]
    assert system_block["cache_control"] == {"type": "ephemeral"}

    openai_params = build_completion_params(
        model="openai/gpt-4o-mini",
        messages=messages,
        max_completion_tokens=64,
        temperature=0.0,
    )
    # OpenAI requests are untouched (string content, no cache_control).
    assert openai_params["messages"] == messages


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


