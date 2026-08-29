"""Tests for the pure helpers in :mod:`deepretro.utils.llm_helpers`."""

import pytest

from deepretro.utils.llm_helpers import (
    MIN_REASONING_OUTPUT_TOKENS,
    PREFERRED_DEEPSEEK_MODEL,
    coerce_response_text,
    extract_tag_content,
    looks_like_anthropic_reasoning_model,
    looks_like_openai_model,
    looks_like_openai_reasoning_model,
    normalize_completion_model,
    resolve_model_selection,
    resolve_output_token_limit,
    strip_code_fences,
    strip_provider_prefix,
)


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("openai/gpt-4o-mini", "gpt-4o-mini"),
        ("anthropic/claude-opus-4-8", "claude-opus-4-8"),
        ("claude-opus-4-8", "claude-opus-4-8"),
        (
            "fireworks_ai/accounts/fireworks/models/deepseek-r1",
            "fireworks_ai/accounts/fireworks/models/deepseek-r1",
        ),
    ],
)
def test_strip_provider_prefix(model, expected):
    assert strip_provider_prefix(model) == expected


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("openai/gpt-4o-mini", True),
        ("gpt-4o", True),
        ("o3-mini", True),
        ("claude-opus-4-8", False),
        ("fireworks/deepseek-r1", False),
    ],
)
def test_looks_like_openai_model(model, expected):
    assert looks_like_openai_model(model) is expected


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("openai/gpt-5", True),
        ("o1-preview", True),
        ("openai/gpt-4o-mini", False),
        ("claude-opus-4-8", False),
    ],
)
def test_looks_like_openai_reasoning_model(model, expected):
    assert looks_like_openai_reasoning_model(model) is expected


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("anthropic/claude-sonnet-4-6", True),
        ("claude-opus-4-8", True),
        ("claude-fable-5", True),
        ("claude-sonnet-5", True),
        ("claude-3-5-haiku-20241022", False),
        ("openai/gpt-5", False),
    ],
)
def test_looks_like_anthropic_reasoning_model(model, expected):
    assert looks_like_anthropic_reasoning_model(model) is expected


def test_normalize_completion_model_passthrough_non_deepseek():
    assert (
        normalize_completion_model("openai/gpt-4o-mini", "openai")
        == "openai/gpt-4o-mini"
    )


def test_normalize_completion_model_rewrites_deepseek_alias():
    assert (
        normalize_completion_model("fireworks/deepseek-v3p2", "deepseek")
        == PREFERRED_DEEPSEEK_MODEL
    )


def test_output_token_limit_raised_for_reasoning_with_thinking():
    selection = resolve_model_selection("openai/gpt-5")
    result = resolve_output_token_limit(selection, 32, enable_thinking=True)
    assert result == MIN_REASONING_OUTPUT_TOKENS


def test_output_token_limit_preserves_request_when_above_floor():
    selection = resolve_model_selection("openai/gpt-5")
    big = MIN_REASONING_OUTPUT_TOKENS + 1000
    assert resolve_output_token_limit(selection, big, enable_thinking=True) == big


def test_output_token_limit_unchanged_without_thinking():
    selection = resolve_model_selection("openai/gpt-5")
    assert resolve_output_token_limit(selection, 32, enable_thinking=False) == 32


@pytest.mark.parametrize(
    ("content", "expected"),
    [
        ("OK", "OK"),
        (None, ""),
        ([{"text": "O"}, {"text": "K"}], "OK"),
        ([{"text": "A"}, {"no_text": "B"}], "A"),
        (123, "123"),
    ],
)
def test_coerce_response_text(content, expected):
    assert coerce_response_text(content) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ('```json\n{"data": []}\n```', '{"data": []}'),
        ("```\nplain\n```", "plain"),
        ('{"data": []}', '{"data": []}'),
    ],
)
def test_strip_code_fences(text, expected):
    assert strip_code_fences(text) == expected


def test_extract_tag_content_found():
    assert extract_tag_content("<json>{}</json>", "json") == "{}"


def test_extract_tag_content_missing_returns_none():
    assert extract_tag_content("no tags here", "json") is None
