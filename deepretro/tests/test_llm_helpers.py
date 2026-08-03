"""Unit tests for the pure helpers in :mod:`deepretro.utils.llm_helpers`.

Focuses on the provider/model classification and response-parsing helpers that
were not already exercised by ``test_llm_interface.py`` / ``test_llm.py``.
All targets here are pure functions, so the tests are fast and deterministic.
"""

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

# ---------------------------------------------------------------------------
# strip_provider_prefix
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("openai/gpt-4o-mini", "gpt-4o-mini"),
        ("anthropic/claude-opus-4-8", "claude-opus-4-8"),
        ("claude-opus-4-8", "claude-opus-4-8"),  # no prefix
        # Unknown providers are left intact.
        ("fireworks_ai/accounts/fireworks/models/deepseek-r1",
         "fireworks_ai/accounts/fireworks/models/deepseek-r1"),
    ],
)
def test_strip_provider_prefix(model, expected):
    """Only known ``openai/``/``anthropic/`` prefixes are stripped."""
    assert strip_provider_prefix(model) == expected


# ---------------------------------------------------------------------------
# OpenAI / Anthropic model classification
# ---------------------------------------------------------------------------


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
    """OpenAI-style identifiers are recognised via their base name."""
    assert looks_like_openai_model(model) is expected


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("openai/gpt-5", True),
        ("o1-preview", True),
        ("openai/gpt-4o-mini", False),  # non-reasoning OpenAI model
        ("claude-opus-4-8", False),
    ],
)
def test_looks_like_openai_reasoning_model(model, expected):
    """Only OpenAI reasoning families (o1/o3/o4, gpt-5) return True."""
    assert looks_like_openai_reasoning_model(model) is expected


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("anthropic/claude-sonnet-4-6", True),
        ("claude-opus-4-8", True),
        ("claude-fable-5", True),
        ("claude-sonnet-5", True),
        ("claude-3-5-haiku-20241022", False),  # Claude 3.x is not a reasoning family
        ("openai/gpt-5", False),
    ],
)
def test_looks_like_anthropic_reasoning_model(model, expected):
    """Claude 4.x/5 reasoning families are recognised; Claude 3.x is not."""
    assert looks_like_anthropic_reasoning_model(model) is expected


# ---------------------------------------------------------------------------
# normalize_completion_model
# ---------------------------------------------------------------------------


def test_normalize_completion_model_passthrough_non_deepseek():
    """Non-DeepSeek providers are returned unchanged."""
    assert (
        normalize_completion_model("openai/gpt-4o-mini", "openai")
        == "openai/gpt-4o-mini"
    )


def test_normalize_completion_model_rewrites_deepseek_alias():
    """A DeepSeek alias is normalised to the preferred LiteLLM identifier."""
    assert (
        normalize_completion_model("fireworks/deepseek-v3p2", "deepseek")
        == PREFERRED_DEEPSEEK_MODEL
    )


# ---------------------------------------------------------------------------
# resolve_output_token_limit
# ---------------------------------------------------------------------------


def test_output_token_limit_raised_for_reasoning_with_thinking():
    """A reasoning model with thinking on gets at least the reasoning floor."""
    selection = resolve_model_selection("openai/gpt-5")
    result = resolve_output_token_limit(selection, 32, enable_thinking=True)
    assert result == MIN_REASONING_OUTPUT_TOKENS


def test_output_token_limit_preserves_request_when_above_floor():
    """A request already above the floor is preserved."""
    selection = resolve_model_selection("openai/gpt-5")
    big = MIN_REASONING_OUTPUT_TOKENS + 1000
    assert resolve_output_token_limit(selection, big, enable_thinking=True) == big


def test_output_token_limit_unchanged_without_thinking():
    """With thinking disabled the requested limit is returned as-is."""
    selection = resolve_model_selection("openai/gpt-5")
    assert resolve_output_token_limit(selection, 32, enable_thinking=False) == 32


# ---------------------------------------------------------------------------
# Response parsing: coerce_response_text / strip_code_fences / extract_tag_content
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("content", "expected"),
    [
        ("OK", "OK"),
        (None, ""),
        ([{"text": "O"}, {"text": "K"}], "OK"),
        ([{"text": "A"}, {"no_text": "B"}], "A"),  # non-text blocks ignored
        (123, "123"),  # fallback to str()
    ],
)
def test_coerce_response_text(content, expected):
    """Assorted LiteLLM content shapes coerce to a plain string."""
    assert coerce_response_text(content) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ('```json\n{"data": []}\n```', '{"data": []}'),
        ("```\nplain\n```", "plain"),
        ('{"data": []}', '{"data": []}'),  # no fence, just stripped
    ],
)
def test_strip_code_fences(text, expected):
    """A single surrounding Markdown fence is removed."""
    assert strip_code_fences(text) == expected


def test_extract_tag_content_found():
    """Content between matching tags is returned."""
    assert extract_tag_content("<json>{}</json>", "json") == "{}"


def test_extract_tag_content_missing_returns_none():
    """A missing tag yields None."""
    assert extract_tag_content("no tags here", "json") is None
