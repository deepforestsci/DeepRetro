"""Unit tests for the AiZynthFinder in-process finder cache (no real models)."""

from __future__ import annotations

import pytest

import deepretro.utils.az as az


class _FakePolicy:
    def select(self, name: str) -> None:  # noqa: D401 - trivial stub
        self.selected = name


class _FakeFinder:
    instances = 0

    def __init__(self, configfile: str) -> None:
        _FakeFinder.instances += 1
        self.configfile = configfile
        self.stock = _FakePolicy()
        self.expansion_policy = _FakePolicy()
        self.filter_policy = _FakePolicy()


@pytest.fixture(autouse=True)
def _clear_cache(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(az, "AiZynthFinder", _FakeFinder)
    az._FINDER_CACHE.clear()
    _FakeFinder.instances = 0
    yield
    az._FINDER_CACHE.clear()


def test_get_finder_builds_once_per_config() -> None:
    first = az._get_finder("configA")
    second = az._get_finder("configA")
    assert first is second
    assert _FakeFinder.instances == 1


def test_get_finder_separate_instance_per_config() -> None:
    a = az._get_finder("configA")
    b = az._get_finder("configB")
    assert a is not b
    assert _FakeFinder.instances == 2


def test_get_finder_selects_policies_once() -> None:
    finder = az._get_finder("configA")
    assert finder.stock.selected == "zinc"
    assert finder.expansion_policy.selected == "uspto"
    assert finder.filter_policy.selected == "uspto"


def test_get_finder_raises_without_aizynthfinder(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(az, "AiZynthFinder", None)
    az._FINDER_CACHE.clear()
    with pytest.raises(ImportError, match="deepretro\\[az\\]"):
        az._get_finder("configA")
