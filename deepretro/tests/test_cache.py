"""Tests for deepretro's explicit cache manager helpers."""

from __future__ import annotations

from deepretro.utils.cache import CacheManager, make_cache_key


def test_make_cache_key_is_stable_and_versioned() -> None:
    """make_cache_key should be deterministic and include versioning."""
    first = make_cache_key("run_az", "CCO", az_model="USPTO", version=1)
    second = make_cache_key("run_az", "CCO", az_model="USPTO", version=1)
    third = make_cache_key("run_az", "CCO", az_model="USPTO", version=2)

    assert first == second
    assert first != third


def test_cache_manager_keeps_values_in_memory(tmp_path, monkeypatch) -> None:
    """CacheManager should not create on-disk state for cached values."""
    cache_dir = tmp_path / "cache"
    monkeypatch.setenv("DEEPRETRO_CACHE_DIR", str(cache_dir))

    cache = CacheManager()
    key = make_cache_key("demo", "CCO", version=1)
    miss = object()

    assert cache.get(key, default=miss) is miss

    cache.set(key, {"smiles": "CCO"}, tag="CCO")

    assert cache.get(key) == {"smiles": "CCO"}
    assert cache.stats().num_entries == 1
    assert not cache_dir.exists()


def test_cache_instances_do_not_share_state() -> None:
    """Separate cache instances should not share cached values implicitly."""
    first = CacheManager()
    second = CacheManager()
    key = make_cache_key("demo", "CCO", version=1)
    miss = object()

    first.set(key, {"smiles": "CCO"})

    assert second.get(key, default=miss) is miss


def test_cache_manager_evict_tag_returns_removed_entries() -> None:
    """evict_tag should remove and count matching in-memory entries."""
    cache = CacheManager()
    first_key = make_cache_key("demo", "CCO", version=1)
    second_key = make_cache_key("demo", "CCN", version=1)
    miss = object()

    cache.set(first_key, {"smiles": "CCO"}, tag="group-a")
    cache.set(second_key, {"smiles": "CCN"}, tag="group-b")

    assert cache.evict_tag("group-a") == 1
    assert cache.get(first_key, default=miss) is miss
    assert cache.get(second_key) == {"smiles": "CCN"}
    assert cache.stats().num_entries == 1
