"""Tests for deepretro's explicit cache manager helpers."""

from __future__ import annotations

from deepretro.utils.cache import CacheEntry, CacheManager, make_args_hash, make_cache_key


def test_make_cache_key_is_stable_and_versioned() -> None:
    """make_cache_key should be deterministic and include versioning."""
    first = make_cache_key("run_az", "CCO", az_model="USPTO", version=1)
    second = make_cache_key("run_az", "CCO", az_model="USPTO", version=1)
    third = make_cache_key("run_az", "CCO", az_model="USPTO", version=2)

    assert first == second
    assert first != third


def test_make_args_hash_is_stable_for_identical_inputs() -> None:
    """make_args_hash should produce stable hashes for the same payload."""
    first = make_args_hash("CCO", az_model="USPTO", include_scores=True)
    second = make_args_hash("CCO", az_model="USPTO", include_scores=True)

    assert first == second


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


def test_cache_manager_purge_if_expired_removes_stale_keys(monkeypatch) -> None:
    """purge_if_expired should drop expired entries and report the removal."""
    cache = CacheManager()
    key = make_cache_key("demo", "CCO", version=1)
    monotonic_values = iter([10.0, 11.5])

    monkeypatch.setattr("deepretro.utils.cache.time.monotonic", lambda: next(monotonic_values))
    cache.set(key, {"smiles": "CCO"}, expire=1.0, tag="alcohol")

    assert cache.purge_if_expired(key) is True
    assert cache.get(key, default=None) is None
    assert cache.stats().num_entries == 0


def test_cache_manager_delete_key_removes_tag_index_entry() -> None:
    """delete_key should remove an entry and its tag association."""
    cache = CacheManager()
    key = make_cache_key("demo", "CCO", version=1)
    cache.set(key, CacheEntry(value={"smiles": "CCO"}, expires_at=None, tag="alcohol"))

    assert cache.delete_key(key) is True
    assert cache.evict_tag("alcohol") == 0
    assert cache.stats().num_entries == 0


def test_cache_manager_evict_tag_returns_removed_entries() -> None:
    """evict_tag should remove every key associated with the same tag."""
    cache = CacheManager()
    first_key = make_cache_key("demo", "CCO", version=1)
    second_key = make_cache_key("demo", "CCN", version=1)
    third_key = make_cache_key("demo", "CCC", version=1)
    miss = object()

    cache.set(first_key, {"smiles": "CCO"}, tag="group-a")
    cache.set(second_key, {"smiles": "CCN"}, tag="group-a")
    cache.set(third_key, {"smiles": "CCC"}, tag="group-b")

    assert cache.evict_tag("group-a") == 2
    assert cache.get(first_key, default=miss) is miss
    assert cache.get(second_key, default=miss) is miss
    assert cache.get(third_key) == {"smiles": "CCC"}
    assert cache.stats().num_entries == 1
