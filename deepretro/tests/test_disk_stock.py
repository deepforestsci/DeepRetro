"""Offline checks for exact, disk-backed stock membership."""

from __future__ import annotations

import builtins
import runpy
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from deepretro.utils import disk_stock
from deepretro.utils.disk_stock import SortedInchiKeyQuery

KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
OTHER = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"


def save_stock(tmp_path: Path, keys: list[str]) -> Path:
    """Write a tiny valid sorted stock."""
    path = tmp_path / "stock.npy"
    np.save(path, np.array(sorted(keys), dtype="S27"))
    return path


def test_exact_membership_and_mmap(tmp_path: Path) -> None:
    """Query full identities, preserve duplicates, and keep storage read-only."""
    with SortedInchiKeyQuery(save_stock(tmp_path, [KEY, KEY, OTHER]), 2) as stock:
        assert isinstance(stock._keys, np.memmap)
        assert not stock._keys.flags.writeable
        assert len(stock) == 2
        assert SimpleNamespace(inchi_key=KEY) in stock
        assert stock.cached_search(SimpleNamespace(inchi_key=OTHER))
        assert SimpleNamespace(inchi_key=KEY[:-1] + "M") not in stock
        assert SimpleNamespace(inchi_key="AAAAAAAAAAAAAA-AAAAAAAAAA-A") not in stock
        assert SimpleNamespace(inchi_key="ZZZZZZZZZZZZZZ-ZZZZZZZZZZ-Z") not in stock
        stock.clear_cache()
        assert SimpleNamespace(inchi_key=KEY) in stock


@pytest.mark.parametrize(
    "key", [None, 1, b"x", "", KEY[:-1], KEY + "N", "é" * 27, "A" * 27, KEY.lower()]
)
def test_malformed_key_is_absent(tmp_path: Path, key: Any) -> None:
    """Malformed values cannot be silently truncated into a match."""
    with SortedInchiKeyQuery(save_stock(tmp_path, [KEY]), 1) as stock:
        assert SimpleNamespace(inchi_key=key) not in stock
        assert object() not in stock


def test_empty_stock(tmp_path: Path) -> None:
    """An empty valid array has no members."""
    with SortedInchiKeyQuery(save_stock(tmp_path, []), 0) as stock:
        assert len(stock) == 0
        assert SimpleNamespace(inchi_key=KEY) not in stock


@pytest.mark.parametrize("length", [-1, 2, True, 1.0])
def test_bad_length(tmp_path: Path, length: Any) -> None:
    """Reject invalid producer cardinality metadata."""
    with pytest.raises((TypeError, ValueError), match="length"):
        SortedInchiKeyQuery(save_stock(tmp_path, [KEY]), length)


@pytest.mark.parametrize(
    "data",
    [np.array([[KEY]], dtype="S27"), np.array([KEY], dtype="U27"), np.array([1])],
)
def test_bad_array_metadata(tmp_path: Path, data: np.ndarray) -> None:
    """Reject incompatible array layout before any membership query."""
    path = tmp_path / "bad.npy"
    np.save(path, data)
    with pytest.raises(ValueError, match="one-dimensional"):
        SortedInchiKeyQuery(path, 1)


def test_archive_rejected(tmp_path: Path) -> None:
    """Compressed archives cannot replace a memory-mapped array."""
    path = tmp_path / "bad.npz"
    np.savez(path, keys=np.array([KEY], dtype="S27"))
    with pytest.raises(ValueError, match="memory-mappable"):
        SortedInchiKeyQuery(path, 1)


def test_context_closes_on_exception(tmp_path: Path) -> None:
    """An interrupted consumer releases the map and later queries fail safely."""
    stock = SortedInchiKeyQuery(save_stock(tmp_path, [KEY]), 1)
    mapping = stock._keys._mmap
    with pytest.raises(LookupError):
        with stock:
            raise LookupError("consumer stopped")
    assert mapping.closed
    stock.close()
    assert len(stock) == 1
    with pytest.raises(RuntimeError, match="closed"):
        stock.cached_search(SimpleNamespace(inchi_key=KEY))
    with pytest.raises(RuntimeError, match="closed"):
        stock.__enter__()


def test_import_without_optional_backend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The standalone lookup remains usable when the optional backend is absent."""
    original = builtins.__import__

    def without_backend(name: str, *args: Any, **kwargs: Any) -> Any:
        if name.startswith("aizynthfinder"):
            raise ModuleNotFoundError(name)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_backend)
    namespace = runpy.run_path(str(Path(disk_stock.__file__)))
    with namespace["SortedInchiKeyQuery"](save_stock(tmp_path, [KEY]), 1) as stock:
        assert stock.cached_search(SimpleNamespace(inchi_key=KEY))


def test_invalid_array_releases_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Constructor validation must release an already opened invalid mapping."""
    path = tmp_path / "invalid.npy"
    np.save(path, np.array([[KEY]], dtype="S27"))
    mapping = np.load(path, mmap_mode="r", allow_pickle=False)
    monkeypatch.setattr(disk_stock.np, "load", lambda *args, **kwargs: mapping)
    with pytest.raises(ValueError, match="one-dimensional"):
        SortedInchiKeyQuery(path, 1)
    assert mapping._mmap.closed


def test_optional_backend_interface(tmp_path: Path) -> None:
    """When installed, the backend recognizes the custom query as its stock type."""
    queries = pytest.importorskip("aizynthfinder.context.stock.queries")
    with SortedInchiKeyQuery(save_stock(tmp_path, [KEY]), 1) as stock:
        assert isinstance(stock, queries.StockQueryMixin)
