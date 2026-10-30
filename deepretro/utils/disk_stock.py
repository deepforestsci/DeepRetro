"""Exact stock membership over sorted, read-only memory-mapped InChIKeys."""

from __future__ import annotations

import re
from pathlib import Path
from types import TracebackType
from typing import Any

import numpy as np

try:
    from aizynthfinder.context.stock.queries import StockQueryMixin  # ty: ignore[unresolved-import]
except ImportError:
    # Keep the query usable without the optional template-search dependency.
    class StockQueryMixin:  # type: ignore[no-redef]
        """Fallback base when template-search integration is unavailable."""


_StockQueryBase: Any = StockQueryMixin

_KEY_PATTERN = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]\Z")


class SortedInchiKeyQuery(_StockQueryBase):
    """Query an exact stock without loading its keys into a Python set.

    Parameters
    ----------
    path : str or pathlib.Path
        Sorted, one-dimensional ``.npy`` array of dtype ``S27``. The producer
        must verify sorting and valid keys; opening does not scan the array.
    length : int
        Producer-verified unique cardinality, between zero and the array size.
        Duplicate keys may remain in the array.

    Notes
    -----
    When AiZynthFinder is installed, this class implements its stock query
    interface. Membership is also usable without that optional dependency.
    Close the query when finished, or use it as a context manager.

    Examples
    --------
    >>> with SortedInchiKeyQuery("stock.npy", length=123) as stock:  # doctest: +SKIP
    ...     present = molecule in stock
    """

    def __init__(self, path: str | Path, length: int) -> None:
        """Open read-only key storage and validate its metadata."""
        self._keys: np.memmap | None = None
        if isinstance(length, bool) or not isinstance(length, int):
            raise TypeError("Stock length must be an integer")
        keys: Any = np.load(path, mmap_mode="r", allow_pickle=False)
        try:
            if not isinstance(keys, np.memmap):
                raise ValueError("Stock must be a memory-mappable .npy array")
            if keys.ndim != 1 or keys.dtype != np.dtype("S27"):
                raise ValueError("Stock must be a one-dimensional array of dtype S27")
            if not 0 <= length <= keys.size:
                raise ValueError("Stock length must be between zero and array size")
        except Exception:
            if isinstance(keys, np.memmap):
                keys._mmap.close()
            elif hasattr(keys, "close"):
                keys.close()
            raise
        self._keys = keys
        self._length = length

    def __len__(self) -> int:
        """Return the producer-verified unique cardinality.

        Returns
        -------
        int
            Number of distinct keys, including after the query is closed.
        """
        return self._length

    def __contains__(self, mol: object) -> bool:
        """Look up a molecule's complete InChIKey using binary search.

        Parameters
        ----------
        mol : object
            Molecule exposing an ``inchi_key`` string.

        Returns
        -------
        bool
            Whether the exact 27-character key occurs in the stock. Missing or
            malformed keys return false; stereochemical suffixes are retained.

        Raises
        ------
        RuntimeError
            If the query has been closed.
        """
        if self._keys is None:
            raise RuntimeError("Stock query is closed")
        key = getattr(mol, "inchi_key", None)
        if not isinstance(key, str) or _KEY_PATTERN.fullmatch(key) is None:
            return False
        encoded = key.encode("ascii")
        index = int(np.searchsorted(self._keys, encoded))
        return index < self._keys.size and bytes(self._keys[index]) == encoded

    def cached_search(self, mol: object) -> bool:
        """Query membership without allocating an additional cache.

        Parameters
        ----------
        mol : object
            Molecule exposing an ``inchi_key`` string.

        Returns
        -------
        bool
            Exact stock membership.
        """
        return mol in self

    def clear_cache(self) -> None:
        """Leave the query unchanged because no search cache is allocated."""

    def close(self) -> None:
        """Release the memory mapping; repeated calls are harmless."""
        if self._keys is not None:
            self._keys._mmap.close()
            self._keys = None

    def __enter__(self) -> SortedInchiKeyQuery:
        """Return this open query for a context-managed stock lookup.

        Returns
        -------
        SortedInchiKeyQuery
            This query.

        Raises
        ------
        RuntimeError
            If the query has been closed.
        """
        if self._keys is None:
            raise RuntimeError("Stock query is closed")
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Close the mapping on both normal and exceptional context exit."""
        self.close()
