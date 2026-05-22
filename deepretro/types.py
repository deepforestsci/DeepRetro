"""Shared type aliases used across the deepretro package."""

from __future__ import annotations

from typing import TypeAlias

ScalarValue: TypeAlias = str | int | float | bool | None
ConditionValue: TypeAlias = str | int | float
JSONValue: TypeAlias = ScalarValue | list["JSONValue"] | dict[str, "JSONValue"]
