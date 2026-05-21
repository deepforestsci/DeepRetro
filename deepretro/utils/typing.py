"""Shared type aliases for the deepretro.utils package."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

RouteNode = Mapping[str, Any]
Step = dict[str, Any]
DependencyMap = dict[str, list[str]]
ParseOutput = dict[str, Any]
