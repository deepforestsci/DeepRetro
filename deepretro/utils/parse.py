"""Utilities for converting retrosynthesis route trees into viewer steps.

The module keeps the historical function API (``parse_step``,
``fix_dependencies``, and ``format_output``) while exposing
``RetrosynthesisRouteParser`` for callers that need dependency injection or a
clearer unit to test.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Iterable, Mapping, MutableMapping
from typing import Any, Optional

from deepretro.utils.parse_metrics import calc_scalability_index
from deepretro.utils.utils_molecule import calc_chemical_formula, calc_mol_wt
from deepretro.utils.variables import BASIC_MOLECULES

RouteNode = Mapping[str, Any]
Step = dict[str, Any]
DependencyMap = dict[str, list[str]]
ParseOutput = dict[str, Any]


class RetrosynthesisRouteParser:
    """
    Convert recursive retrosynthesis route trees into viewer-ready dictionaries.

    Examples
    --------
    >>> parser = RetrosynthesisRouteParser(
    ...     basic_molecules=set(),
    ...     chemical_formula_calculator=lambda smiles: smiles,
    ...     mass_calculator=lambda smiles: 1.0,
    ...     scalability_calculator=lambda reactant, product: "N/A",
    ... )
    >>> output = parser.format_output({
    ...     "smiles": "CCO",
    ...     "children": [{
    ...         "children": [{"smiles": "CC"}],
    ...     }],
    ... })
    >>> output["steps"][0]["products"][0]["smiles"]
    'CCO'
    """

    def __init__(
        self,
        basic_molecules: Optional[Collection[str]] = None,
        chemical_formula_calculator: Optional[Callable[[str], str]] = None,
        mass_calculator: Optional[Callable[[str], float]] = None,
        scalability_calculator: Optional[Callable[[str, str], str]] = None,
    ) -> None:
        """
        Parameters
        ----------
        basic_molecules : Collection[str], optional
            Molecules that should be emitted as reagents instead of reactants.
            Defaults to ``deepretro.utils.variables.BASIC_MOLECULES``.
        chemical_formula_calculator : Callable[[str], str], optional
            Function used to calculate formula metadata for each molecule.
        mass_calculator : Callable[[str], float], optional
            Function used to calculate molecular mass metadata.
        scalability_calculator : Callable[[str, str], str], optional
            Function used to score a precursor/product pair for scalability.
        """
        self.basic_molecules = (
            BASIC_MOLECULES if basic_molecules is None else basic_molecules
        )
        self.chemical_formula_calculator = (
            calc_chemical_formula
            if chemical_formula_calculator is None
            else chemical_formula_calculator
        )
        self.mass_calculator = (
            calc_mol_wt if mass_calculator is None else mass_calculator
        )
        self.scalability_calculator = (
            calc_scalability_index
            if scalability_calculator is None
            else scalability_calculator
        )

    def parse_step(
        self,
        data: RouteNode,
        include_metadata: Optional[bool] = None,
        step_list: Optional[list[Step]] = None,
        dependency_list: Optional[DependencyMap] = None,
        parent_id: Optional[int] = None,
    ) -> ParseOutput:
        """
        Parse a route subtree into dependency and step collections.

        Parameters
        ----------
        data : Mapping[str, Any]
            Route node containing a product ``smiles`` value and optional
            children.
        include_metadata : bool, optional
            Compatibility argument retained for the historical function API.
            Metadata is always included by the current output schema.
        step_list : list[dict[str, Any]], optional
            Existing step accumulator used by recursive callers.
        dependency_list : dict[str, list[str]], optional
            Existing dependency accumulator used by recursive callers.
        parent_id : int, optional
            One-based step id of the parent reaction.

        Returns
        -------
        dict[str, Any]
            Dictionary with ``dependencies`` and ``steps`` keys.

        Examples
        --------
        >>> parser = RetrosynthesisRouteParser(
        ...     basic_molecules={"O"},
        ...     chemical_formula_calculator=lambda smiles: smiles,
        ...     mass_calculator=lambda smiles: 1.0,
        ...     scalability_calculator=lambda reactant, product: "N/A",
        ... )
        >>> output = parser.parse_step(
        ...     {
        ...         "smiles": "CCO",
        ...         "children": [{"children": [{"smiles": "CC"}, {"smiles": "O"}]}],
        ...     }
        ... )
        >>> output["dependencies"]
        {'1': []}
        >>> [molecule["smiles"] for molecule in output["steps"][0]["reactants"]]
        ['CC']
        >>> [molecule["smiles"] for molecule in output["steps"][0]["reagents"]]
        ['O']
        """
        steps = [] if step_list is None else step_list
        dependencies = {} if dependency_list is None else dependency_list
        self._parse_node(data, steps, dependencies, parent_id)
        return {"dependencies": dependencies, "steps": steps}

    def format_output(self, data: RouteNode) -> ParseOutput:
        """
        Parse a route tree and rebuild dependencies from product/reactant links.

        Parameters
        ----------
        data : Mapping[str, Any]
            Root retrosynthesis route tree.

        Returns
        -------
        dict[str, Any]
            Viewer-ready route data with ``steps`` and ``dependencies``.
        """
        output_data = self.parse_step(data)
        output_data["dependencies"] = self.fix_dependencies(output_data["steps"])
        return output_data

    def fix_dependencies(self, step_list: Iterable[Step]) -> DependencyMap:
        """
        Rebuild step dependencies from matching product and reactant SMILES.

        Parameters
        ----------
        step_list : Iterable[dict[str, Any]]
            Parsed reaction steps. Each step should contain one product and
            zero or more reactants.

        Returns
        -------
        dict[str, list[str]]
            Mapping from each step id to the step ids that produce its
            reactants.
        """
        steps = list(step_list)
        product_to_step = {
            step["products"][0]["smiles"]: step["step"]
            for step in steps
            if step.get("products")
        }

        return {
            step["step"]: [
                product_to_step[reactant["smiles"]]
                for reactant in step.get("reactants", [])
                if reactant.get("smiles") in product_to_step
            ]
            for step in steps
        }

    def _parse_node(
        self,
        data: RouteNode,
        steps: list[Step],
        dependencies: MutableMapping[str, list[str]],
        parent_id: Optional[int],
    ) -> None:
        """
        Parse a single route node into steps and dependency links.

        Parameters
        ----------
        data : Mapping[str, Any]
            Route node to parse. Nodes may represent a product-bearing molecule
            with optional precursor children, or a leaf molecule without
            ``children``.
        steps : list[dict[str, Any]]
            Mutable step accumulator populated in traversal order.
        dependencies : MutableMapping[str, list[str]]
            Mutable dependency accumulator keyed by step id.
        parent_id : int, optional
            One-based step id of the parent reaction consuming this node as a
            precursor. ``None`` indicates the root node.
        """
        step = self._create_step(data, len(steps) + 1)
        self._attach_to_parent(data, steps, parent_id)

        if step is None:
            if parent_id is not None:
                dependencies.setdefault(str(parent_id), [])
            return

        steps.append(step)
        if parent_id is not None:
            dependencies.setdefault(str(parent_id), []).append(step["step"])

        for child in self._iter_reaction_children(data):
            self._parse_node(child, steps, dependencies, int(step["step"]))

    def _create_step(self, data: RouteNode, step_id: int) -> Optional[Step]:
        """Create a reaction step for route nodes that produce a molecule."""
        if "children" not in data:
            return None

        product_smiles = self._smiles(data)
        return {
            "step": str(step_id),
            "reactants": [],
            "reagents": [],
            "products": [self._molecule_entry(product_smiles, "product_metadata")],
            "conditions": [],
            "reactionmetrics": [
                {
                    "scalabilityindex": "",
                    "closestliterature": "",
                }
            ],
        }

    def _attach_to_parent(
        self,
        data: RouteNode,
        steps: list[Step],
        parent_id: Optional[int],
    ) -> None:
        """
        Attach a precursor molecule node to its parent reaction step.

        Parameters
        ----------
        data : Mapping[str, Any]
            Molecule or reaction-wrapper node being attached as an input to the
            parent reaction.
        steps : list[dict[str, Any]]
            Parsed steps accumulated so far. The parent step is updated in
            place when this node is a precursor molecule.
        parent_id : int, optional
            One-based identifier of the parent step that consumes this
            precursor. ``None`` means there is no parent step to update.

        Returns
        -------
        None
            The method mutates ``steps`` in place by appending a reactant or
            reagent entry and updating the parent step's scalability metric.
        """
        if parent_id is None or data.get("is_reaction", False):
            return

        smiles = self._smiles(data)
        parent_step = steps[parent_id - 1]
        molecule_type = "reagents" if smiles in self.basic_molecules else "reactants"
        metadata_key = f"{molecule_type[:-1]}_metadata"
        parent_step[molecule_type].append(self._molecule_entry(smiles, metadata_key))
        parent_step["reactionmetrics"][0]["scalabilityindex"] = (
            self.scalability_calculator(smiles, parent_step["products"][0]["smiles"])
        )

    def _molecule_entry(self, smiles: str, metadata_key: str) -> dict[str, Any]:
        """
        Build a molecule entry with the metadata key required by the schema.

        Parameters
        ----------
        smiles : str
            SMILES string for the molecule being emitted into the viewer
            schema.
        metadata_key : str
            Metadata field name expected by the caller, such as
            ``"product_metadata"`` or ``"reactant_metadata"``.

        Returns
        -------
        dict[str, Any]
            Molecule payload containing the original SMILES string and
            calculated formula and mass metadata under ``metadata_key``.
        """
        return {
            "smiles": smiles,
            metadata_key: {
                "name": "",
                "chemical_formula": self.chemical_formula_calculator(smiles),
                "mass": self.mass_calculator(smiles),
            },
        }

    @staticmethod
    def _iter_reaction_children(data: RouteNode) -> Iterable[RouteNode]:
        """Yield precursor molecule nodes from route reaction wrappers."""
        for child in data.get("children", []):
            yield from child.get("children", [child])

    @staticmethod
    def _smiles(data: RouteNode) -> str:
        """Return the required SMILES value from a route node."""
        try:
            return str(data["smiles"])
        except KeyError as exc:
            raise ValueError("Route node is missing required 'smiles' value.") from exc


def parse_step(
    data: RouteNode,
    include_metadata: Optional[bool] = None,
    step_list: Optional[list[Step]] = None,
    dependency_list: Optional[DependencyMap] = None,
    parent_id: Optional[int] = None,
) -> ParseOutput:
    """
    Parse a retrosynthesis route tree into dependencies and steps.
    This compatibility wrapper delegates to ``RetrosynthesisRouteParser``.
    
    Parameters
    ----------
    data : Mapping[str, Any]
        Route node containing a product ``smiles`` value and optional
        children.
    include_metadata : bool, optional
        Compatibility argument retained for the historical function API.
        Metadata is always included by the current output schema.
    step_list : list[dict[str, Any]], optional
        Existing step accumulator used by recursive callers.
    dependency_list : dict[str, list[str]], optional
        Existing dependency accumulator used by recursive callers.
    parent_id : int, optional
        One-based step id of the parent reaction.

    Returns
    -------
    dict[str, Any]
        Dictionary with ``dependencies`` and ``steps`` keys.

    Examples
    --------
    >>> output = parse_step(
    ...     {
    ...         "smiles": "CCO",
    ...         "children": [{"children": [{"smiles": "CC"}, {"smiles": "O"}]}],
    ...     }
    ... )
    >>> output["dependencies"]
    {'1': []}
    >>> output["steps"][0]["products"][0]["smiles"]
    'CCO'
    """
    parser = RetrosynthesisRouteParser()
    return parser.parse_step(
        data,
        include_metadata=include_metadata,
        step_list=step_list,
        dependency_list=dependency_list,
        parent_id=parent_id,
    )


def fix_dependencies(
    dependencies: Mapping[str, list[str]],
    step_list: Iterable[Step],
) -> DependencyMap:
    """
    Rebuild dependencies from parsed route steps.

    Parameters
    ----------
    dependencies : Mapping[str, list[str]]
        Compatibility argument retained from the historical API. The mapping is
        ignored because dependencies are derived from ``step_list``.
    step_list : Iterable[dict[str, Any]]
        Parsed route steps.

    Returns
    -------
    dict[str, list[str]]
        Rebuilt dependency mapping.
    """
    return RetrosynthesisRouteParser().fix_dependencies(step_list)


def format_output(data: RouteNode) -> ParseOutput:
    """
    Format retrosynthesis route data for visualization.

    Parameters
    ----------
    data : Mapping[str, Any]
        Root route tree returned by the retrosynthesis pipeline.

    Returns
    -------
    dict[str, Any]
        Viewer-ready route data with ``steps`` and ``dependencies``.
    """
    return RetrosynthesisRouteParser().format_output(data)
