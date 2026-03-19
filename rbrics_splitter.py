"""
r-BRICS Balanced Molecular Splitter

Splits large molecules into roughly equal fragments for retrosynthesis.

Logic:
  1. r-BRICS finds all candidate bonds (more than BRICS — handles fused rings, chains)
  2. Classify bonds: macrocycle-ring / small-ring / non-ring
  3. Pool = macrocycle bonds + non-ring bonds (small rings NEVER cut)
  4. Auto-determine fragment count from molecule size
  5. Try all combos, pick most balanced split

Usage:
    from rbrics_splitter import RBricsSplitter
    splitter = RBricsSplitter()
    result = splitter.split("SMILES_HERE")          # auto fragment count
    result = splitter.split("SMILES_HERE", n_fragments=4)  # force 4
"""

import sys, os, math
from itertools import combinations
from dataclasses import dataclass, field
from pathlib import Path

# r-BRICS import
RBRICS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "r-BRICS")
if os.path.exists(RBRICS_PATH):
    sys.path.insert(0, RBRICS_PATH)

from rdkit import Chem
from rdkit.Chem import FragmentOnBonds, GetMolFrags

try:
    from rBRICS_public import FindrBRICSBonds
    HAS_RBRICS = True
except ImportError:
    HAS_RBRICS = False
    from rdkit.Chem import BRICS


@dataclass
class SplitResult:
    smiles: str
    n_heavy_atoms: int
    fragments: list[str]
    fragment_sizes: list[int]
    bonds_cut: list[dict]
    balance_cv: float = 0.0
    n_fragments: int = 0
    auto_n: bool = False
    strategy: str = ""

    def __post_init__(self):
        self.n_fragments = len(self.fragments)

    def __repr__(self):
        lines = [
            f"SplitResult(strategy='{self.strategy}')",
            f"  {self.n_heavy_atoms} atoms -> {self.n_fragments} fragments "
            f"(balance CV={self.balance_cv:.3f})",
        ]
        total = sum(self.fragment_sizes) or 1
        for i, (size, frag) in enumerate(zip(self.fragment_sizes, self.fragments)):
            pct = size * 100 / total
            bar = "#" * max(1, int(pct / 2.5))
            lines.append(f"  [{i}] {size:3d} atoms ({pct:4.1f}%) {bar}")
        if self.bonds_cut:
            lines.append(f"  Bonds cut:")
            for b in self.bonds_cut:
                lines.append(f"    {b['symbols']} {b['type']} ({b['category']})")
        return "\n".join(lines)


class RBricsSplitter:
    """
    Fragment molecules using r-BRICS + balance optimization.
    Preserves small rings (pyrans, cyclopentanes, etc).
    Only cuts macrocycle backbone bonds and non-ring bonds.
    """

    def __init__(self, min_ring_size_to_cut: int = 10):
        """
        Args:
            min_ring_size_to_cut: Only cut bonds in rings with this many
                atoms or more. Default 10 = macrolactones, macrocycles.
                Smaller rings (5,6-membered pyrans etc) are never cut.
        """
        self.min_ring_size_to_cut = min_ring_size_to_cut
        if not HAS_RBRICS:
            print("WARNING: r-BRICS not found, falling back to standard BRICS")

    def _find_bonds(self, mol):
        """Find all r-BRICS bonds, return as list of dicts."""
        if HAS_RBRICS:
            raw_bonds = list(FindrBRICSBonds(mol))
        else:
            from rdkit.Chem import BRICS
            raw_bonds = list(BRICS.FindBRICSBonds(mol))

        results = []
        seen = set()
        for (a1, a2), (t1, t2) in raw_bonds:
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None:
                continue
            bid = bond.GetIdx()
            if bid in seen:
                continue
            seen.add(bid)

            results.append({
                "bond_idx": bid,
                "atoms": (a1, a2),
                "symbols": (
                    f"{mol.GetAtomWithIdx(a1).GetSymbol()}({a1})-"
                    f"{mol.GetAtomWithIdx(a2).GetSymbol()}({a2})"
                ),
                "type": f"L{t1}-L{t2}",
                "in_ring": bond.IsInRing(),
            })
        return results

    def _classify_bonds(self, mol, bonds):
        """
        Classify bonds into:
          - macrocycle: in a ring with >= min_ring_size_to_cut atoms,
                        but NOT also in a small ring (shared/fused bond)
          - small_ring: in a ring with < min_ring_size_to_cut atoms
          - nonring: not in any ring
        """
        ri = mol.GetRingInfo()
        atom_rings = list(ri.AtomRings())

        macro_bonds = []
        small_ring_bonds = []
        nonring_bonds = []

        for b in bonds:
            if not b["in_ring"]:
                b["category"] = "non-ring"
                nonring_bonds.append(b)
                continue

            a1, a2 = b["atoms"]
            in_macro = False
            in_small = False

            for ring_atoms in atom_rings:
                if a1 in ring_atoms and a2 in ring_atoms:
                    if len(ring_atoms) >= self.min_ring_size_to_cut:
                        in_macro = True
                    else:
                        in_small = True

            if in_small:
                b["category"] = "small-ring (protected)"
                small_ring_bonds.append(b)
            elif in_macro:
                b["category"] = "macrocycle"
                macro_bonds.append(b)
            else:
                b["category"] = "small-ring (protected)"
                small_ring_bonds.append(b)

        return macro_bonds, small_ring_bonds, nonring_bonds

    def _balance_score(self, sizes):
        """Coefficient of variation. 0 = perfect balance."""
        if len(sizes) <= 1:
            return float("inf")
        mean = sum(sizes) / len(sizes)
        if mean == 0:
            return float("inf")
        variance = sum((s - mean) ** 2 for s in sizes) / len(sizes)
        return math.sqrt(variance) / mean

    def _auto_fragment_count(self, n_atoms, target_size=15):
        """
        Decide how many fragments based on molecule size.

        target_size: ideal heavy atoms per fragment (default 15).

        Rules:
          - < 20 atoms: no split (1)
          - 20-35 atoms: 2 fragments
          - 35-50 atoms: 3 fragments
          - 50-70 atoms: 4 fragments
          - 70-100 atoms: 5 fragments
          - 100+: 6-7 fragments

        Or simply: round(n_atoms / target_size), clamped to [1, 7]
        """
        n = max(1, round(n_atoms / target_size))
        return max(1, min(n, 7))

    def split(
        self,
        smiles: str,
        n_fragments: int | None = None,
        target_fragment_size: int = 15,
        min_fragment_size: int = 4,
        max_combinations: int = 50000,
    ) -> SplitResult:
        """
        Split a molecule into balanced fragments.

        Args:
            smiles: input SMILES
            n_fragments: target fragment count (None = auto from molecule size)
            target_fragment_size: ideal atoms per fragment (used for auto n)
            min_fragment_size: reject splits producing fragments smaller than this
            max_combinations: safety cap on combinatorial search
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        canonical = Chem.MolToSmiles(mol)
        n_atoms = mol.GetNumHeavyAtoms()

        auto_n = n_fragments is None
        if auto_n:
            n_fragments = self._auto_fragment_count(n_atoms, target_fragment_size)

        if n_fragments <= 1 or n_atoms <= target_fragment_size:
            return SplitResult(
                smiles=canonical, n_heavy_atoms=n_atoms,
                fragments=[canonical], fragment_sizes=[n_atoms],
                bonds_cut=[], balance_cv=0.0, auto_n=auto_n,
                strategy=f"no_split_needed({n_atoms} atoms)",
            )

        all_bonds = self._find_bonds(mol)
        macro, small, nonring = self._classify_bonds(mol, all_bonds)

        pool = macro + nonring

        if not pool:
            return SplitResult(
                smiles=canonical, n_heavy_atoms=n_atoms,
                fragments=[canonical], fragment_sizes=[n_atoms],
                bonds_cut=[], balance_cv=0.0, auto_n=auto_n,
                strategy="no_cuttable_bonds",
            )

        n_cuts = n_fragments - 1

        if n_cuts > len(pool):
            n_cuts = len(pool)

        n_combos = math.comb(len(pool), n_cuts)

        if n_combos > max_combinations:
            priority_pool = macro + sorted(
                nonring,
                key=lambda b: 0 if "O" in b["symbols"] or "N" in b["symbols"] else 1,
            )
            pool = priority_pool[: n_cuts + 8]
            n_combos = math.comb(len(pool), n_cuts)

        best_cv = float("inf")
        best_result = None

        for combo in combinations(range(len(pool)), n_cuts):
            bond_indices = [pool[i]["bond_idx"] for i in combo]
            try:
                frag_mol = FragmentOnBonds(mol, bond_indices, addDummies=True)
                frag_mols = GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
                sizes = [fm.GetNumHeavyAtoms() for fm in frag_mols]

                if any(s < min_fragment_size for s in sizes):
                    continue

                cv = self._balance_score(sizes)

                if cv < best_cv:
                    best_cv = cv
                    smiles_list = [Chem.MolToSmiles(fm) for fm in frag_mols]
                    paired = sorted(zip(sizes, smiles_list), reverse=True)

                    best_result = SplitResult(
                        smiles=canonical,
                        n_heavy_atoms=n_atoms,
                        fragments=[p[1] for p in paired],
                        fragment_sizes=[p[0] for p in paired],
                        bonds_cut=[pool[i] for i in combo],
                        balance_cv=cv,
                        auto_n=auto_n,
                        strategy=(
                            f"rbrics_balanced(n={n_cuts},"
                            f"pool={len(pool)},combos={n_combos})"
                        ),
                    )
            except Exception:
                continue

        if best_result is None:
            return SplitResult(
                smiles=canonical, n_heavy_atoms=n_atoms,
                fragments=[canonical], fragment_sizes=[n_atoms],
                bonds_cut=[], balance_cv=0.0, auto_n=auto_n,
                strategy=f"no_valid_split(min_size={min_fragment_size})",
            )

        return best_result

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    def visualize_split(
        self,
        result: "SplitResult",
        save_path: str | None = None,
        show: bool = True,
        img_size: tuple[int, int] = (950, 500),
    ):
        """
        Draw the original molecule with cut-bonds highlighted in red,
        candidate-but-uncut bonds in blue, and the resulting fragments
        colour-coded underneath.

        Args:
            result:    A SplitResult returned by .split()
            save_path: If given, save the figure as PNG to this path.
            show:      Call plt.show()  (set False in scripts).
            img_size:  (w, h) in pixels for each RDKit sub-image.
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
        from rdkit.Chem import AllChem, Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        from io import BytesIO
        from PIL import Image

        # ---- palette for fragments ----
        FRAG_COLORS = [
            (0.22, 0.56, 0.87, 0.35),   # blue
            (0.90, 0.30, 0.24, 0.35),   # red
            (0.30, 0.75, 0.40, 0.35),   # green
            (0.80, 0.55, 0.10, 0.35),   # amber
            (0.58, 0.30, 0.75, 0.35),   # purple
            (0.10, 0.70, 0.70, 0.35),   # teal
            (0.85, 0.45, 0.65, 0.35),   # pink
        ]
        CUT_BOND_COLOR = (0.90, 0.15, 0.15)   # red
        CAND_BOND_COLOR = (0.40, 0.60, 0.85)  # light blue

        # ---- rebuild the mol and find bonds again ----
        mol = Chem.MolFromSmiles(result.smiles)
        if mol is None:
            raise ValueError(f"Cannot parse SMILES: {result.smiles}")
        AllChem.Compute2DCoords(mol)

        all_bonds = self._find_bonds(mol)
        macro, small, nonring = self._classify_bonds(mol, all_bonds)
        cuttable_pool = macro + nonring

        # bond-idx sets
        cut_bids = {b["bond_idx"] for b in result.bonds_cut}
        cand_bids = {b["bond_idx"] for b in cuttable_pool} - cut_bids

        # atoms adjacent to cut bonds (for highlight)
        cut_atoms = set()
        for b in result.bonds_cut:
            cut_atoms.update(b["atoms"])

        # ---- draw original molecule ----
        def _draw_mol(m, w, h, highlight_atoms=None, highlight_bonds=None,
                       atom_colors=None, bond_colors=None, legend=""):
            d = rdMolDraw2D.MolDraw2DCairo(w, h)
            opts = d.drawOptions()
            opts.useBWAtomPalette()
            opts.additionalAtomLabelPadding = 0.15
            opts.bondLineWidth = 2.5
            opts.setBackgroundColour((1, 1, 1, 1))

            d.DrawMolecule(
                m,
                highlightAtoms=highlight_atoms or [],
                highlightBonds=highlight_bonds or [],
                highlightAtomColors=atom_colors or {},
                highlightBondColors=bond_colors or {},
                legend=legend,
            )
            d.FinishDrawing()
            return Image.open(BytesIO(d.GetDrawingText()))

        # Colour maps for original mol
        bond_highlight = list(cut_bids | cand_bids)
        bond_color_map = {}
        for bid in cut_bids:
            bond_color_map[bid] = CUT_BOND_COLOR
        for bid in cand_bids:
            bond_color_map[bid] = CAND_BOND_COLOR

        atom_color_map = {}
        for a in cut_atoms:
            atom_color_map[a] = CUT_BOND_COLOR

        orig_img = _draw_mol(
            mol, img_size[0], img_size[1],
            highlight_atoms=list(cut_atoms),
            highlight_bonds=bond_highlight,
            atom_colors=atom_color_map,
            bond_colors=bond_color_map,
            legend="Original (red = cut, blue = candidate)",
        )

        # ---- draw fragments ----
        frag_imgs = []
        for i, smi in enumerate(result.fragments):
            fm = Chem.MolFromSmiles(smi)
            if fm is None:
                continue
            AllChem.Compute2DCoords(fm)
            c = FRAG_COLORS[i % len(FRAG_COLORS)]
            # highlight every heavy atom in that fragment's colour
            atoms = list(range(fm.GetNumAtoms()))
            ac = {a: c[:3] for a in atoms}
            sz = result.fragment_sizes[i]
            img = _draw_mol(
                fm, img_size[0] // 2, img_size[1] // 2,
                highlight_atoms=atoms, atom_colors=ac,
                legend=f"Frag {i}  ({sz} atoms)",
            )
            frag_imgs.append(img)

        # ---- assemble with matplotlib ----
        n_frags = len(frag_imgs)
        cols = min(n_frags, 4)
        rows = math.ceil(n_frags / cols) if n_frags else 1
        fig_h = 5 + 3 * rows

        fig, axes = plt.subplots(
            1 + rows, max(cols, 1),
            figsize=(4 * max(cols, 1), fig_h),
            gridspec_kw={"height_ratios": [2] + [1] * rows},
        )
        fig.suptitle(
            f"r-BRICS Split -- {result.n_heavy_atoms} atoms -> "
            f"{result.n_fragments} fragments   (CV={result.balance_cv:.3f})",
            fontsize=13, fontweight="bold", y=0.98,
        )

        # make axes 2-D always
        if axes.ndim == 1:
            axes = axes.reshape(-1, 1) if cols == 1 else axes.reshape(1 + rows, -1)

        # top row: original molecule spans all columns
        for c_idx in range(cols):
            axes[0, c_idx].axis("off")
        # place original image centered in top row
        ax_orig = fig.add_subplot(1 + rows, 1, 1)
        ax_orig.imshow(orig_img)
        ax_orig.set_title(
            "Red bonds = cut  |  Blue bonds = other candidates",
            fontsize=10, color="grey",
        )
        ax_orig.axis("off")

        # bottom rows: fragments
        for i, img in enumerate(frag_imgs):
            r = 1 + i // cols
            c_idx = i % cols
            axes[r, c_idx].imshow(img)
            axes[r, c_idx].axis("off")
        # hide unused axes
        for i in range(n_frags, rows * cols):
            r = 1 + i // cols
            c_idx = i % cols
            axes[r, c_idx].axis("off")

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight",
                        facecolor="white")
            print(f"Saved visualisation -> {save_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    def analyze(self, smiles: str):
        """Show all cuttable bonds and ring classification."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        n_atoms = mol.GetNumHeavyAtoms()
        ri = mol.GetRingInfo()

        print(f"Heavy atoms: {n_atoms}")
        print(f"Rings:")
        for i, r in enumerate(ri.AtomRings()):
            label = "MACROCYCLE" if len(r) >= self.min_ring_size_to_cut else "small"
            print(f"  Ring {i}: {len(r)} atoms ({label})")

        all_bonds = self._find_bonds(mol)
        macro, small, nonring = self._classify_bonds(mol, all_bonds)

        auto_n = self._auto_fragment_count(n_atoms)

        print(f"\nBond pool:")
        print(f"  Macrocycle bonds (cuttable): {len(macro)}")
        for b in macro:
            print(f"    {b['symbols']} {b['type']}")
        print(f"  Non-ring bonds (cuttable):   {len(nonring)}")
        for b in nonring:
            print(f"    {b['symbols']} {b['type']}")
        print(f"  Small-ring bonds (PROTECTED): {len(small)}")
        print(f"\n  Total cuttable: {len(macro) + len(nonring)}")
        print(f"  Auto fragment count: {auto_n} (target ~15 atoms each)")

        return macro, small, nonring


# =============================================================================
if __name__ == "__main__":
    splitter = RBricsSplitter()

    molecules = {
        "bryostatin_1": (
            "[H][C@]12C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@]3([H])"
            r"C\C(C[C@@]([H])(O3)\C=C\C(C)(C)[C@]3(O)O[C@@]([H])"
            r"(C\C(=C/C(=O)OC)[C@@H]3OC(=O)\C=C\C=C\CCC)"
            r"C[C@@]([H])(OC(=O)C[C@H](O)C1)[C@@H](C)O)=C\C(=O)OC)O2"
        ),
        "sorafenib": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",
        "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "erythromycin": (
            "CC[C@@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)"
            "[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)"
            "C[C@@H](N(C)C)[C@@H]2O)[C@](C)(O)C[C@@H](C)C(=O)"
            "[C@H](C)[C@@H](O)[C@]1(C)O"
        ),
        "small_drug": "c1ccc(NC(=O)c2ccccc2)cc1",
    }

    out_dir = Path(__file__).parent / "split_visualisations"
    out_dir.mkdir(exist_ok=True)

    for name, smi in molecules.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"\nFAILED TO PARSE: {name}")
            continue

        print(f"\n{'='*65}")
        print(f"{name.upper()} ({mol.GetNumHeavyAtoms()} heavy atoms)")
        print(f"{'='*65}")

        splitter.analyze(smi)

        print(f"\n--- AUTO SPLIT ---")
        result = splitter.split(smi)
        print(result)

        # Save visualisation of the auto split
        splitter.visualize_split(
            result,
            save_path=str(out_dir / f"{name}_auto.png"),
            show=False,
        )

        if mol.GetNumHeavyAtoms() > 40:
            print(f"\n--- FORCED 4 FRAGMENTS ---")
            result4 = splitter.split(smi, n_fragments=4)
            print(result4)
            splitter.visualize_split(
                result4,
                save_path=str(out_dir / f"{name}_4frag.png"),
                show=False,
            )
