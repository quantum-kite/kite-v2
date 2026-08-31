"""Native, pybinding-free lattice templates, mirroring `pybinding.repository`.

This module ports a handful of standard tight-binding models from
`pybinding.repository` (graphene, phosphorene, group-6 TMDs) to KITE's native
`kite.lattice.Lattice` construction, so example scripts can eventually stop
importing `pybinding` just to get a lattice definition. Module/function/constant
names deliberately mirror pybinding's own (`graphene.a`, `graphene.monolayer`,
`phosphorene.monolayer_4band`, `group6_tmd.monolayer_3band`, ...) per explicit
direction from the user -- this is a *port*, not a re-export: no `pybinding`
import anywhere in this file.

Because this lives in a single flat module (not a package), pybinding's
`graphene`/`phosphorene`/`group6_tmd` *modules* are represented here as plain
namespace classes (never instantiated) holding the same constants and
functions as `staticmethod`s -- e.g. `kite.repository.graphene.monolayer()`
and `kite.repository.graphene.a` play the same role as
`pybinding.repository.graphene.monolayer()` and `pybinding.repository.graphene.a`.

Units: all lengths are in nm and all hopping/onsite energies are in eV,
matching the convention already used by pybinding's own repository modules
and by this repo's existing pybinding-free examples (`examples/dos_graphene.py`,
`examples/dccond_phosphorene.py`).

Multi-orbital sites: `kite.lattice.Lattice.add_one_sublattice` only accepts a
single scalar onsite energy, unlike pybinding's `Lattice.add_sublattices`
(which accepts a small onsite matrix for a multi-orbital site, e.g. TMDs'
3-orbital metal site). The TMD model below therefore follows the pattern
already used natively in this repo (`examples/pybinding/arpes_tmd.py`'s
`tmd_monolayer`, which itself documents this as "the recommended construction
for multi-orbital sites in KITE"): one named sublattice per orbital, connected
by explicit scalar hoppings, rather than bundling orbitals into one
sublattice with a dense hopping matrix.
"""

import math
import re

import numpy as np

from kite.lattice import Lattice


class graphene:
    """Monolayer graphene (honeycomb lattice), nearest + optional next-nearest
    neighbor hopping.

    Constants match `pybinding.repository.graphene.constants` exactly (verified
    numerically against the installed pybinding package -- see the
    verification script referenced in the implementing session's report).

    Deferred (not ported in this round -- see module docstring and the
    session report for reasoning): `monolayer_alt` (alternative lattice
    vectors, same physics), `monolayer_4atom` (needs `Lattice.add_aliases`,
    which `kite.lattice.Lattice` does not have), `bilayer` (same alias
    dependency), and the 3rd-nearest-neighbor `t_nnn` term.
    """

    a = 0.24595       #: [nm] unit cell length
    a_cc = 0.142       #: [nm] carbon-carbon distance
    t = -2.8           #: [eV] nearest-neighbor hopping
    t_nn = 0.1         #: [eV] next-nearest-neighbor hopping

    @staticmethod
    def monolayer(nearest_neighbors=1, onsite=(0.0, 0.0), t=None, t_nn=None):
        """Monolayer graphene lattice with 1 or 2 nearest-neighbor shells.

        Parameters
        ----------
        nearest_neighbors : int
            1: nearest-neighbor only. 2: also add next-nearest-neighbor
            hopping `t_nn` (onsite energies are shifted by `+3*t_nn` to keep
            the Dirac point at the same energy as the nearest-neighbor-only
            model, exactly matching pybinding's own compensation). No other
            value is supported (`t_nnn`/3rd-neighbor is out of scope for now
            -- raises `ValueError` rather than silently ignoring the request).
        onsite : tuple of float
            Onsite energy for sublattices A and B (before the `nearest_neighbors
            == 2` compensation shift, if any).
        t, t_nn : float, optional
            Override the default hopping constants (`graphene.t`, `graphene.t_nn`).

        Returns
        -------
        kite.lattice.Lattice
        """
        if nearest_neighbors not in (1, 2):
            raise ValueError(
                "graphene.monolayer only supports nearest_neighbors in (1, 2) "
                "in this version (3rd-neighbor t_nnn is deferred)"
            )

        t_val = graphene.t if t is None else t
        t_nn_val = graphene.t_nn if t_nn is None else t_nn
        a, a_cc = graphene.a, graphene.a_cc

        lat = Lattice(a1=[a, 0.0], a2=[a / 2, a / 2 * math.sqrt(3)])

        onsite_offset = 0.0 if nearest_neighbors < 2 else 3 * t_nn_val
        lat.add_sublattices(
            ("A", [0.0, -a_cc / 2], onsite[0] + onsite_offset),
            ("B", [0.0, a_cc / 2], onsite[1] + onsite_offset),
        )

        lat.register_hopping_energies({"t": t_val, "t_nn": t_nn_val})

        lat.add_hoppings(
            ([0, 0], "A", "B", "t"),
            ([1, -1], "A", "B", "t"),
            ([0, -1], "A", "B", "t"),
        )

        if nearest_neighbors >= 2:
            lat.add_hoppings(
                ([0, -1], "A", "A", "t_nn"),
                ([0, -1], "B", "B", "t_nn"),
                ([1, -1], "A", "A", "t_nn"),
                ([1, -1], "B", "B", "t_nn"),
                ([1, 0], "A", "A", "t_nn"),
                ([1, 0], "B", "B", "t_nn"),
            )

        return lat


class phosphorene:
    """Monolayer black phosphorus (phosphorene), 4-band nearest-neighbor model.

    Constants match `pybinding.repository.phosphorene` exactly (verified
    numerically). This is the same model already used natively in
    `examples/dccond_phosphorene.py`'s `phosphorene_lattice`; this port is the
    library-level equivalent (mirroring pybinding's module naming) rather than
    an example-embedded copy.

    Deferred: pybinding's `phosphorene` module only offers this single
    (`monolayer_4band`) model -- nothing else to port here.
    """

    a = 0.222                        #: [nm]
    ax = 0.438                        #: [nm]
    ay = 0.332                        #: [nm]
    theta = 96.79 * (math.pi / 180)   #: [rad]
    phi = 103.69 * (math.pi / 180)    #: [rad]

    @staticmethod
    def monolayer_4band(num_hoppings=5):
        """Monolayer phosphorene lattice using the four-band model.

        Parameters
        ----------
        num_hoppings : int
            Number of hopping shells to include, from `t1`/`t2` (minimum,
            `num_hoppings >= 2`) up to `t5` (`num_hoppings == 5`, the
            default -- matches pybinding's own default).

        Returns
        -------
        kite.lattice.Lattice
        """
        if num_hoppings < 2:
            raise ValueError("phosphorene.monolayer_4band: t1 and t2 must be included (num_hoppings >= 2)")
        if num_hoppings > 5:
            raise ValueError("phosphorene.monolayer_4band: t5 is the last shell (num_hoppings <= 5)")

        a, ax, ay = phosphorene.a, phosphorene.ax, phosphorene.ay
        theta, phi = phosphorene.theta, phosphorene.phi

        h = a * math.sin(phi - math.pi / 2)
        s = 0.5 * ax - a * math.cos(theta / 2)

        lat = Lattice(a1=[ax, 0.0], a2=[0.0, ay])

        lat.add_sublattices(
            ("A", [-s / 2, -ay / 2, h], 0.0),
            ("B", [s / 2, -ay / 2, 0.0], 0.0),
            ("C", [-s / 2 + ax / 2, 0.0, 0.0], 0.0),
            ("D", [s / 2 + ax / 2, 0.0, h], 0.0),
        )

        lat.register_hopping_energies({
            "t1": -1.22,
            "t2": 3.665,
            "t3": -0.205,
            "t4": -0.105,
            "t5": -0.055,
        })

        if num_hoppings >= 2:
            lat.add_hoppings(
                ([-1, 0], "A", "D", "t1"),
                ([-1, -1], "A", "D", "t1"),
                ([0, 0], "B", "C", "t1"),
                ([0, -1], "B", "C", "t1"),
            )
            lat.add_hoppings(
                ([0, 0], "A", "B", "t2"),
                ([0, 0], "C", "D", "t2"),
            )
        if num_hoppings >= 3:
            lat.add_hoppings(
                ([0, 0], "A", "D", "t3"),
                ([0, -1], "A", "D", "t3"),
                ([1, 1], "C", "B", "t3"),
                ([1, 0], "C", "B", "t3"),
            )
        if num_hoppings >= 4:
            lat.add_hoppings(
                ([0, 0], "A", "C", "t4"),
                ([0, -1], "A", "C", "t4"),
                ([-1, 0], "A", "C", "t4"),
                ([-1, -1], "A", "C", "t4"),
                ([0, 0], "B", "D", "t4"),
                ([0, -1], "B", "D", "t4"),
                ([-1, 0], "B", "D", "t4"),
                ([-1, -1], "B", "D", "t4"),
            )
        if num_hoppings >= 5:
            lat.add_hoppings(
                ([-1, 0], "A", "B", "t5"),
                ([-1, 0], "C", "D", "t5"),
            )

        return lat


class group6_tmd:
    """Monolayer group-6 transition metal dichalcogenides, nearest-neighbor
    3-band (d_z2, d_xy, d_x2-y2) model of Liu et al., Phys. Rev. B 88, 085433
    (2013).

    Parameter table matches `pybinding.repository.group6_tmd._default_3band_params`
    exactly (verified numerically) -- same 6 materials (MoS2, WS2, MoSe2,
    WSe2, MoTe2, WTe2), same units (lattice constant in nm).

    Multi-orbital construction: unlike pybinding (one sublattice with a 3x3
    onsite/hopping matrix), this uses one sublattice per orbital
    (`{metal}_0`, `{metal}_1`, `{metal}_2` for d_z2, d_xy, d_x2-y2
    respectively) with scalar hoppings -- see module docstring. This matches
    the pattern already used natively in `examples/pybinding/arpes_tmd.py`'s
    `tmd_monolayer` (which itself hard-codes MoS2 only); this port
    generalizes that pattern to all 6 species pybinding offers.
    """

    #: name -> [a (nm), eps1, eps2, t0, t1, t2, t11, t12, t22]
    #: from https://doi.org/10.1103/PhysRevB.88.085433
    _default_3band_params = {
        "MoS2":  [0.3190, 1.046, 2.104, -0.184, 0.401, 0.507, 0.218, 0.338, 0.057],
        "WS2":   [0.3191, 1.130, 2.275, -0.206, 0.567, 0.536, 0.286, 0.384, -0.061],
        "MoSe2": [0.3326, 0.919, 2.065, -0.188, 0.317, 0.456, 0.211, 0.290, 0.130],
        "WSe2":  [0.3325, 0.943, 2.179, -0.207, 0.457, 0.486, 0.263, 0.329, 0.034],
        "MoTe2": [0.3557, 0.605, 1.972, -0.169, 0.228, 0.390, 0.207, 0.239, 0.252],
        "WTe2":  [0.3560, 0.606, 2.102, -0.175, 0.342, 0.410, 0.233, 0.270, 0.190],
    }

    @staticmethod
    def _hopping_matrices(t0, t1, t2, t11, t12, t22):
        """The 3 nearest-neighbor 3x3 hopping matrices (one per lattice
        direction [1,0], [0,-1], [1,-1]), in the (d_z2, d_xy, d_x2-y2) basis.
        Identical formulas to pybinding's `group6_tmd.monolayer_3band` (and
        to `examples/pybinding/arpes_tmd.py`'s `_hopping_matrices`)."""
        rt3 = math.sqrt(3)
        h1 = [[t0, -t1, t2],
              [t1, t11, -t12],
              [t2, t12, t22]]
        h2 = [[t0, 0.5 * t1 + rt3 / 2 * t2, rt3 / 2 * t1 - 0.5 * t2],
              [-0.5 * t1 + rt3 / 2 * t2, 0.25 * t11 + 0.75 * t22, rt3 / 4 * (t11 - t22) - t12],
              [-rt3 / 2 * t1 - 0.5 * t2, rt3 / 4 * (t11 - t22) + t12, 0.75 * t11 + 0.25 * t22]]
        h3 = [[t0, -0.5 * t1 - rt3 / 2 * t2, rt3 / 2 * t1 - 0.5 * t2],
              [0.5 * t1 - rt3 / 2 * t2, 0.25 * t11 + 0.75 * t22, rt3 / 4 * (t22 - t11) + t12],
              [-rt3 / 2 * t1 - 0.5 * t2, rt3 / 4 * (t22 - t11) - t12, 0.75 * t11 + 0.25 * t22]]
        return [h1, h2, h3]

    @staticmethod
    def monolayer_3band(name, override_params=None):
        """Monolayer of a group-6 TMD, nearest-neighbor 3-band model.

        Parameters
        ----------
        name : str
            One of "MoS2", "WS2", "MoSe2", "WSe2", "MoTe2", "WTe2".
        override_params : dict, optional
            Replace or add entries in the parameter table, in the same
            `[a, eps1, eps2, t0, t1, t2, t11, t12, t22]` format.

        Returns
        -------
        kite.lattice.Lattice
        """
        params = dict(group6_tmd._default_3band_params)
        if override_params:
            params.update(override_params)

        if name not in params:
            raise ValueError(
                f"Unknown TMD '{name}'; available: {sorted(params)} "
                "(or supply it via override_params)"
            )

        a, eps1, eps2, t0, t1, t2, t11, t12, t22 = params[name]
        rt3 = math.sqrt(3)

        lat = Lattice(a1=[a, 0.0], a2=[0.5 * a, rt3 / 2 * a])

        metal_match = re.findall("[A-Z][a-z]*", name)
        if not metal_match:
            raise ValueError(f"Could not parse a metal element name out of '{name}'")
        metal = metal_match[0]
        subs = [f"{metal}_0", f"{metal}_1", f"{metal}_2"]  # d_z2, d_xy, d_x2-y2

        lat.add_sublattices(
            (subs[0], [0.0, 0.0], eps1),
            (subs[1], [0.0, 0.0], eps2),
            (subs[2], [0.0, 0.0], eps2),
        )

        directions = [[1, 0], [0, -1], [1, -1]]
        matrices = group6_tmd._hopping_matrices(t0, t1, t2, t11, t12, t22)
        for direction, matrix in zip(directions, matrices):
            for i in range(3):
                for j in range(3):
                    val = matrix[i][j]
                    if abs(val) > 1e-8:
                        lat.add_hoppings((direction, subs[i], subs[j], val))

        return lat
