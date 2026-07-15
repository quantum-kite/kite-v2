import numpy as np
from scipy.spatial import Voronoi

class Sublattice:
    def __init__(self, position, energy, alias_id):
        self.position = position
        self.energy = energy
        self.alias_id = alias_id

class HoppingFamily:
    def __init__(
        self,
        relative_index,
        from_sub: str,
        to_sub: str,
        energy,
        from_id: int,
        to_id: int,
    ):
        self.relative_index = relative_index
        self.from_sub = from_sub
        self.to_sub = to_sub
        self.energy = energy
        self.terms = [HoppingTerms(relative_index, from_id, to_id)]

class HoppingTerms:
    def __init__(self, relative_index: list, from_id: int, to_id: int):
        self.relative_index = relative_index
        self.from_id = from_id
        self.to_id = to_id

# Class that replaces the pybinding for lattice generation
class Lattice:
    def __init__(self, a1, a2=None, a3=None):
        self._vectors = [v for v in (a1, a2, a3) if v is not None]
        self.sublattices = {}
        self.hoppings = {}
        self._hopping_energies_map = {}
        self.nsub = 0
        self.nhop = 0

    @property
    def vectors(self) -> list:
        """Primitive lattice vectors"""
        return self._vectors

    def add_one_hopping(
        self, relative_index: list[int], from_sub: str, to_sub: str, hop_name_or_energy
    ) -> None:
        """Does not add the complex conjugate hopping"""
        name = f"Hop{self.nhop}"

        if type(hop_name_or_energy) == str:
            if hop_name_or_energy in self._hopping_energies_map:
                hop_val = self._hopping_energies_map[hop_name_or_energy]
            else:
                raise TypeError(f"The hopping label has not been defined")
        else:
            hop_val = complex(hop_name_or_energy)

        self.hoppings[name] = HoppingFamily(
            np.asarray(relative_index),
            from_sub,
            to_sub,
            np.array([[hop_val]]),
            self.sublattices[from_sub].alias_id,
            self.sublattices[to_sub].alias_id,
        )
        self.nhop += 1

    def add_hoppings(self, *hoppings) -> None:
        """Does not add the complex conjugate hopping"""
        for hop in hoppings:
            self.add_one_hopping(*hop)

    def add_one_sublattice(
        self, name: str, position: list[float], onsite_energy: float = 0.0
    ) -> None:
        if name in self.sublattices:
            raise ValueError(f"Sublattice '{name}' already exists")

        self.sublattices[name] = Sublattice(
            position, np.array([[onsite_energy]]), self.nsub
        )
        self.nsub += 1

    def add_sublattices(self, *sublattices) -> None:
        for sub in sublattices:
            self.add_one_sublattice(*sub)

    def register_hopping_energies(self, mapping: dict[str, float]) -> None:
        """Register a mapping of user-friendly names to hopping energies"""
        self._hopping_energies_map = mapping

    def reciprocal_vectors(self) -> list:
        """Reciprocal-space primitive vectors, dimension-agnostic.

        Let A be the m x D matrix whose rows are the m in {1, 2, 3} real-space
        primitive vectors (D >= m; D > m happens when e.g. a 2D lattice is defined
        using 3-component position vectors). Requiring b_i to lie in the row space
        of A (the only physically sensible choice when D > m) together with the
        standard duality a_i . b_j = 2*pi*delta_ij gives

            B = 2*pi * (A @ A.T)^-1 @ A

        which reduces to the familiar 2*pi * A^-T when m == D. This is the exact
        relation used by KITE's own k-vector fractional/Cartesian conversion
        (see tools/Src/Spectral/arpes.cpp), so downstream k-vectors line up with
        the rest of the codebase's convention.

        Returns a list of m numpy arrays (rows of B), matching the list-of-arrays
        style of the `vectors` property.
        """
        A = np.array(self._vectors, dtype=float)
        m = A.shape[0]
        B = 2 * np.pi * np.linalg.inv(A @ A.T) @ A
        return [B[i] for i in range(m)]

    def brillouin_zone(self) -> np.ndarray:
        """Brillouin zone: the Voronoi cell of the reciprocal lattice centered at k=0.

        Only m = 1 and m = 2 (1D and 2D lattices) are supported for now:
          - m == 1: returns the two boundary points [-b1/2, +b1/2] as a (2, 1) array.
          - m == 2: builds a small shell of candidate reciprocal lattice points
            (integer combinations n_i in {-2,...,2} of the reciprocal vectors),
            runs scipy.spatial.Voronoi on them, and extracts the (finite, convex)
            region belonging to the point at the origin. Vertices are returned
            sorted counter-clockwise by polar angle (arctan2), matching pybinding's
            convention, as a numpy array of shape (n_vertices, 2).

            The {-2,...,2} shell is enough to correctly bound the origin's Voronoi
            cell for the square/rectangular/hexagonal lattices this has been
            checked against. It is not proven to converge for arbitrarily oblique
            (near-degenerate b1/b2 angle) lattices -- if this is ever used for a
            strongly oblique lattice, widen the shell and confirm the extracted
            region stays bounded (no -1 entries in the raw Voronoi region).

        m == 3 raises NotImplementedError: 3D Voronoi-cell construction is
        straightforward with the same call, but rendering a closed, correctly
        wound polyhedron from Qhull's ridge output is real additional work that is
        out of scope for now (pybinding itself does not support 3D BZ plotting
        either).
        """
        b_vectors = self.reciprocal_vectors()
        m = len(b_vectors)

        if m == 1:
            b1 = b_vectors[0]
            return np.array([-b1 / 2, b1 / 2])

        if m == 2:
            b1, b2 = b_vectors
            shell = range(-2, 3)
            points = []
            origin_index = None
            for n1 in shell:
                for n2 in shell:
                    if n1 == 0 and n2 == 0:
                        origin_index = len(points)
                    points.append(n1 * b1 + n2 * b2)
            points = np.array(points)

            vor = Voronoi(points)
            region_index = vor.point_region[origin_index]
            region = vor.regions[region_index]
            if len(region) == 0 or -1 in region:
                raise RuntimeError(
                    "Brillouin zone construction failed: the origin's Voronoi "
                    "region is unbounded or empty within the current candidate "
                    "shell; the reciprocal lattice may be too oblique for the "
                    "{-2,...,2} shell used here -- widen it and retry."
                )

            vertices = vor.vertices[region]
            angles = np.arctan2(vertices[:, 1], vertices[:, 0])
            order = np.argsort(angles)
            return vertices[order]

        raise NotImplementedError("3D Brillouin zone construction is not yet supported")
