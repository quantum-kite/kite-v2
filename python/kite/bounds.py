"""Spectral bounds for the Chebyshev rescaling, without pybinding.

`kite.Configuration` needs an interval that contains the spectrum of the
Hamiltonian it is about to export. The estimate here is built additively, term
by term, which is legitimate because of **Weyl's inequality**: for Hermitian A
and B, `lambda_max(A+B) <= lambda_max(A) + lambda_max(B)`. Each contribution
therefore only has to bound the extreme eigenvalue of its own term, and the sum
bounds the total. Because the inequality is attained at a different eigenvector
for each term, the result over-estimates -- which is the safe direction.

The clean-lattice term is exact rather than estimated: `hamiltonian_k` builds
the Bloch Hamiltonian from the lattice's own hoppings, so diagonalising it on a
dense grid gives the true band extrema. That is better than the bound this
function replaced, which read pybinding's KPM scaling factors.

Contributions, following the project's own rules:

    clean bands          BZ scan of H(k), exact
    with magnetic field  Gershgorin instead -- Hofstadter broadening exceeds
                         the zero-field band edges, so the scan is invalid
    uniform disorder     mu +- sqrt(3)*sigma   (sigma is a standard deviation;
                         the box half-width is sqrt(3) times larger)
    gaussian disorder    mu +- lambda*sigma, lambda = 8 by default
    deterministic        mu
    vacancies            zero -- a principal submatrix, which by Cauchy
                         interlacing can only narrow the spectrum
    structural defects   max over nodes of |dU| + sum|dt|, since the supplied
                         values are increments on top of the clean hopping

The returned interval carries a margin so the spectrum does not sit on +-1,
where the Chebyshev reconstruction weight 1/sqrt(1-x^2) diverges.
"""

import numpy as np

__all__ = ["estimate_bounds", "clean_band_extrema", "gershgorin_bound"]

#: Gaussian tails are unbounded, so a probabilistic cut is needed. lambda = 8
#: covers any run KITE can execute (it corresponds to N_eff beyond 1e13 at
#: p = 1e-3) and costs one extra sigma of energy scale.
DEFAULT_GAUSSIAN_LAMBDA = 8.0

#: Fractional headroom left between the spectrum and +-1.
DEFAULT_MARGIN = 0.01


def clean_band_extrema(lattice, grid=24):
    """Exact band extrema of the clean lattice, by scanning the Brillouin zone.

    `grid` points per reciprocal direction. The extrema of a tight-binding band
    structure are attained at high-symmetry points, which a uniform grid over
    the reciprocal cell hits or comes close to; the margin absorbs the rest.
    """
    from .visualize import hamiltonian_k

    vectors = np.asarray(lattice.vectors, dtype=float)
    dim = vectors.shape[0]

    # reciprocal vectors: B = 2*pi*(A^-1)^T for the sub-block actually spanned
    square = vectors[:, :dim]
    reciprocal = 2 * np.pi * np.linalg.inv(square).T

    axes = [np.linspace(0, 1, grid, endpoint=False)] * dim
    lo, hi = np.inf, -np.inf
    for fractional in np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, dim):
        k = fractional @ reciprocal
        eigenvalues = np.linalg.eigvalsh(hamiltonian_k(lattice, k))
        lo = min(lo, eigenvalues[0])
        hi = max(hi, eigenvalues[-1])
    return float(lo), float(hi)


def gershgorin_bound(lattice):
    """Row-sum bound: max over rows of |onsite| + sum of |hopping|.

    Valid whatever the Hamiltonian does to the bands -- used when a magnetic
    field makes the zero-field band edges meaningless.
    """
    nsub = lattice.nsub
    onsite = np.zeros(nsub)
    for sub in lattice.sublattices.values():
        onsite[sub.alias_id] = np.max(np.abs(np.asarray(sub.energy)))

    row_sum = np.zeros(nsub)
    for family in lattice.hoppings.values():
        magnitude = np.max(np.abs(np.asarray(family.energy)))
        for term in family.terms:
            # each stored term also implies its conjugate
            row_sum[term.from_id] += magnitude
            row_sum[term.to_id] += magnitude
    return float(np.max(onsite + row_sum))


def _disorder_spread(disorder, gaussian_lambda):
    """Half-width contributed by on-site disorder, and its mean shift."""
    if disorder is None:
        return 0.0, 0.0
    types = np.asarray(getattr(disorder, "_type_id", []), dtype=object).ravel()
    means = np.asarray(getattr(disorder, "_mean", []), dtype=float).ravel()
    stdvs = np.asarray(getattr(disorder, "_stdv", []), dtype=float).ravel()
    if means.size == 0:
        return 0.0, 0.0

    # type ids are set in Disorder.__init__: Gaussian 1, Uniform 2, Deterministic 3
    GAUSSIAN, UNIFORM, DETERMINISTIC = 1, 2, 3

    spread = 0.0
    for i in range(means.size):
        sigma = float(stdvs[i]) if i < stdvs.size else 0.0
        try:
            kind = int(types[i]) if i < types.size else UNIFORM
        except (TypeError, ValueError):
            kind = UNIFORM
        if kind == GAUSSIAN:
            spread = max(spread, gaussian_lambda * sigma)
        elif kind == DETERMINISTIC:
            pass                                # a shift only, carried by the mean
        else:                                   # uniform, and the safe default
            spread = max(spread, np.sqrt(3.0) * sigma)
    return spread, float(np.max(np.abs(means)))


def _structural_spread(disorder_structural):
    """Largest node increment: |dU| + sum|dt|, per RULES 8.1.

    The values supplied to `add_structural_disorder` are increments added on
    top of the clean hopping, so the bound is a max over nodes, not a sum --
    concentration does not enter it.
    """
    if disorder_structural is None:
        return 0.0
    items = disorder_structural if isinstance(disorder_structural, (list, tuple)) \
        else [disorder_structural]
    worst = 0.0
    for item in items:
        onsite = np.abs(np.asarray(getattr(item, "_disorder_onsite", []), dtype=float)).ravel()
        hopping = np.abs(np.asarray(getattr(item, "_disorder_hopping", []),
                                    dtype=complex)).ravel()
        worst = max(worst, float(onsite.sum() + hopping.sum()) if onsite.size or hopping.size else 0.0)
    return worst


def estimate_bounds(lattice, disorder=None, disorder_structural=None,
                    magnetic_field=False, gaussian_lambda=DEFAULT_GAUSSIAN_LAMBDA,
                    margin=DEFAULT_MARGIN, grid=24):
    """Return `(E_min, E_max)` containing the spectrum, with margin.

    Parameters
    ----------
    lattice : kite.lattice.Lattice
    disorder : kite.Disorder, optional
    disorder_structural : kite.StructuralDisorder or list, optional
    magnetic_field : bool
        Set when a field will be applied. The clean band edges are then not a
        bound -- Hofstadter broadening exceeds them -- so a Gershgorin row sum
        is used instead.
    gaussian_lambda : float
        Tail cut for Gaussian disorder. See the module docstring.
    margin : float
        Fractional headroom, so the spectrum does not land on +-1.
    """
    if magnetic_field:
        radius = gershgorin_bound(lattice)
        lo, hi = -radius, radius
    else:
        lo, hi = clean_band_extrema(lattice, grid=grid)

    spread, shift = _disorder_spread(disorder, gaussian_lambda)
    spread += _structural_spread(disorder_structural)

    lo -= spread + shift
    hi += spread + shift

    centre, half = 0.5 * (hi + lo), 0.5 * (hi - lo)
    half *= 1.0 + margin
    return centre - half, centre + half
