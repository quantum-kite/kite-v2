## Néel-gapped zigzag ribbon: operator-weighted LDOS and spin-resolved DOS

A worked example of `#!python ldos_map(operators=[...])` and per-spin `#!python custom_one()`
on a zigzag-terminated graphene ribbon with a fixed (not self-consistently solved) staggered
on-site mass — opposite sign per sublattice, opposite sign per spin. It exercises three related
API features together: operator-weighted local spectral density, structural vacancy disorder,
and spin-resolved total DOS.

### The model

`#!python examples/afm_zigzag_ribbon.py`'s `#!python afm_zigzag_lattice()` builds a 4-sublattice
honeycomb lattice (`Aup, Bup, Adn, Bdn`) with on-site energies $+\Delta$ on `Aup`/`Bdn` and
$-\Delta$ on `Bup`/`Adn` ($\Delta=0.3\,t$ here), periodic along $\mathbf a_1$ and open along
$\mathbf a_2$ — i.e. a ribbon with two zigzag-terminated edges. This mass pattern gaps the bulk
Dirac cone and pins two flat, chiral-symmetry-protected edge bands at $E=\pm\Delta$, each
localized on one sublattice at one edge and (because the mass flips sign between spins) on one
spin as well — see the script's own docstring for the full derivation and literature references.

### Mapping the local spin density: `ldos_map(operators=[...])`

`#!python ldos_map()` normally computes a single, plain quantity: the local density of states at
every site, at one target energy, via a stochastic (random-vector) KPM estimate. Its
`#!python operators=` parameter generalizes this to $\rho_O(r,E)=-\tfrac1\pi\text{Im}\,
\text{Tr}[O\cdot G(r,r,E)]$ — the same map, but weighted at every site by a Hermitian on-site
operator $O$ instead of the identity, so instead of "how many states" you get "what expectation
value of $O$ do the states at this site/energy carry." Passing `#!python operators=None` (the
default) reproduces the plain LDOS map unchanged; each label in the list adds one *extra* map,
computed at no extra propagation cost as a byproduct of the same stochastic run.

$O$ itself is never passed as a matrix. It's built the same way `#!python custom_one()`/
`#!python custom_two()`/`#!python gaussian_wave_packet()` build *their* operators — two steps:

``` python
def register_sz_operators(calculation):
    # Step 1: give each orbital a name. This just creates a name -> index lookup
    # table on `calculation`; it does not touch the Hamiltonian.
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)

    # Step 2: add_orbital_coupling(row, col, value, label) sets one matrix entry
    # O[row, col] = value in the operator registered under `label`. Call it once
    # per nonzero entry; unset entries default to zero.
    calculation.add_orbital_coupling('Aup', 'Aup', 0.5, 'l0')   # -> l0 = diag(+0.5, 0, -0.5, 0)
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')  #    = Sz restricted to sublattice A
    calculation.add_orbital_coupling('Bup', 'Bup', 0.5, 'l1')   # -> l1 = diag(0, +0.5, 0, -0.5)
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l1')  #    = Sz restricted to sublattice B
    return ['l0', 'l1']  # hand the label list straight to ldos_map's operators=
```

One combined $S_z=\text{diag}(+\tfrac12,+\tfrac12,-\tfrac12,-\tfrac12)$ (uniform across both
sublattices) would still be a valid registration, but it would sum the two sublattices' opposite
edge polarizations together at every site and hide the checkerboard/edge pattern this example is
about — hence two separate labels, `l0` (A-sublattice-only) and `l1` (B-sublattice-only), each
zero on the sublattice it doesn't cover.

With that registration done, the actual call is the same `#!python ldos_map()` you'd write for a
plain map, plus the label list:

``` python
operators = register_sz_operators(calculation)          # -> ['l0', 'l1']
calculation.ldos_map(energy_=energy, sigma_=sigma, vectors_=vectors, operators=operators)
```

`energy_`/`sigma_`/`vectors_` mean exactly what they mean for a plain `#!python ldos_map()` call
(target energy, Gaussian broadening width, number of random vectors) — `operators=` doesn't
change any of that, it only adds extra output. After the run, KITEx writes one array per label:
`/Calculation/ldos_map/Map_Operators/l0` and `.../l1` in the output HDF5 file, each one value per
unit cell (reshape to `(ly, lx)` to recover the 2D grid — see `afm_zigzag_ribbon_process.py`).
Comparing a clean run (`vacancy_concentration=0.0`) against a ~5% vacancy run (both spin channels
removed at the same site, via two `#!python add_vacancy()` calls on one shared
`#!python kite.StructuralDisorder` instance) shows how a defect locally disturbs the edge spin
pattern:

<figure>
    <img src="../../../assets/images/custom_vertex_operators/afm_zigzag_ribbon.png" style="width: 40em;" />
    <figcaption>Sz_A (circles) and Sz_B (squares) at every site, clean (left) vs. 5% vacancy
    concentration (right). Signed, symmetric-log color scale shared across both panels.</figcaption>
</figure>

### Mapping the global spin balance: `custom_one()` with spin projectors

`#!python ldos_map` is a *local* probe (one value per site) — it doesn't by itself tell you
whether the two edges' opposite polarizations cancel when you add up the *whole* ribbon.
`afm_zigzag_dos.py` answers that different question with a different tool:
`#!python calculation.custom_one()` computes $\text{Tr}[A\cdot T_n(\tilde H)]$ — a single number
per Chebyshev moment $n$, summed over the *entire* system, which reconstructs into
$\text{Tr}[A\,\delta(E-H)]$: the **total, energy-resolved density of states weighted by operator
$A$** (no $r$ dependence at all, unlike `ldos_map`).

Here $A$ is deliberately the **projector** onto one spin (weight 1 on that spin's two
sublattices, 0 on the other), not $S_z$: $\text{Tr}[S_z\,\delta(E-H)]$ would only give you
spin-up-minus-spin-down (their *difference*), whereas registering each spin's projector
separately lets you compare spin-up DOS against spin-down DOS *directly*, side by side:

``` python
def register_spin_projector(calculation, spin):
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)
    if spin == "up":
        calculation.add_orbital_coupling('Aup', 'Aup', 1.0, 'l0')   # -> l0 = diag(1,0,0,0)
        calculation.add_orbital_coupling('Bup', 'Bup', 1.0, 'l0')   #    + diag(0,1,0,0)
    else:                                                            #    = P_up
        calculation.add_orbital_coupling('Adn', 'Adn', 1.0, 'l0')   # -> l0 = P_down instead
        calculation.add_orbital_coupling('Bdn', 'Bdn', 1.0, 'l0')
    return ['l0']
```

Calling it looks like this — `#!python custom.Vertex` packages the registered label (with a
coefficient, here `1.0`) into the "vertex" `#!python custom_one()` expects, and `num_moments`
here plays the same role as it does for `#!python ldos()`/`#!python dos()` (Chebyshev expansion
order, i.e. energy resolution):

``` python
from kite import custom

vertex = custom.Vertex(num_moments, [[1.0, "l0"]])
calculation.custom_one(stream_=vertex, num_random_=num_random, num_disorder_=num_disorder)
```

`#!python calculation.dos()` is also called in the same script, unrelated to the operator — it's
just the ordinary, un-weighted total DOS, kept alongside as a reference curve. Because
KITE-tools' `--CustomOne` reconstruction only handles one `#!python custom_one()` vertex per
HDF5 file, spin-up and spin-down are run as two separate KITEx jobs (two calls to
`register_spin_projector`, two output files), each producing its own `.h5`/`.dat` output —
`afm_zigzag_dos_process.py` then loads all four `.dat` files (clean/vacancy × up/down) and plots
them together.

<figure>
    <img src="../../../assets/images/custom_vertex_operators/afm_zigzag_dos.png" style="width: 40em;" />
    <figcaption>Spin-up vs. spin-down total DOS, clean (left) vs. 5% vacancy (right). Clean: the
    two curves coincide at every energy — the opposite edge polarizations cancel globally.
    Vacancy: they visibly split, most strongly at E=&plusmn;&Delta;, the edge-band energy.</figcaption>
</figure>

### Running it

``` bash
# 1. Real-space Sz map: clean vs. 5% vacancy
python afm_zigzag_ribbon.py 0.0   # -> afm_zigzag_ribbon_clean-output.h5
python afm_zigzag_ribbon.py 0.05  # -> afm_zigzag_ribbon_vac5-output.h5
../build/KITEx afm_zigzag_ribbon_clean-output.h5
../build/KITEx afm_zigzag_ribbon_vac5-output.h5
python afm_zigzag_ribbon_process.py afm_zigzag_ribbon_clean-output.h5 afm_zigzag_ribbon_vac5-output.h5

# 2. Spin-resolved total DOS: clean vs. 5% vacancy, up vs. down (4 runs)
python afm_zigzag_dos.py up 0.0
python afm_zigzag_dos.py down 0.0
python afm_zigzag_dos.py up 0.05
python afm_zigzag_dos.py down 0.05
../build/KITEx afm_zigzag_dos_clean_up-output.h5
# ... (repeat KITEx for the other three .h5 files)
../build/KITE-tools afm_zigzag_dos_clean_up-output.h5 --DOS -N dos.dat \
    --CustomOne -E -4 4 1000 -N custom_clean_up.dat
# ... (repeat KITE-tools for the other three, matching custom_<tag>.dat names)
python afm_zigzag_dos_process.py
```

`afm_zigzag_bands.py` is a standalone, plain-numpy exact diagonalization of the same lattice
(bulk and ribbon bands) — no KITEx run needed — useful as a quick independent check of the
gap and edge-band energies before spending time on a full KPM run.

!!! example

    Run [the full script family][afm-ribbon-example] yourself, and try a larger `ly` (more ribbon
    width) or a different `vacancy_concentration` to see how the edge signal and its disorder
    sensitivity scale.

[afm-ribbon-example]: https://github.com/quantum-kite/kite-v2/tree/master/examples/afm_zigzag_ribbon.py
