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

Registering **one operator per sublattice** (rather than a single combined $S_z$) keeps the
checkerboard/edge pattern resolved instead of averaging it away:

``` python
def register_sz_operators(calculation):
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)
    calculation.add_orbital_coupling('Aup', 'Aup', 0.5, 'l0')   # Sz on sublattice A
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bup', 0.5, 'l1')   # Sz on sublattice B
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l1')
    return ['l0', 'l1']

calculation.ldos_map(energy_=energy, sigma_=sigma, vectors_=vectors,
                      operators=register_sz_operators(calculation))
```

This writes `/Calculation/ldos_map/Map_Operators/l0` and `.../l1`, each one value per unit cell
(reshape to `(ly, lx)` to recover the 2D grid — see `afm_zigzag_ribbon_process.py`). Comparing a
clean run (`vacancy_concentration=0.0`) against a ~5% vacancy run (both spin channels removed at
the same site, via two `#!python add_vacancy()` calls on one shared
`#!python kite.StructuralDisorder` instance) shows how a defect locally disturbs the edge spin
pattern:

<figure>
    <img src="../../../assets/images/custom_vertex_operators/afm_zigzag_ribbon.png" style="width: 40em;" />
    <figcaption>Sz_A (circles) and Sz_B (squares) at every site, clean (left) vs. 5% vacancy
    concentration (right). Signed, symmetric-log color scale shared across both panels.</figcaption>
</figure>

### Mapping the global spin balance: `custom_one()` with spin projectors

`#!python ldos_map` is a *local* probe — it doesn't by itself tell you whether the two edges'
opposite polarizations cancel in the total DOS. `afm_zigzag_dos.py` answers that with the
**projector** onto each spin (not $S_z$, which would only give their difference), run through
`#!python calculation.dos()` and `#!python calculation.custom_one()` in the same script:

``` python
def register_spin_projector(calculation, spin):
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)
    if spin == "up":
        calculation.add_orbital_coupling('Aup', 'Aup', 1.0, 'l0')
        calculation.add_orbital_coupling('Bup', 'Bup', 1.0, 'l0')
    else:
        calculation.add_orbital_coupling('Adn', 'Adn', 1.0, 'l0')
        calculation.add_orbital_coupling('Bdn', 'Bdn', 1.0, 'l0')
    return ['l0']
```

Because KITE-tools only reconstructs one `#!python custom_one()` vertex per HDF5 file, spin-up
and spin-down are run as two separate KITEx jobs, each producing its own `.h5`/`.dat` output.

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
