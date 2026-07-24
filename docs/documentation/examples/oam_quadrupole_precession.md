## Orbital angular momentum and quadrupole precession: coupled real-time dynamics

`#!python examples/paper/Section_4_E_spintronics`'s wave-packet examples track real-time spin
precession $\langle S(t)\rangle$ under `#!python calculation.gaussian_wave_packet()`. This page
uses the same real-time Chebyshev propagation of $e^{-iHt}$, but tracks orbital angular momentum
$L$ and orbital quadrupoles $Q_{ab}=\tfrac12\{L_a,L_b\}$ instead — a different physics regime from
the disorder-averaged Kubo-Bastin transport examples elsewhere on this site: this is coherent,
single-particle real-time dynamics, not a DC transport coefficient.

### Tracking an arbitrary operator, not just spin

`#!python gaussian_wave_packet()` has no built-in notion of spin: any operator is tracked through
[`#!python add_orbital_coupling()`][calculation-add_orbital_coupling]'s registration mechanism,
the same one `#!python custom_one()`/`#!python custom_two()` use, passed in as
`#!python operators=[...]`. This has one hard restriction: the operator must be **on-site** — a
single matrix acting on one site's orbitals, applied identically at every site (see the warning
under [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling]). Spin, orbital angular
momentum, and orbital quadrupoles all satisfy this; an operator that couples different lattice
sites cannot be tracked this way.

### The model: an sp³ manifold

$p_x,p_y,p_z$ orbitals at a single site are exactly the $l=1$ angular-momentum representation, so
an sp³ Slater-Koster square lattice (`#!python s, px, py, pz` on-site orbitals, spin-orbit coupling
off) is the natural minimal model — the same lattice used elsewhere for an orbital Hall effect
calculation. $L_x, L_y, L_z$ are built as on-site $4\times4$ matrices (zero on/with the $s$
orbital) via

$$
(L_k)_{jl} = -i\,\epsilon_{kjl}
$$

on the $(p_x,p_y,p_z)$ block — verified Hermitian, satisfying $[L_i,L_j]=i\epsilon_{ijk}L_k$ and
$L^2=2\mathbb 1$ (the standard $l=1$ algebra). The six $Q_{ab}=\tfrac12\{L_a,L_b\}$ follow directly
by matrix multiplication. (This differs by an overall sign from the $L_z$-like operator used in the
orbital Hall example — as with spin, the overall sign of an angular momentum operator is a
convention, not a physical ambiguity.)

### Why L and Q are coupled, in general

$Q_{ab}$'s Heisenberg torque, $\mathcal T_{Q_{ab}}=\tfrac{i}{\hbar}[H,Q_{ab}]$, is generally nonzero
whenever the Hamiltonian's hopping is anisotropic between orbital characters. Here
$V_{pp\sigma}\neq V_{pp\pi}$ breaks the continuous rotational symmetry that would otherwise decouple
$L$ and $Q$. Whether that shows up as fast driven precession or something much subtler depends
entirely on which part of the band structure the wave packet actually samples — there is no Fermi
energy in this calculation (`#!python gaussian_wave_packet()` propagates one coherent single-particle
state, not a thermal/Fermi-averaged quantity); the analogous choice is the seed **k-vector**.

### Validation: a clean k=0 seed, and why almost nothing moves

The wave packet here is seeded at $k=0$ with a **pure** $L_z=+1$ state, $(p_x+ip_y)/\sqrt2$ — no
$s$ or $p_z$ admixture. This is deliberately the least ambiguous starting point, but at exactly
$k=0$ this model's Hamiltonian is diagonal in the $(s,p_x,p_y,p_z)$ basis (the odd-parity $s$-$p$
hopping cancels there, and $p$-$p$ hopping doesn't mix orbitals), with $p_x,p_y$ **exactly
degenerate** — each orbital sees one $\sigma$-type and one $\pi$-type bond direction, just swapped.
A combination of degenerate states is itself stationary, so a true $k=0$ plane wave in this state
would show *no* dynamics at all.

What's actually observed: $\langle L_z(t)\rangle$ decays slowly (~10% over the propagation window)
while $\langle L_x\rangle,\langle L_y\rangle$ and every off-diagonal $\langle Q_{ab}\rangle$
($Q_{xy}, Q_{xz}, Q_{yz}$) stay pinned at **exactly zero** throughout — both facts expected, not
bugs (checked by `cmt-physicist`):

- The off-diagonal operators all connect the populated $m=+1$ sector to the $m=0$ ($p_z$) or $m=-1$
  sectors — both completely unpopulated here — so none of them can pick up any expectation value,
  at any time.
- The slow $L_z$ decay is a finite-size-in-$k$-space effect: a spatially localized wave packet is
  necessarily a superposition of many $k$-vectors around $k=0$, not a single Bloch state, and away
  from $k=0$ the $p_x/p_y$ degeneracy lifts. So the packet slowly dephases out of the pure
  $L_z=+1$ eigenspace even though a true $k=0$ plane wave would conserve it exactly.

Seeding the wave packet at a $k$-vector away from $\Gamma$ (where $p_x,p_y$ are genuinely split by
the hopping anisotropy) would show faster, driven precession instead of this slow dephasing — a
natural next experiment, not done on this page since the point here is the cleanest,
minimal-assumption case.

<figure>
    <img src="../../../assets/images/custom_vertex_operators/oam_quadrupole_precession.png" style="width: 34em;" />
    <figcaption>A pure L_z=+1 seed at k=0: L_z (and correspondingly Q_zz, Q_xx=Q_yy) decay slowly
    from k-space dephasing, while L_x, L_y and every off-diagonal Q_ab stay pinned at exactly
    zero -- selection rules, not bugs (see Validation above).</figcaption>
</figure>

!!! example

    Get more familiar with KITE: run [`#!python examples/oam_quadrupole_precession.py`][oam-example]
    and its post-processing yourself, and try seeding the wave packet at a k-vector away from
    $\Gamma$ to see faster, driven L-Q precession instead of slow dephasing.

[calculation-add_orbital_coupling]: ../../api/kite.md#calculation-add_orbital_coupling
[oam-example]: https://github.com/quantum-kite/kite-v2/tree/master/examples/oam_quadrupole_precession.py
