## Rashba-Edelstein effect: spin density from a bare-density vertex

The [custom-vertex spin Hall page][custom-vertex-example] builds the **spin current** vertex
$\tfrac12\{\hat v_x,\hat s_z\}$. This page uses the same `#!python custom_two()` Kubo-Bastin
machinery with a **different** vertex — the bare spin density $\hat s_x$, with no velocity
symmetrization — to compute the Rashba-Edelstein effect (REE): the spin-density response
$\chi_{yx}=\partial\langle\hat s_x\rangle/\partial E_y$ to an in-plane electric field, present
only when inversion symmetry is broken.

### The two spin-orbit terms, and why only one gives REE

The lattice is the same spin-doubled honeycomb model as
[the spin Hall page][custom-vertex-example] (`Aup, Bup, Adn, Bdn`), with **two** spin-orbit
terms that can be switched on independently:

- **Kane-Mele intrinsic SOC** ($\lambda_I$, next-nearest-neighbor, spin-diagonal, opposite phase
  per spin) — **inversion-symmetric**. It opens a gap and drives the quantized spin Hall
  response, but by itself must contribute **zero** to any spin-density response to a field: it
  cannot distinguish "spin up here" from "spin up there," so it cannot polarize.
- **Rashba SOC** ($\lambda_R$, nearest-neighbor, **spin-flip**, direction-dependent) —
  **inversion-breaking**. This is the actual source of REE.

$$
\hat H_R = i\lambda_R\sum_{\langle ij\rangle}\hat c_i^\dagger\,(\boldsymbol\sigma\times\hat
d_{ij})_z\,\hat c_j
$$

### The vertex

`#!python examples/rashba_edelstein_graphene.py` defines:

``` python
A = custom.Vertex(moments, [[1.0, "l0"]])   # bare s_x, no velocity symmetrization
B = custom.Vertex(moments, [[1.0, "vy"]])   # charge velocity
```

$\hat s_x=\tfrac12\sigma_x$ is **off-diagonal** in the up/down sublattice basis (unlike $\hat
s_z$, which is diagonal on the spin Hall page) — registered as

``` python
calculation.add_orbital_coupling('Adn', 'Aup', 0.5, 'l0')
calculation.add_orbital_coupling('Aup', 'Adn', 0.5, 'l0')
calculation.add_orbital_coupling('Bdn', 'Bup', 0.5, 'l0')
calculation.add_orbital_coupling('Bup', 'Bdn', 0.5, 'l0')
```

!!! Info "Custom-operator labels are positional, not arbitrary names"

    The trailing digit in a label like `#!python "l0"` is read back by KITEx with `#!cpp
    std::stoi` and used directly as a **0-indexed lookup** into the vector of registered operator
    matrices (`#!cpp act_with_stream`, `#!bash Src/Simulation/Custom/SimulationRankOne.cpp`).
    Since only **one** custom operator is registered in this example, it must be called
    `#!python "l0"` — calling it `#!python "l1"` (as if it were an independently-chosen name)
    segfaults, an out-of-bounds vector access, not merely a cosmetic misnaming.

### Units: chi_yx is in e/h, not e^2/h

KITE's Kubo-Bastin/Gamma2D reconstruction machinery — the numerical prefactors inside
[`#!python conductivity_dc()`][calculation-conductivity_dc] and everything `#!python custom_two()`
shares with it — is calibrated to output $e^2/h$ **only for a genuine two-velocity correlator**
(`#!python NumVelocities=2`, e.g. `#!python conductivity_dc()` itself, or the spin-Hall vertex
pair). Verified directly (not assumed): the Hall conductivity of a $C=1$ Haldane Chern insulator
(`#!python examples/dos_dccond_haldane.py`'s own parameters, disorder included) reconstructs to
$\sigma_{xy}\approx1.01$–$1.02$ at mid-gap — a topological plateau that can only be quantized at
$C\cdot e^2/h$, ruling out $e^2/\hbar$, $e^2/(\pi h)$, and other common conventions.

This page's vertex pair has `#!python NumVelocities=1` (bare $\hat s_x$ contributes zero, $\hat
v_y$ contributes one) — one fewer velocity leg than a genuine conductivity, and that missing leg
is exactly where the charge $e$ entered the Kubo-Bastin current operator. So $\chi_{yx}$ carries
one fewer power of $e$ than $\sigma_{xy}$: **its unit is $e/h$, not $e^2/h$.** Mechanically, this
means `#!python rashba_edelstein_graphene_process.py`'s `#!python edelstein()` needs one additional
factor of `#!python 1/energy_scale` beyond what `#!python kane_mele_spin_hall_process.py`'s
`#!python spin_hall()` applies (that function's `#!python NumVelocities=2` needs no such
correction) — the general rule being an extra factor of `#!bash EnergyScale`$^{N_v-2}$ for a
vertex pair with combined velocity-token count $N_v$.

### The Rashba hopping term: adapted, not re-derived

Rather than re-deriving the honeycomb Rashba bond phases from scratch, this reuses the exact
numerical convention already validated and shipped elsewhere in this codebase —
`#!bash examples/paper/Section_4_B_topology/Input/qahe_disorder.py` and
`#!bash examples/paper/Section_4_E_spintronics/Input/gaussian_wavepacket_only_anderson.py` both
use $\rho=\lambda_R\cdot2i/3$ as the Rashba magnitude, with per-bond coefficients that are a pure
phase times $\rho$, determined by the bond's unit direction $\hat d=(d_x,d_y)$:

$$
t(\hat A_\uparrow\!\to\!\hat B_\downarrow,\,\hat d) = \rho\,(d_y+i d_x), \qquad
t(\hat A_\downarrow\!\to\!\hat B_\uparrow,\,\hat d) = \rho\,(d_y-i d_x)
$$

— cross-checked to reproduce **every one** of both reference files' hard-coded coefficients
exactly, despite the two files using different unit-cell orientations from each other and from
`#!python kane_mele_spin_hall.py`. The formula above (not either file's literal numbers) is what's
reused, applied to `#!python kane_mele_spin_hall.py`'s own, already-verified bond geometry: the
three NN bonds (`relative_index` `[0,0]`, `[1,-1]`, `[0,-1]`) point at $90°,-30°,-150°$.

### Validation

Two structural checks, cheap and independent of any absolute-normalization ambiguity:

1. $\lambda_R\neq0,\ \lambda_I=0$: REE response should be nonzero.
2. $\lambda_R=0,\ \lambda_I\neq0$: REE response should vanish — Kane-Mele alone is
   inversion-symmetric and cannot polarize spin.

<figure>
    <img src="../../../assets/images/custom_vertex_operators/rashba_edelstein_graphene.png" style="width: 34em;" />
    <figcaption>chi_yx(mu) with Rashba SOC on (solid) vs. a Kane-Mele-only control with
    lambda_R=0 (dashed) — the control stays near zero throughout, while turning on Rashba SOC
    produces a clear plateau.</figcaption>
</figure>

Check 2 is a **near-exact cancellation**, not a trivially-zero one: when $\lambda_R=0$ the
up/down spin sectors are exactly decoupled (orthogonal subspaces), so $\Gamma_{mn}$ should vanish
identically — but only in expectation over the stochastic trace, not per random-vector sample. At
`#!python num_random_=1`, the residual sampling noise from this cancellation was comparable in
size to the actual $\lambda_R\neq0$ signal (both were, at that statistics level, artifacts of a
single random vector). Raising `#!python num_random_` to 50 shrinks the $\lambda_R=0$ control
roughly 20-fold while leaving the $\lambda_R\neq0$ plateau essentially unchanged — confirming
genuine averaging convergence, not a wrong vertex or lattice.

!!! example

    Get more familiar with KITE: run [`#!python examples/rashba_edelstein_graphene.py`][ree-example]
    and its post-processing yourself, and try sweeping $\lambda_R$ at fixed $\lambda_I$ to see the
    plateau height scale with the Rashba strength.

[custom-vertex-example]: custom_vertex_operators.md
[ree-example]: https://github.com/quantum-kite/kite-v2/tree/master/examples/rashba_edelstein_graphene.py
[calculation-conductivity_dc]: ../../api/kite.md#calculation-conductivity_dc
