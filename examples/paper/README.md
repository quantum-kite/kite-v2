!!! Info

    This file describes the contents of the [`#!bash kite/examples/paper/`][examples-paper-github]-folder.

## Examples from the KITE paper

In the subdirectories, you can find the scripts to generate the figures from the [KITE-paper].

!!! Warning "Not re-run since the paper"

    These scripts build their lattices with `#!python kite.lattice` and
    `#!python kite.repository` rather than pybinding, and plot with matplotlib
    rather than seaborn, so they run in the environment `#!bash environment.yml`
    describes. They have **not** been re-run: the results stored here are the
    ones behind the published figures, and the parameters are the published
    ones, which for some scripts means a very long calculation.

    One exception: `#!python Section_4_A_twisted_bilayer/Input/tblg_DOS.py`
    still needs pybinding, because it loads a pre-built lattice that carries the
    interlayer hoppings.

    `#!python Section_4_E_spintronics` stores its wave-packet results under
    `#!bash Sx`/`#!bash Sy`/`#!bash Sz`, the names KITE wrote when the paper ran.
    A run today registers operators and writes `#!bash l0`/`#!bash l1`/`#!bash l2`
    instead; `#!python plot_spin.py` reads either.

[KITE-paper]: rsos.191809.pdf
[examples-paper-github]: https://github.com/quantum-kite/kite/tree/master/examples/paper