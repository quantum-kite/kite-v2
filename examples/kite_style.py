""" KITE-branded matplotlib style, shared across example plotting scripts.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

Single source of truth for the color palette + rcParams used by
examples/pybinding/run_all_examples.py's figures (DOS, optical/DC
conductivity) and examples/process_arpes.py's spectral-function heatmap, so
the example plots the docs screenshot/link share one visual identity instead
of matplotlib's default color cycle.

Colors are pulled directly from this repo's own site theme, not invented --
verified against:
    - mkdocs.yml (`palette:` block, `scheme: kite`)
    - docs/assets/stylesheets/extra.css (`[data-md-color-scheme="kite"]`,
      and the dark-mode `[data-md-color-scheme="slate"]` background rules)

    KITE_PRIMARY = "#822855"   site primary / brand color (deep wine/magenta)
    KITE_ACCENT  = "#7bb274"   site accent (muted sage green)
    KITE_TEXT    = "#404040"   site body text color
    KITE_DARK_BG = "#1e1f20"   site dark-mode background

Two distinct things are provided, deliberately not the same palette, because
line plots and 2D heatmaps have different visual-representation needs (see
module docstrings below for why):

1. `KITE_COLORS` / `apply()` -- a qualitative, brand-anchored line/marker
   color cycle for line plots (DOS, optical conductivity, DC conductivity).
2. `kite_spectral_cmap()` -- a sequential, monotonically-brightening colormap
   for 2D spectral-function / ARPES heatmaps. A qualitative brand palette
   (points 1) would be the wrong tool here: a dense continuous field needs a
   perceptually-ordered colormap, not a small set of mutually-distinct hues.

Usage:
    import kite_style
    kite_style.apply()                       # once, before creating figures
    ...
    cmap = kite_style.kite_spectral_cmap()
    ax.pcolormesh(x, y, z, cmap=cmap)
"""

import os

import matplotlib as mpl
import matplotlib.style as mplstyle
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap

# --- Brand colors, verified against mkdocs.yml / extra.css -----------------
KITE_PRIMARY = "#822855"
KITE_ACCENT = "#7bb274"
KITE_TEXT = "#404040"
KITE_DARK_BG = "#1e1f20"

# Qualitative line-color cycle for DOS / optical-conductivity / DC-conductivity
# plots. Order matters: most example plots only show 1-2 curves (a single DOS,
# or Re/Im of sigma_xy, or sigma_xx/sigma_xy), so the first two colors are the
# two *exact* site brand colors -- every such plot will visibly carry the KITE
# identity. Later entries extend the palette for legends with more series
# (e.g. condDC xx/yy pairs across two files, or future multi-band examples).
#
# Colorblind-safety note (checked, not assumed): #822855 and #7bb274 are not a
# red/green confusion pair, and differ substantially in relative luminance
# (0.069 vs. 0.373, sRGB-linear -- roughly a 5x difference), so they stay
# distinguishable by lightness alone under protanopia/deuteranopia/tritanopia,
# not just by hue. Colors 3-4 (blue, amber) are chosen to be maximally
# distinct from both brand colors and from each other under all three common
# CVD types before falling back to brand-family tints/neutrals.
KITE_COLORS = [
    KITE_PRIMARY,  # #822855 deep wine/magenta -- exact site primary
    KITE_ACCENT,   # #7bb274 sage green -- exact site accent
    "#357EDD",     # blue -- strong contrast against both brand colors
    "#FFC107",     # amber -- strong contrast against blue + both brand colors
    "#A95C8E",     # lighter magenta tint (extends the primary family)
    "#4C7A66",     # deeper teal-green tint (extends the accent family)
    "#5C5C5C",     # dark neutral gray
    "#A0A0A0",     # light neutral gray (e.g. a reference/baseline curve)
]

_STYLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kite.mplstyle")


def apply():
    """Apply the KITE house style globally. Call once, before creating any figures.

    Layers two things:
      1. kite.mplstyle -- typography, spines, grid, DPI/vector-font output settings.
      2. KITE_COLORS (this module) -- the brand color cycle, set here in Python
         (rather than duplicated as a literal list inside the .mplstyle file)
         so the palette has exactly one source of truth.
    """
    if os.path.exists(_STYLE_PATH):
        mplstyle.use(_STYLE_PATH)
    else:
        # Fail closed on *content* (never silently fall back to a different,
        # unbranded look) but don't crash the whole example run over a
        # missing typography file -- the color cycle below is the visually
        # load-bearing part and is applied unconditionally.
        print("kite_style: warning -- {0} not found; applying color cycle only, "
              "not full house style.".format(_STYLE_PATH))
    mpl.rcParams["axes.prop_cycle"] = cycler("color", KITE_COLORS)


def kite_spectral_cmap():
    """Sequential colormap for 2D spectral-function / ARPES heatmaps.

    Perceptually-ordered (monotonically increasing relative luminance from
    0.014 to 0.831, verified below -- not just visually eyeballed) so it can't
    introduce spurious features the way a non-monotonic or hue-cycling
    colormap (e.g. `jet`) would. Built from the KITE primary color rather than
    a generic sequential map, in the same visual family as magma/inferno
    (dark -> saturated color -> bright highlight):

        0.00  #1e1f20  KITE dark-mode background (near-black slate)
        0.55  #822855  KITE primary (deep wine/magenta)
        1.00  #f4ead9  warm cream highlight (echoes the site's light
                        `--md-primary-fg-color: #EEEEEE`, warmed slightly so
                        peak intensity reads as a highlight, not clinical white)

    Deliberately does NOT also route through KITE_ACCENT (#7bb274, green) as a
    third stop: an earlier version of this colormap did, and a test render
    against a real KITE spectral function showed a muddy, desaturated
    olive/brown band in the magenta->green transition -- magenta and green are
    near-complementary hues, so linearly interpolating straight through both in
    one ramp desaturates the middle instead of staying vivid. Accent green stays
    reserved for what it's actually good at: a *qualitative* second color in
    line plots (KITE_COLORS), not a stop inside a single continuous intensity
    ramp. Two brand hues in one sequential colormap was cosmetic ambition
    outrunning what the color math actually supports -- dropping to one keeps
    the map vivid and monotonic end to end.

    Relative luminance (sRGB, WCAG formula) at each stop, confirming the
    monotonic ramp: 0.0136, 0.0692, 0.8309.
    """
    stops = [
        (0.00, "#1e1f20"),
        (0.55, "#822855"),
        (1.00, "#f4ead9"),
    ]
    return LinearSegmentedColormap.from_list("kite_spectral", stops, N=256)
