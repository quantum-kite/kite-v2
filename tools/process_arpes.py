"""       
        ##############################################################################      
        #                        KITE | Release  1.1                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2022                 #      
        #                                                                            #      
        ##############################################################################      
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# This file is normally accessed via the examples/process_arpes.py symlink
# (see tools/README.md / examples/arpes_*.py), so cwd is typically examples/
# and `import kite_style` would just work off sys.path[0]. But it can also be
# invoked with its real tools/ path directly (e.g. from a script running with
# cwd == tools/), where examples/ would not otherwise be on sys.path. Resolve
# the shared style module's location relative to *this file's own* location
# (like run_all_examples.py's _this_dir pattern) so both invocation styles work.
_examples_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples")
if _examples_dir not in sys.path:
    sys.path.insert(0, _examples_dir)
import kite_style

kite_style.apply()


def process_arpes(name):
    # open file and read from it
    f = open(name,"r")
    whole = f.readlines()


    # inside that file, find where is each field. There are three fields: 
    # list of energies, list of k-vectors and the arpes matrix
    print("Processing file: finding location of each dataset. ", end="")
    cutarpes, cutenergies, cutvectors = -1, -1, -1
    count = 0
    found = 0
    for i in whole:
      if("ARPES:" in i):
        cutarpes = count + 1
        found += 1
      elif("Energies:" in i):
        cutenergies = count + 1
        found += 1
      elif("k-vectors:" in i):
        cutvectors = count + 1
        found += 1

      count +=1
      if(found == 3):
        break
    print(f"Done. Locations: {cutarpes} {cutenergies} {cutvectors}")

    # process the energies and arpes matrix into numpy arrays
    print("Fetching energies. ", end="")
    energies = np.array([float(i[:-1]) for i in whole[cutenergies:cutarpes - 1]])
    print(f"Done. Number of energies: {len(energies)}")

    # process the arpes matrix
    print("Fetching ARPES matrix. ", end="")
    arpes = np.array([[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]])
    print(f"Done. Arpes matrix dimensions: {len(arpes)} {len(arpes[0])}")

    # find the number of k-vectors
    print("Fetching k vectors. ", end="")
    kvectors = whole[cutvectors:cutenergies - 1]
    num_vectors = len(kvectors)
    nvectors = np.linspace(1, num_vectors, num_vectors)
    print(f"Done. Number of k vectors: {num_vectors}")
    # print(nvectors)

    # plot the matrix
    print("Making meshgrid. ", end="")
    X,Y = np.meshgrid(nvectors, energies)
    print("Done.")

    print("Plotting 2D color map. ", end="")
    fig, axs = plt.subplots(1, 1, figsize=(20, 15))
    # A(k,E) is a dense continuous field, not a small set of discrete series --
    # it needs a perceptually-ordered sequential colormap, not the qualitative
    # brand-color line cycle used elsewhere in these examples (kite_style.KITE_COLORS).
    # kite_spectral_cmap() is KITE-branded in the same spirit (dark -> brand
    # color -> bright highlight, echoing magma/inferno's structure) while
    # staying monotonic in luminance -- unlike the previous 'hot', it won't
    # visually flatten real intensity variation, and unlike 'jet' it can't
    # introduce false banding/edges. See kite_style.py for the luminance check.
    mesh = axs.pcolormesh(X, Y, arpes, cmap=kite_style.kite_spectral_cmap(), rasterized=True, shading="auto")
    cbar = fig.colorbar(mesh, ax=axs)
    cbar.set_label(r"$A(\mathbf{k},E)$ (arb. units)", fontsize=26)
    cbar.ax.tick_params(labelsize=20)
    axs.set_xlabel("k-point index (along path)", fontsize=30)
    axs.set_ylabel("Energy (eV)", fontsize=30)
    print("Done.")

    print("Saving figure.", end="")
    # Figure is saved as a vector PDF (publication-quality), but the mesh
    # itself is rasterized (rasterized=True above) to keep file size sane for
    # a dense k x E grid -- dpi below controls only that rasterized layer's
    # resolution. Raised from the previous 50 dpi (visibly blocky on a
    # 20x15in canvas) to 150: a real, deliberate size/sharpness trade-off,
    # not the 300+ dpi this project otherwise defaults to for raster output --
    # 300 dpi at this figure size (20x15in) produces a much heavier PDF for
    # marginal visual gain at typical (screen/docs) viewing sizes. Bump to
    # 300 in `plt.savefig` below if a print-quality copy is specifically needed.
    plt.savefig("arpes.pdf", format="pdf", dpi=150)
    plt.close(fig)
    print("Done.")

if __name__ == "__main__":
    name = sys.argv[1]
    process_arpes(name)
