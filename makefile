PREFIX  := build
CXX     := g++
OBJDIR  := $(PREFIX)/.obj

# Eigen-Dir: uses the bundled copy in third_party/eigen3 by default (see
# third_party/eigen3/VENDORED.md). To use your own system Eigen3 instead,
# override on the command line, e.g.:
#   make EIGEN_INC=-I/usr/local/include/eigen3
EIGEN_INC := -Ithird_party/eigen3

# HDF5-Flags
HDF5_LIBS   := -lhdf5_cpp -lhdf5

# FFTW3 -- only float/double are linked; long-double FFT/spectral support
# throws a clear runtime error instead (see Src/FFT/TraitsFFTW.hpp), since
# the long-double variant isn't built by default by common package managers.
FFT_LIBS   := -lfftw3f -lfftw3

CPPFLAGS := $(EIGEN_INC) \
	-Itools/Src \
	-ISrc \
	-ISrc/Tools \
	-ISrc/Lattice \
	-ISrc/Hamiltonian \
	-ISrc/Vector \
	-ISrc/Simulation \
	-DCOMPILE_WAVEPACKET=0 \
	$(HDF5_FLAGS)

CXXFLAGS := -std=c++17 -O3 -fopenmp
LDFLAGS  := -fopenmp $(HDF5_LIBS) $(FFT_LIBS)

SRC_KITEX := $(shell find Src -name '*.cpp')
SRC_TOOLS := $(shell find tools/Src -name '*.cpp')

OBJ_KITEX := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_KITEX))
OBJ_TOOLS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_TOOLS))

all: KITEx KITE-tools

KITEx: $(OBJ_KITEX)
	@mkdir -p $(PREFIX)
	$(CXX) $(LDFLAGS) $^ -o $(PREFIX)/$@

KITE-tools: $(OBJ_TOOLS)
	@mkdir -p $(PREFIX)
	$(CXX) $(LDFLAGS) $^ -o $(PREFIX)/$@

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ_KITEX:.o=.d) $(OBJ_TOOLS:.o=.d)

clean:
	rm -rf $(PREFIX)

.PHONY: all clean
