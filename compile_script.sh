echo "Compiling KITEx and KITE-tools. Requires a C++17 compiler (GCC 8 or newer)."
WAVEPACKET=1

# create the directory structure
mkdir -p build
mkdir -p tools/build

# Eigen ships with KITE under third_party/eigen3; these are fallbacks in case that's missing
eigen="-I../third_party/eigen3 -I/opt/local/include/eigen3 -I$HOME/include/eigen3 -I/usr/include/eigen3"

# Library and header locations
libraries="-L$HOME/lib -L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5"
headers="-I$HOME/include -I/opt/local/include $eigen" 

# Compile KITEx
include="-I../Src -I../Src/Hamiltonian -I../Src/Tools -I../Src/Vector -I../Src/Simulation -I ../Src/Lattice"
sources="../Src/*.cpp ../Src/*/*.cpp"

echo "Compiling KITEx"
cd build
for i in $sources; do
  echo "Compiling $i"
  g++ -DCOMPILE_WAVEPACKET=$WAVEPACKET -std=c++17 -O2 $i  $include $headers  -fopenmp -c
done 
echo "Linking"
g++ -std=c++17 -O2 *.o $libraries -fopenmp -o ../KITEx
echo "Done."
cd ..

# Compile KITE-tools
include="-I../Src -I../Src/OptCond_1order -I../Src/OptCond_2order -I../Src/Spectral -I../Src/Tools -I ../Src/CondDC"
sources="../Src/*.cpp ../Src/*/*.cpp"

echo "Compiling KITE-tools"
cd tools/build
for i in $sources; do
  echo "Compiling $i"
  g++ -std=c++17 -O2 $i $include $headers -fopenmp -c
done 
echo "Linking"
g++ -std=c++17 -O2 *.o $libraries -fopenmp -o ../../KITE-tools
echo "Done."
cd ../..
