mkdir build_gcc_libint_openmp
cp -r Test build_gcc_libint_openmp
cd build_gcc_libint_openmp
cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_FLAGS='-w -O3 -std=c++11 -I $(PWD)/../deps/eigen -I $(PWD)/../deps/BTAS' -DUSE_LIBINT=ON -DUSE_OMP=ON ..
make VERBOSE=1 -j4
cp -r ../Test/ .
cd ..
