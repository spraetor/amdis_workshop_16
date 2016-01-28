# AMDiS Workshop 2016

AMDiS is a finite element framework for adaptive multidimensional
simulations.

This repository contains all necessary data for the workshop on 
AMDiS programming in August 2016.

## Build the examples

We provide a CMake configuration file CMakeLists.txt and we recommend 
an configuration in a separate directory:

```
cd build
cmake -DAMDIS_DIR=/usr/local/amdis/share/amdis ..
make
cd ..
```

# Run the examples

Run the examples by passing the init-file as first argument.

```
./build/PROG init/INITFILE
```


