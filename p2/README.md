# Sample Code

## Graph class
`graph.hpp` and `graph.cpp` contain a simple class to model unweighed undirected graphs that you may use if you wish. For convenience, the graph already supports output in the DIMACS format.

## Main routine
`main.cpp` contains a toy `main` routine that, for demonstration purposes,
creates a graph that represents a path and outputs it to `stdout`.

## Makefiles
`Makefile` (and `Make.config`) provide a convenient build environment that you may use if you want to.
It automatically compiles all files ending with `.cpp`, `.C` and `.CPP` using the C++ compiler specified in `CXX`, and compiles all files ending with `.c` using the C compiler specified in `CC` and finally links everything together into a binary which is referenced by the symlink `build/main`.

You can create a non-optimized build with debug info (`-O0 -g`) using `make debug`, and create an optimized (`-O3`) build (in order to test runtime) using `make opt`.
You can override all relevant variables (e.g. `CC`, `CXX`, etc.) in `Make.config`.
