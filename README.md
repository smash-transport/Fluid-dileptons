# Fluid-dileptons

This codebase can be used to calculate the spectra of thermal dielectrons, given some thermodynamic parameters.

It can be used either as a library or as a standalone code.

## Compilation and usage

From the overarching repository run:

```
mkdir build && cd build
cmake ..
cmake --build .
```

This will compile the examplary `main.cpp` and produce the `dilepton_example` executable, which uses the configuration file `dilepton_config.txt`