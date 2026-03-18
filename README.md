# Fluid-dileptons

This codebase can be used to calculate the spectra of thermal dielectrons, given some thermodynamic parameters.

It can be used either as a library or as a standalone code.

## Compilation and usage

From the overarching repository run:

```
mkdir build
cd build
cmake ..
cmake --build .
```

This will compile the examplary `tools/example.cpp` and produce the `dilepton_example` executable (in the build folder), which uses the configuration file `tools/dilepton_config.txt`

## Configuration file

Takes in pairs of keys and values, parsed by a `:` character. The available keys and value types are:

- `masses` (3 doubles): dilepton invariant mass grid.
- `abs_momenta` (3 doubles): dilepton absolute momentum grid.
- `x_range` (2 doubles): x-position acceptance.
- `y_range` (2 doubles): y-position acceptance.
- `eta_range` (2 doubles): spatial rapidity acceptance.
- `pT_range` (2 doubles): transverse momentum acceptance.
- `yrap_range` (2 doubles): rapidity acceptance.
- `oversample` (int > 0): number of oversampled dileptons.

The grid values (3 doubles) are ordered are (min, max, step), and the range values (2 doubles) are ordered (min, max).