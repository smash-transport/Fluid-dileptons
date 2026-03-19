# Fluid-dileptons

This codebase can be used to calculate the spectra of thermal dielectrons, given some thermodynamic parameters.

It can be used either as a library or as a standalone code.

## Compilation and usage

To build and compile the code run, from the overarching repository:

```
mkdir build
cd build
cmake ..
cmake --build .
```

This will produce the `libfluid_dileptons.a` library. Alternatively, to link it directly via a CMake project, simply use the same commands as above from the `example` folder. This will compile the examplary source code and compile the `dilepton_example` executable (in the build folder), which uses the configuration file `dilepton_config.txt` and creates the `example_spectra.dat` file.

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

| Key              |    Type   | Default | Description |
| :--------------- | :-------: | :----: | -------: |
| `masses`         | 3 doubles | 0, 2, 0.01 | dilepton invariant mass grid

The grid values (3 doubles) are ordered are (min, max, step), and the range values (2 doubles) are ordered (min, max).