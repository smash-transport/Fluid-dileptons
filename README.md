# Fluid-dileptons

This codebase can be used to calculate spectra of thermal dielectrons in a volume, given the thermodynamic parameters of temperature, baryochemical potential, the 4-volume and fraction of QGP in it, besides the position and velocity of the volume with respect to a calculation frame.

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

To include the library in a codebase, simply add `#include "fluiddileptons.h"` to your main source code.

## Configuration file

Takes in pairs of keys and values, parsed by a `:` character. The available keys and value types are:

| Key           |    Type    |   Default    | Description                          |
| :------------ | :--------: | :----------: | :----------------------------------: |
| `masses`          | 3 doubles  | `0, 2, 0.01` | dilepton invariant mass grid         |
| `abs_momenta`     | 3 doubles  | `0, 3, 0.05` | dilepton absolute momentum grid      |
| `x_range`         | 2 doubles  |    no cut    | x-position acceptance                |
| `y_range`         | 2 doubles  |    no cut    | y-position acceptance                |
| `eta_range`       | 2 doubles  |    no cut    | spatial rapidity acceptance          |
| `pT_range`        | 2 doubles  |    no cut    | transverse momentum acceptance       |
| `yrap_range`      | 2 doubles  |    no cut    | rapidity acceptance                  |
| `oversample`      | `int > 0`  |     `1`      | number of oversampled dileptons      |
| `suppress_output` |   string   |    `none`    | what to NOT calculate: `dilepton`, `spectra`, or `none` |


The grid values (3 doubles) are ordered are (min, max, step), and the range values (2 doubles) are ordered (min, max).

## Output

The code can output either a file with dileptons (mass, q, weight, source), or a file with the spectra histogrammed in mass and momentum, with block for each source.