# FSAE MPC

FYP project to design an MPC controller for a student driverless FSAE racecar.

## Installation

As there are submodules in this repository, make sure to clone it using the `--recursive` tag:
```bash
git clone https://github.com/kerry-he/crane-vision.git --recursive
```
If the repository has already been cloned, then the submodules can be fetched using:
```bash
git submodule update --init
```

Step 1: Follow instructions to install [`qpOASES`](https://github.com/coin-or/qpOASES). The manual for the package can be found [here](https://www.coin-or.org/qpOASES/doc/3.0/manual.pdf).

**MATLAB:**
In the MATLAB terminal, run the following commands from within the repository directory.

```matlab
mex -setup C++ % Chooses a C++ compiler. Note you may need to install a compiler if one isn't installed already.
cd packages/qpOASES/interfaces/matlab
make % Runs the compilation script
cd ../../../../
```