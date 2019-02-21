# Test optimisation

On Mahuika:

```
# load modules
module load CMake
module load GSL/2.3-gimkl-2017a

# build directory
mkdir build
cd build

# compile
cmake ..
make

# run
./testopt
```

## Example output

```
mahuika01:build csco212$ ./testopt
============= Calling find_max_original =============
Have bracket: (4.76837, 5.96046)
Num function calls after bracketing: 20
Optimum found at p = 4.8, f(p) = -1.4499e-11
Total function calls after opt: 51
============= Calling find_max_brent =============
starting points: 0.8, 0.88, 1.00944
Found bracket: 2.99331, 4.42885, 6.75161
             : 3.26412, 0.137749, 3.80878
Num function calls after bracketing: 9
Converged
Found maximum at 4.8, function value -8.5704e-22
Number of iterations: 5
Total GSL Func evals: 16
============= Summary =============
Bisection optimiser used: 51 force evaluations
Brent     optimiser used: 16 force evaluations
```

Over 3x less function evaluations for a simple quadratic function.
