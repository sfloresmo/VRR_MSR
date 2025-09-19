Vector Resonant Relaxation (VRR) and the MSR one-loop closure
-
This project explores statistical closure theory using stellar dynamics as a testbed, focusing on Vector Resonant Relaxation (VRR),
a long-range, non-linear, and correlated relaxation mechanism that drives the reorientation of stellar orbital planes around a supermassive black hole.
We implement the Martin–Siggia–Rose (MSR) formalism at bare and one-loop order via an iterative fixed-point approach. 
The code computes
(i) two-point two-time correlation functions;
(ii) the renormalized three-point interaction vertex;
(iii) three-point three-time correlation functions.

Packages
-
This code is written in Julia. To install Julia: https://julialang.org/downloads/platform/

List of Julia packages required to run the code:
- ArgParse
- BenchmarkTools 
- HDF5 
- OffsetArrays 
- Memoize
- CGcoefficient
- Plots

Parameters
-
The parameters of the run are specified in the file Args.jl.
In particular:
- lmax — cutoff in the harmonic number
- Tc_over_dt — time step of integration
- Tmax_over_Tc — time depth of correlations
- order — order of approximation: order=1 for the bare MSR prediction, order=2 for the one-loop MSR prediction
- init — initial condition for the correlation in the fixed-point search
- iter_init — iteration number used to initialize the fixed-point search
- nb_iter — total number of iterations performed starting from iter_init
- coupling — type of coupling (here, quadrupolar with ℓ = 2)

For the parameters used to compute the numerical predictions in Sec. V of the paper, we used:
lmax=7, Tc_over_dt=80, Tmax_over_Tc=1, order=1 and 2, init=Gaussian, and coupling=Quad.

For the one-loop predictions (order=2), each iteration takes about 62 hours of computing time on 128 cores.
The code is structured such that each iteration outputs a data file.
When iter_init=0, the correlation is initialized as a Gaussian.
After each iteration, iter_init is updated: iter_init = 1, 2, …, nb_iter.

In the following section, we show the command to plot the two-point correlations from the generated data files.
For demonstration purposes, we use downgraded parameters so the computation is faster and can be run locally, in particular we use order=1.

Run
-
I. Open the run folder in the terminal:
```sh
cd run
```
II. Generate the data of the prediction by writing the command (computing time ~ 1 min): 

```sh
julia -t 8 Run.jl --parallel true --lmax 7 --coupling Quad --Tc_over_dt 10 --Tmax_over_Tc 1 --order 1 --init Gaussian --iter_init 0 --nb_iter 10 
```
III. Generate the figure from the data:
```sh
julia Run_Create_Figs.jl --lmax 7 --coupling Quad --Tc_over_dt 10 --Tmax_over_Tc 1 --order 1 --init Gaussian --iter_init 0 --nb_iter 10 
```

Figure
-
