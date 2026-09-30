# OCP Hessian benchmark report

Results from a Release build on an Intel Core i7-10850H, 30 September 2026.

## Build and method

- 30 problems, one class and header per problem; exact, block BFGS, and L-BFGS.
- GCC 11.4, `-O3 -DNDEBUG -std=c++17`, Release Fatrop and BLASFEO; Intel Core i7-10850H.
- The build-time generator differentiates each distinct stage expression and emits ordinary C++ kernels. The default benchmark uses no runtime AD graph or derivative sweeps.
- Jacobians use forward differentiation; cost gradients use adjoints; Hessians use forward differentiation of the adjoint expressions. Controls and states are differentiated; multipliers and objective scale remain parameters.
- Model configuration is fixed at generation time. Rebuild after changing model equations or constants. Every stage is traced; identical stage kernels share compiled code.
- One warmup and seven measured solves per successful row; the sample with median total time is reported. Failed rows report the first failed solve.
- Eval includes all OCP value and derivative callbacks. Fatrop is total solve time minus Eval. Both include timer overhead. A clock pair costs about 30 ns on this machine, which is material for the smallest generated kernels.
- Iteration cap: 1,000, except 3,000 for minimum-time bicycle obstacle.

```sh
cmake -S . -B build/benchmarks -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_BENCHMARKS=ON -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF
cmake --build build/benchmarks --target \
  fatrop_ocp_bfgs_benchmark fatrop_ocp_derivative_check -j
build/benchmarks/benchmarks/fatrop_ocp_derivative_check
build/benchmarks/benchmarks/fatrop_ocp_bfgs_benchmark 7
```

The benchmark accepts an optional problem name after the repeat count, such as
`pendulum`. `BUILD_BENCHMARKS` is off by default. The separate Ipopt comparison
requires Ipopt's pkg-config package and is enabled with
`-DBENCHMARKS_WITH_IPOPT=ON`; build and run
`fatrop_ocp_vs_ipopt_benchmark` from the same build directory. It compares the
solvers on a configurable large-state OCP, separate from the 30 scenarios below.

## Validation

- Suite: 148,085 finite-difference checks, zero failures; largest scaled error 1.355e-09.
- Expanded test: 444,255 finite-difference checks at three perturbed points per stage, zero failures. Generated Hessians agree with runtime AD within 1.90e-16 scaled error with zero, positive, and negative objective scales and varying multipliers.

## Comparison with runtime AD

These are medians of per-problem speedup ratios, restricted to cases successful with both backends. Callback comparisons use time per call to account for iteration changes. Whole-solve comparisons use total time. Ratios above one mean the compiled backend is faster.

| Method | Hessian callback | Dynamics Jacobian | Cost gradient | All evaluation time | Total solve time |
|---|---:|---:|---:|---:|---:|
| Exact | 34.64× | 2.53× | 6.12× | 6.06× | 1.47× |
| Block BFGS | — | 2.54× | 5.76× | 1.95× | 1.07× |
| L-BFGS | — | 2.60× | 6.08× | 1.79× | 1.05× |

Exact and L-BFGS converged on all 30 problems. Block BFGS reached the iteration cap on minimum-time unicycle and minimum-time quadrotor obstacle, and failed its line search on cart-pole swing-up. The cart-pole failure is new relative to the runtime baseline. Symbolic simplification changes floating-point evaluation order; several sensitive trajectories and iteration counts changed despite the derivative checks passing.

Callback speedups do not imply a proportional end-to-end speedup; solver work, numerical trajectories, and run-to-run timing variation matter.

## Optional Ipopt comparison

`fatrop_ocp_vs_ipopt_benchmark` uses a separate 60-stage large-state OCP with
the same model and starting point in both solvers. Ipopt and Fatrop use L-BFGS
for the comparison; Fatrop exact is shown as a reference. Times are median
milliseconds from five measured runs after two warmups.

| States | Fatrop exact | Fatrop L-BFGS | Ipopt L-BFGS | Iterations: exact / L-BFGS / Ipopt |
|---:|---:|---:|---:|---:|
| 16 | 2.237 | 23.991 | 420.849 | 7 / 24 / 24 |
| 64 | 32.633 | 175.525 | 1048.432 | 11 / 24 / 24 |

## Results

Times are milliseconds. Objective and violation for failed solves describe the final candidate.

| Problem | Method | Status | Iterations | Total | Eval | Fatrop | Eval calls | Objective | Max violation |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| pendulum | Exact | OK | 11 | 0.908 | 0.191 | 0.717 | 6,575 | 14.067379884 | 1.176e-08 |
| pendulum | Block BFGS | OK | 16 | 1.523 | 0.197 | 1.327 | 8,305 | 14.067379884 | 1.184e-08 |
| pendulum | L-BFGS | OK | 32 | 6.005 | 0.394 | 5.611 | 16,192 | 14.067379884 | 1.184e-08 |
| coupled-chain | Exact | OK | 8 | 2.016 | 0.160 | 1.856 | 4,131 | 193.235108790 | 5.274e-16 |
| coupled-chain | Block BFGS | OK | 11 | 2.512 | 0.167 | 2.345 | 4,925 | 193.235108790 | 2.498e-15 |
| coupled-chain | L-BFGS | OK | 27 | 12.999 | 0.397 | 12.602 | 11,691 | 193.235108790 | 1.110e-16 |
| saturating-actuator | Exact | OK | 6 | 0.647 | 0.158 | 0.489 | 3,235 | 4.491738098 | 1.217e-11 |
| saturating-actuator | Block BFGS | OK | 13 | 1.645 | 0.267 | 1.378 | 5,721 | 4.491738098 | 9.659e-15 |
| saturating-actuator | L-BFGS | OK | 27 | 7.659 | 0.552 | 7.106 | 11,691 | 4.491738098 | 1.221e-15 |
| minimum-time-transfer | Exact | OK | 57 | 4.706 | 0.316 | 4.389 | 16,988 | 1.937388422 | 3.543e-13 |
| minimum-time-transfer | Block BFGS | OK | 288 | 18.003 | 1.312 | 16.691 | 71,555 | 1.944071301 | 4.924e-09 |
| minimum-time-transfer | L-BFGS | OK | 98 | 11.072 | 0.478 | 10.595 | 26,537 | 1.937939646 | 4.693e-09 |
| minimum-time-unicycle | Exact | OK | 49 | 5.431 | 0.373 | 5.058 | 13,911 | 2.067734414 | 3.593e-09 |
| minimum-time-unicycle | Block BFGS | Iteration cap | 1001 | 83.456 | 5.874 | 77.582 | 246,953 | 2.084089877 | 2.220e-16 |
| minimum-time-unicycle | L-BFGS | OK | 84 | 13.014 | 0.528 | 12.486 | 21,863 | 2.081759722 | 4.009e-09 |
| van-der-pol | Exact | OK | 13 | 1.889 | 0.141 | 1.747 | 7,651 | 20.527667945 | 2.220e-16 |
| van-der-pol | Block BFGS | OK | 13 | 1.293 | 0.124 | 1.170 | 6,871 | 20.527667943 | 2.220e-16 |
| van-der-pol | L-BFGS | OK | 32 | 6.246 | 0.315 | 5.931 | 17,148 | 20.527667945 | 2.220e-16 |
| nonlinear-path-limit | Exact | OK | 7 | 0.497 | 0.066 | 0.431 | 3,683 | 6.871959878 | 7.572e-09 |
| nonlinear-path-limit | Block BFGS | OK | 12 | 1.080 | 0.092 | 0.988 | 5,323 | 6.871959876 | 7.548e-09 |
| nonlinear-path-limit | L-BFGS | OK | 11 | 1.620 | 0.086 | 1.534 | 4,925 | 6.871959876 | 7.548e-09 |
| predator-prey-harvest | Exact | OK | 11 | 1.643 | 0.141 | 1.502 | 7,675 | 0.778658206 | 1.110e-16 |
| predator-prey-harvest | Block BFGS | OK | 16 | 1.710 | 0.171 | 1.540 | 9,695 | 0.778658206 | 8.771e-15 |
| predator-prey-harvest | L-BFGS | OK | 27 | 5.944 | 0.308 | 5.637 | 16,112 | 0.778658206 | 1.887e-15 |
| minimum-time-powered-ascent | Exact | OK | 89 | 11.437 | 0.765 | 10.672 | 33,437 | 3.457894721 | 7.696e-09 |
| minimum-time-powered-ascent | Block BFGS | OK | 116 | 8.539 | 0.816 | 7.723 | 38,264 | 3.457841551 | 1.554e-08 |
| minimum-time-powered-ascent | L-BFGS | OK | 57 | 7.844 | 0.414 | 7.430 | 19,682 | 3.457881411 | 7.696e-09 |
| duffing-well-transfer | Exact | OK | 38 | 4.394 | 0.558 | 3.835 | 30,775 | 8.603645950 | 8.601e-09 |
| duffing-well-transfer | Block BFGS | OK | 160 | 17.498 | 1.925 | 15.573 | 105,863 | 8.603645985 | 6.144e-09 |
| duffing-well-transfer | L-BFGS | OK | 63 | 9.983 | 0.775 | 9.207 | 43,524 | 8.603645950 | 8.601e-09 |
| lorenz-stabilization | Exact | OK | 8 | 0.960 | 0.096 | 0.864 | 5,044 | 16.710067192 | 2.776e-17 |
| lorenz-stabilization | Block BFGS | OK | 8 | 0.955 | 0.082 | 0.873 | 4,556 | 16.710067192 | 1.749e-15 |
| lorenz-stabilization | L-BFGS | OK | 13 | 2.612 | 0.129 | 2.483 | 6,986 | 16.710067192 | 2.776e-17 |
| exothermic-reactor | Exact | OK | 6 | 1.054 | 0.118 | 0.936 | 5,250 | 9.086529004 | 1.854e-14 |
| exothermic-reactor | Block BFGS | OK | 11 | 1.634 | 0.164 | 1.470 | 7,994 | 9.086529004 | 2.220e-16 |
| exothermic-reactor | L-BFGS | OK | 19 | 4.927 | 0.276 | 4.651 | 13,162 | 9.086529004 | 2.220e-16 |
| moon-lander | Exact | OK | 7 | 1.514 | 0.120 | 1.394 | 5,977 | -0.058215619 | 1.110e-16 |
| moon-lander | Block BFGS | OK | 25 | 5.483 | 0.323 | 5.160 | 17,038 | -0.058215619 | 3.331e-16 |
| moon-lander | L-BFGS | OK | 10 | 3.690 | 0.140 | 3.549 | 7,348 | -0.058215619 | 1.155e-14 |
| minimum-time-moon-lander | Exact | OK | 98 | 16.021 | 0.845 | 15.176 | 37,881 | 2.007676910 | 9.072e-13 |
| minimum-time-moon-lander | Block BFGS | OK | 177 | 18.199 | 1.192 | 17.007 | 58,313 | 2.010095070 | 1.136e-09 |
| minimum-time-moon-lander | L-BFGS | OK | 170 | 47.196 | 1.571 | 45.625 | 80,970 | 2.007676832 | 3.551e-09 |
| moon-lander-obstacle | Exact | OK | 17 | 3.973 | 0.275 | 3.698 | 13,247 | -0.045674933 | 6.173e-09 |
| moon-lander-obstacle | Block BFGS | OK | 42 | 8.632 | 0.557 | 8.075 | 28,020 | -0.045674937 | 8.612e-09 |
| moon-lander-obstacle | L-BFGS | OK | 18 | 6.506 | 0.247 | 6.259 | 12,516 | -0.045674933 | 6.173e-09 |
| cart-pole-swing-up | Exact | OK | 27 | 4.146 | 0.512 | 3.634 | 15,680 | 113.624036098 | 1.642e-08 |
| cart-pole-swing-up | Block BFGS | Line search failed | 868 | 157.094 | 16.828 | 140.266 | 566,204 | 113.630678731 | 1.569e-12 |
| cart-pole-swing-up | L-BFGS | OK | 47 | 12.315 | 0.781 | 11.535 | 25,211 | 113.624036076 | 1.870e-08 |
| minimum-time-moon-lander-obstacle | Exact | OK | 134 | 26.927 | 1.307 | 25.620 | 57,939 | 2.058927126 | 8.466e-09 |
| minimum-time-moon-lander-obstacle | Block BFGS | OK | 174 | 19.734 | 1.217 | 18.517 | 57,824 | 2.063755315 | 8.849e-09 |
| minimum-time-moon-lander-obstacle | L-BFGS | OK | 236 | 74.935 | 2.669 | 72.266 | 130,685 | 2.056337998 | 8.858e-09 |
| quadrotor-transfer | Exact | OK | 45 | 14.828 | 0.966 | 13.862 | 33,926 | 81.025230435 | 9.549e-09 |
| quadrotor-transfer | Block BFGS | OK | 77 | 17.658 | 1.338 | 16.320 | 50,630 | 81.025230331 | 1.019e-08 |
| quadrotor-transfer | L-BFGS | OK | 41 | 20.382 | 0.730 | 19.652 | 28,020 | 81.025230435 | 9.549e-09 |
| quadrotor-obstacle | Exact | OK | 37 | 12.852 | 0.812 | 12.040 | 28,110 | 83.843049108 | 9.881e-09 |
| quadrotor-obstacle | Block BFGS | OK | 85 | 21.063 | 1.527 | 19.535 | 55,798 | 83.843048985 | 1.043e-08 |
| quadrotor-obstacle | L-BFGS | OK | 46 | 22.996 | 0.853 | 22.143 | 30,927 | 83.843048985 | 1.043e-08 |
| minimum-time-quadrotor | Exact | OK | 97 | 36.285 | 2.449 | 33.836 | 84,650 | 5.078520055 | 6.074e-09 |
| minimum-time-quadrotor | Block BFGS | OK | 791 | 204.077 | 13.680 | 190.397 | 514,135 | 5.078648761 | 6.106e-09 |
| minimum-time-quadrotor | L-BFGS | OK | 292 | 181.147 | 6.416 | 174.731 | 243,461 | 5.078619306 | 6.080e-09 |
| minimum-time-quadrotor-obstacle | Exact | OK | 59 | 22.811 | 1.452 | 21.359 | 50,564 | 5.152004885 | 6.313e-09 |
| minimum-time-quadrotor-obstacle | Block BFGS | Iteration cap | 1001 | 281.768 | 17.547 | 264.220 | 648,180 | 5.178596762 | 4.752e-12 |
| minimum-time-quadrotor-obstacle | L-BFGS | OK | 255 | 171.763 | 6.241 | 165.522 | 238,616 | 5.152153910 | 8.642e-09 |
| quadrotor-hourglass | Exact | OK | 52 | 17.427 | 1.126 | 16.301 | 39,338 | 82.819969270 | 9.969e-09 |
| quadrotor-hourglass | Block BFGS | OK | 99 | 24.710 | 2.075 | 22.635 | 79,054 | 82.819969270 | 9.969e-09 |
| quadrotor-hourglass | L-BFGS | OK | 46 | 22.705 | 0.805 | 21.900 | 30,604 | 82.819969358 | 9.914e-09 |
| tethered-payload | Exact | OK | 12 | 2.800 | 0.207 | 2.593 | 9,612 | 20.093477131 | 1.467e-08 |
| tethered-payload | Block BFGS | OK | 14 | 2.527 | 0.227 | 2.300 | 9,932 | 20.093477126 | 1.807e-08 |
| tethered-payload | L-BFGS | OK | 16 | 5.284 | 0.230 | 5.054 | 11,224 | 20.093477126 | 1.807e-08 |
| minimum-time-tethered-payload | Exact | OK | 111 | 32.682 | 2.016 | 30.665 | 93,113 | 4.275589059 | 6.807e-09 |
| minimum-time-tethered-payload | Block BFGS | OK | 329 | 55.912 | 4.409 | 51.503 | 213,745 | 4.271539388 | 7.325e-09 |
| minimum-time-tethered-payload | L-BFGS | OK | 79 | 21.618 | 1.126 | 20.492 | 54,183 | 4.272119334 | 2.691e-09 |
| goddard-rocket | Exact | OK | 30 | 6.622 | 0.505 | 6.117 | 22,698 | -7.135605085 | 9.937e-09 |
| goddard-rocket | Block BFGS | OK | 84 | 9.888 | 1.201 | 8.687 | 56,444 | -7.135604956 | 9.827e-09 |
| goddard-rocket | L-BFGS | OK | 138 | 17.216 | 1.922 | 15.294 | 90,036 | -7.135604949 | 9.827e-09 |
| goddard-max-q | Exact | OK | 25 | 5.471 | 0.469 | 5.003 | 19,063 | -6.359206911 | 9.766e-09 |
| goddard-max-q | Block BFGS | OK | 43 | 5.608 | 0.663 | 4.945 | 29,635 | -6.359206911 | 9.766e-09 |
| goddard-max-q | L-BFGS | OK | 151 | 18.783 | 2.320 | 16.463 | 98,434 | -6.359206911 | 9.766e-09 |
| minimum-time-bicycle-obstacle | Exact | OK | 294 | 72.903 | 5.336 | 67.567 | 179,511 | 3.804213012 | 9.911e-09 |
| minimum-time-bicycle-obstacle | Block BFGS | OK | 864 | 132.685 | 11.366 | 121.319 | 434,666 | 3.696793759 | 9.960e-09 |
| minimum-time-bicycle-obstacle | L-BFGS | OK | 471 | 192.113 | 9.446 | 182.668 | 379,748 | 3.801665086 | 9.968e-09 |
| minimum-time-quadrotor-hourglass | Exact | OK | 54 | 20.380 | 1.367 | 19.013 | 45,960 | 5.248729087 | 9.525e-09 |
| minimum-time-quadrotor-hourglass | Block BFGS | OK | 486 | 115.999 | 8.652 | 107.347 | 317,105 | 5.250590748 | 9.543e-09 |
| minimum-time-quadrotor-hourglass | L-BFGS | OK | 171 | 108.356 | 3.449 | 104.907 | 128,150 | 5.251214954 | 9.574e-09 |
| race-car-speed-limit | Exact | OK | 40 | 7.465 | 1.227 | 6.238 | 55,563 | 0.078563504 | 1.110e-16 |
| race-car-speed-limit | Block BFGS | OK | 86 | 14.296 | 1.850 | 12.445 | 85,592 | 0.078563502 | 1.110e-16 |
| race-car-speed-limit | L-BFGS | OK | 39 | 9.877 | 1.173 | 8.704 | 54,261 | 0.078563503 | 1.110e-16 |
| minimum-time-race-car | Exact | OK | 47 | 9.560 | 0.772 | 8.788 | 31,237 | 2.146834312 | 6.698e-09 |
| minimum-time-race-car | Block BFGS | OK | 407 | 49.375 | 4.681 | 44.694 | 199,442 | 2.146828059 | 9.179e-10 |
| minimum-time-race-car | L-BFGS | OK | 122 | 28.388 | 1.513 | 26.875 | 65,792 | 2.145490054 | 7.270e-09 |

## Evaluation share of total solve time

Sums across all 30 problems, including failed solves.

| Method | Evaluation (ms) | Total (ms) | Evaluation share |
|---|---:|---:|---:|
| Exact | 24.846 | 350.249 | 7.1% |
| Block BFGS | 100.524 | 1284.556 | 7.8% |
| L-BFGS | 46.284 | 1055.195 | 4.4% |

## Inspiration

Models draw on [CasADi’s race-car example](https://github.com/casadi/casadi/blob/main/docs/examples/python/race_car.py), [acados’s cart-pole example](https://github.com/acados/acados/blob/main/examples/acados_matlab_octave/getting_started/get_pendulum_on_cart_model.m), [acados’s hourglass example](https://github.com/acados/acados/blob/main/examples/acados_python/unconstrained_ocps/hour_glass_p2p_motion/hour_glass_time_optimal_p2p_motion.py), and [CasADi’s rocket example](https://github.com/casadi/casadi/blob/main/docs/examples/python/rocket.py). Equations and parameters are specific to this suite.
