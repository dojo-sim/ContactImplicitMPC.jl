# ContactControl.jl
[![CI](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml)
## Features
- [x] Contact dynamics computation
- [x] Differentiable contact simulator
- [x] MPC contact controller

## Improvements
### Algorithm
- [ ] Duals initialization and resetting strategy
- [ ] Regularization of the friction forces in the simulator
- [x] Improve simulator interior point solver (predictor-corrector method)

### Experiments
- [x] baseline comparison Raibert hopper: jump over large obstacle, test anything that is not periodic hop.
- [ ] baseline comparison with Robin Deits explicit MPC: inverted pendulum with walls, cartpole with walls,
      https://github.com/TobiaMarcucci/pympc/blob/humanoids2017/inverted_pendulum_with_wall.ipynb

### API
- [ ] Decide over the use of SVector vs SizedVector
- [ ] Ensure that the symbolic-generated functions are allocation free and fast

### Visualization
- [ ] Add friction cone visualization
- [x] Add contact force vectors visualization for 3D systems
- [x] Fix visualization of the trajectories
- [x] Fix quadruped mesh visualization

## Timings:
- [x] quadruped faster than real-time
- [x] biped faster than real-time
- [x] hopper2D faster than real-time
- [x] hopper3D faster than real-time
- [x] pushbot faster than real-time

## Examples
