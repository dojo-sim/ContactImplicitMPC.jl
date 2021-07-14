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

# Interior Point Roadmap
## full friction cone
### models 
- [ ] create particle models with/without linear friction cone
- [x] create box with/without linear friction cone
- [ ] check with Zac how to finite difference MRPs
### algorithm
- [ ] implement ECOS-like LP solver
- [ ] implement ECOS-like QP solver
- [ ] derive Nesterov Todd scaling for contact problem
- [ ] derive Nesterov Todd scaling for conic contact problem
- [ ] derive analytical line-search for conic contact problem
- [ ] implement ECOS-like contact solver
### tuning
- [ ] create benchmark with multiple models
- [ ] test and tune heuristics from ECOS, CVXOPT, Larry Biegler's IPOPT
- [ ] ensure that condition number remains reasonable

## quaternions
### models 
- [ ] create box models with/without quaternion, with/without linear friction cone
### algorithm
- [ ] deal with non-euclidean space
### tuning
- [ ] create benchmark with multiple models
- [ ] test and tune heuristics from ECOS, CVXOPT, Larry Biegler's IPOPT
- [ ] ensure that condition number remains reasonable
