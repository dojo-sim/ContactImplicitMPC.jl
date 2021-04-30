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
- [ ] Improve simulator interior point solver (predictor-corrector method)

### API
- [ ] Decide over the use of SVector vs SizedVector
- [ ] Ensure that the symbolic-generated functions are allocation free and fast

### Visualization
- [ ] Add friction cone visualization 
- [ ] Add contact force vectors visualization for 3D systems
- [ ] fix visualization of the trajectories

## Timings:
- [x] quadruped 0.5x slower
- [x] biped 1.1x slower
- [x] hopper2D 0.2x slower
- [x] hopper3D 1.4x slower
- [x] pushbot 0.5x slower
