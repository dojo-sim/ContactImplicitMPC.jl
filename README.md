# ContactControl.jl
[![CI](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml)

- Dynamics computation
- Contact simulator
- Contact-implicit controller

- [ ] Decide on the use of SVector vs SizedVector
- [ ] design quadruped gait so that the four feet are in contact with the ground for a few times step between each step, instead of having a 'trotting' gait
- [ ] Enforce no slip conditions on nominal gaits
- [ ] Improve visualizations: visualize the torques at the revolute joints and the forces at the prismatic joints
- [ ] Add friction cone visualization and contact force vectors
- [ ] fix visualization of the trajectories
- [ ] add terrain visualization

#Timings:
-  quadruped 1.5x slower
-  hopper2D 0.2x slower
-  
