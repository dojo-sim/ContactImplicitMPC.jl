# ContactControl.jl
[![CI](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/simon-lc/ContactControl.jl/actions/workflows/CI.yml)

- Dynamics computation
- Contact simulator
- Contact-implicit controller

- [ ] Decide on the use of SVector vs SizedVector
- [ ] design quadruped gait so that the four feet are in contact with the ground for a few times step between each step, instead of having a 'trotting' gait
- [ ] Enforce no slip conditions on nominal gaits
- [ ] 

/!\ Changed biped friction coefficient to be 0.7
/!\ Some tests on the fast dynamics are not passing 

#Timings:
-  quadruped 20x slower
-  hopper2D 2x slower
-  
