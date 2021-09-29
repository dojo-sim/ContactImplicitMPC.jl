# ContactImplicitMPC.jl
[![CI](https://github.com/thowell/ContactImplicitMPC.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/thowell/ContactImplicitMPC.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/thowell/ContactImplicitMPC.jl/branch/main/graph/badge.svg?token=3J4VOJ0VCH)](https://codecov.io/gh/thowell/ContactImplicitMPC.jl)

This repository contains the algorithms and examples from our paper: [Fast Contact-Implicit Model-Predictive Control](https://arxiv.org/abs/2107.05616).

We are currently cleaning up the examples for public consumption. At this stage, please try: [flamingo](examples/flamingo_flat.jl), [pushbot](examples/pushbot.jl), [hopper](examples/hopper_flat.jl), and [quadruped](examples/quadruped.jl).

## Installation
```
Pkg.add("ContactImplicitMPC")
```

## Flamingo
<img src="examples/animations/flamingo.gif" alt="drawing" width="400"/>

## PushBot
<img src="examples/animations/pushbot.gif" alt="drawing" width="400"/>

## Hopper Parkour
<img src="examples/animations/hopper_parkour.gif" alt="drawing" width="250"/>

## Quadruped with Payload
<img src="examples/animations/quadruped_payload.gif" alt="drawing" width="400"/>

## Hopper Monte Carlo
<img src="examples/animations/hopper_monte_carlo.gif" alt="drawing" width="400"/>

## Quadruped Monte Carlo
<img src="examples/animations/quadruped_monte_carlo.gif" alt="drawing" width="400"/>
