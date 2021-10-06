# ContactImplicitMPC.jl
[![CI](https://github.com/thowell/ContactImplicitMPC.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/thowell/ContactImplicitMPC.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/thowell/ContactImplicitMPC.jl/branch/main/graph/badge.svg?token=3J4VOJ0VCH)](https://codecov.io/gh/thowell/ContactImplicitMPC.jl)

This repository contains algorithms and examples from our paper: [Fast Contact-Implicit Model-Predictive Control](https://arxiv.org/abs/2107.05616).

A collection of examples are pre-generated in notebooks with the package, please try: [flamingo](examples/flamingo/flat.jl), [pushbot](examples/pushbot/push_recovery.jl), [hopper](examples/hopper/flat.jl), and [quadruped](examples/quadruped/flat.jl). Additional notebooks with examples from the paper can be [generated](examples/README.md).

## Installation
- `ContactImplicitMPC` can be added via the Julia package manager (type `]`): `pkg> add ContactImplicitMPC#main` (this will grab the latest version of the package from this repo)

## Flamingo
<img src="animations/flamingo.gif" alt="drawing" width="400"/>

## PushBot
<img src="animations/pushbot.gif" alt="drawing" width="400"/>

## Hopper Parkour
<img src="animations/hopper_parkour.gif" alt="drawing" width="250"/>

## Quadruped with Payload
<img src="animations/quadruped_payload.gif" alt="drawing" width="400"/>

## Hopper Monte Carlo
<img src="animations/hopper_monte_carlo.gif" alt="drawing" width="400"/>

## Quadruped Monte Carlo
<img src="animations/quadruped_monte_carlo.gif" alt="drawing" width="400"/>

## Walled Cart-Pole
<img src="animations/walledcartpole.gif" alt="drawing" width="400"/>

## Reference Trajectories
The trajectories we track in the examples are generated using [contact-implicit trajectory optimization](https://journals.sagepub.com/doi/10.1177/0278364919849235) and can be run [here](https://github.com/thowell/motion_planning/tree/main/examples/contact_implicit).
