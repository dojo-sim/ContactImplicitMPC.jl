###############################################################################
# Particle (flat)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/particle")
dir_sim   = joinpath(pwd(), "src/simulation/particle")
model = deepcopy(particle)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle (quadratic)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/particle")
dir_sim   = joinpath(pwd(), "src/simulation/particle")
model = deepcopy(particle)
env = deepcopy(quadratic_bowl_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "quadratic/residual.jld2")
path_jac = joinpath(dir_sim, "quadratic/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle 2D (flat)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/particle_2D")
dir_sim   = joinpath(pwd(), "src/simulation/particle_2D")
model = deepcopy(particle_2D)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle 2D (slope)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/particle_2D")
dir_sim   = joinpath(pwd(), "src/simulation/particle_2D")
model = deepcopy(particle_2D)
env = deepcopy(slope1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "slope/residual.jld2")
path_jac = joinpath(dir_sim, "slope/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (2D)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(pwd(), "src/simulation/hopper_2D")
model = deepcopy(hopper_2D)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (2D) (sinusoidal)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(pwd(), "src/simulation/hopper_2D")
model = deepcopy(hopper_2D)
env = deepcopy(sine2_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (3D)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/hopper_3D")
dir_sim   = joinpath(pwd(), "src/simulation/hopper_3D")
model = deepcopy(hopper_3D)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (3D) (sinusoidal)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/hopper_3D")
dir_sim   = joinpath(pwd(), "src/simulation/hopper_3D")
model = deepcopy(hopper_3D)
env = deepcopy(sine2_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped")
model = deepcopy(quadruped)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped Payload
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped")
model = deepcopy(quadruped_payload)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "payload/residual.jld2")
path_jac = joinpath(dir_sim, "payload/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped Sinusoidal
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped")
model = deepcopy(quadruped)
env = deepcopy(sine1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped (Piecewise)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped")
model = deepcopy(quadruped)
env = deepcopy(piecewise1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "piecewise/residual.jld2")
path_jac = joinpath(dir_sim, "piecewise/jacobians.jld2")

instantiate_base!(sim.model, path_base)
expr_dyn = generate_dynamics_expressions(sim.model, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(sim.model, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, jacobians = :approx)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac, jacobians = :approx)

################################################################################
# Quadruped (3D)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped_3D")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped_3D")
model = deepcopy(quadruped_3D)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped (simple)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/quadruped_simple")
dir_sim   = joinpath(pwd(), "src/simulation/quadruped_simple")
model = deepcopy(quadruped_simple)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Biped
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/biped")
dir_sim   = joinpath(pwd(), "src/simulation/biped")
model = deepcopy(biped)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Biped (sinusoidal)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/biped")
dir_sim   = joinpath(pwd(), "src/simulation/biped")
model = deepcopy(biped)
env = deepcopy(sine1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (flat)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/flamingo")
dir_sim   = joinpath(pwd(), "src/simulation/flamingo")
model = deepcopy(flamingo)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (sinusoidal)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/flamingo")
dir_sim   = joinpath(pwd(), "src/simulation/flamingo")
model = deepcopy(flamingo)
env = deepcopy(sine3_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (smooth slope)
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/flamingo")
dir_sim   = joinpath(pwd(), "src/simulation/flamingo")
model = deepcopy(flamingo)
env = deepcopy(slope_smooth_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "slope/residual.jld2")
path_jac = joinpath(dir_sim, "slope/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# PushBot
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/pushbot")
dir_sim   = joinpath(pwd(), "src/simulation/pushbot")
model = deepcopy(pushbot)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# PlanarPush
################################################################################
dir_model = joinpath(pwd(), "src/dynamics/planarpush")
dir_sim   = joinpath(pwd(), "src/simulation/planarpush")
model = deepcopy(planarpush)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)
