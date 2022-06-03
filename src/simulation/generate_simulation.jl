################################################################################
# Particle (flat LC)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle")
dir_sim   = joinpath(module_dir(), "src/simulation/particle")
dir_sim_env   = joinpath(dir_sim, "flat_lc")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(particle)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat_lc/residual.jld2")
path_jac = joinpath(dir_sim, "flat_lc/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle (flat + nonlinear cone)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle")
dir_sim   = joinpath(module_dir(), "src/simulation/particle")
dir_sim_env   = joinpath(dir_sim, "flat_nc")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(particle)
env = deepcopy(flat_3D_nc)
s = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat_nc/residual.jld2")
path_jac = joinpath(dir_sim, "flat_nc/jacobians.jld2")

instantiate_base!(s.model, path_base)
instantiate_dynamics!(s.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(s.model, s.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(s, path_res, path_jac)

################################################################################
# Particle (quadratic)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle")
dir_sim   = joinpath(module_dir(), "src/simulation/particle")
dir_sim_env   = joinpath(dir_sim, "quadratic")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle 2D (flat + LC)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/particle_2D")
dir_sim_env   = joinpath(dir_sim, "flat_lc")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(particle_2D)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat_lc/residual.jld2")
path_jac = joinpath(dir_sim, "flat_lc/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle 2D (flat + nonlinear cone)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/particle_2D")
dir_sim_env   = joinpath(dir_sim, "flat_nc")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(particle_2D)
env = deepcopy(flat_2D_nc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat_nc/residual.jld2")
path_jac = joinpath(dir_sim, "flat_nc/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Particle 2D (slope)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/particle_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/particle_2D")
dir_sim_env   = joinpath(dir_sim, "slope")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (2D)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_2D")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (2D) (sinusoidal)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_2D")
dir_sim_env   = joinpath(dir_sim, "sinusoidal")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper 2D (Piecewise)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_2D")
dir_sim_env   = joinpath(dir_sim, "piecewise")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(hopper_2D)
env = deepcopy(piecewise1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "piecewise/residual.jld2")
path_jac = joinpath(dir_sim, "piecewise/jacobians.jld2")

instantiate_base!(sim.model, path_base)
expr_dyn = generate_dynamics_expressions(sim.model, derivs=true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(sim.model, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, jacobians = :approx)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac, jacobians = :approx)

################################################################################
# Hopper (2D) (stairs)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_2D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_2D")
dir_sim_env   = joinpath(dir_sim, "stairs")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(hopper_2D)
env = deepcopy(stairs3_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "stairs/residual.jld2")
path_jac = joinpath(dir_sim, "stairs/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (3D)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_3D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_3D")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Hopper (3D) (sinusoidal)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/hopper_3D")
dir_sim   = joinpath(module_dir(), "src/simulation/hopper_3D")
dir_sim_env   = joinpath(dir_sim, "sinusoidal")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/quadruped")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped Payload
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/quadruped")
dir_sim_env   = joinpath(dir_sim, "payload")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(quadruped_payload)
env = deepcopy(flat_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics_payload/base.jld2")
path_dyn = joinpath(dir_model, "dynamics_payload/dynamics.jld2")
path_res = joinpath(dir_sim, "payload/residual.jld2")
path_jac = joinpath(dir_sim, "payload/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped Sinusoidal
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/quadruped")
dir_sim_env   = joinpath(dir_sim, "sinusoidal")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(quadruped)
env = deepcopy(sine1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "sinusoidal/residual.jld2")
path_jac = joinpath(dir_sim, "sinusoidal/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, jacobians = :full)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Quadruped (Piecewise)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/quadruped")
dir_sim_env   = joinpath(dir_sim, "piecewise")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac, jacobians = :approx)

################################################################################
# Flamingo (flat)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/flamingo")
dir_sim   = joinpath(module_dir(), "src/simulation/flamingo")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (sinusoidal)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/flamingo")
dir_sim   = joinpath(module_dir(), "src/simulation/flamingo")
dir_sim_env   = joinpath(dir_sim, "sinusoidal")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (smooth slope)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/flamingo")
dir_sim   = joinpath(module_dir(), "src/simulation/flamingo")
dir_sim_env   = joinpath(dir_sim, "slope")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Flamingo (piecewise)
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/flamingo")
dir_sim   = joinpath(module_dir(), "src/simulation/flamingo")
dir_sim_env   = joinpath(dir_sim, "piecewise")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(flamingo)
env = deepcopy(piecewise1_2D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "piecewise/residual.jld2")
path_jac = joinpath(dir_sim, "piecewise/jacobians.jld2")

instantiate_base!(sim.model, path_base)
expr_dyn = generate_dynamics_expressions(sim.model, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(sim.model, path_dyn, derivs=true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, jacobians = :approx)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac, jacobians = :approx)

################################################################################
# PushBot
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/pushbot")
dir_sim   = joinpath(module_dir(), "src/simulation/pushbot")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# WalledCartpole
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/walledcartpole")
dir_sim   = joinpath(module_dir(), "src/simulation/walledcartpole")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(walledcartpole)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Rigid body
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/rigidbody")
dir_sim   = joinpath(module_dir(), "src/simulation/rigidbody")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(rigidbody)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")
instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env, mapping = Gz_func, nv = sim.model.nq - 1)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Centroidal Quadruped
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/centroidal_quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/centroidal_quadruped")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(centroidal_quadruped)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Centroidal Quadruped Box
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/centroidal_quadruped_box")
dir_sim   = joinpath(module_dir(), "src/simulation/centroidal_quadruped_box")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(centroidal_quadruped_box)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Centroidal Quadruped Wall
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/centroidal_quadruped_wall")
dir_sim   = joinpath(module_dir(), "src/simulation/centroidal_quadruped_wall")
dir_sim_env   = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(centroidal_quadruped_wall)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)

################################################################################
# Centroidal Quadruped undamped
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/centroidal_quadruped")
dir_sim   = joinpath(module_dir(), "src/simulation/centroidal_quadruped")
dir_sim_env   = joinpath(dir_sim, "flat_undamped")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(centroidal_quadruped_undamped)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics_undamped/base.jld2")
path_dyn = joinpath(dir_model, "dynamics_undamped/dynamics.jld2")
path_res = joinpath(dir_sim, "flat_undamped/residual.jld2")
path_jac = joinpath(dir_sim, "flat_undamped/jacobians.jld2")

instantiate_base!(sim.model, path_base)
instantiate_dynamics!(sim.model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)


################################################################################
# Point Foot Quadruped
################################################################################
dir_model = joinpath(module_dir(), "src/dynamics/point_foot_quadruped")
dir_sim = joinpath(module_dir(), "src/simulation/point_foot_quadruped")
dir_sim_env = joinpath(dir_sim, "flat")
!isdir(dir_sim_env) && mkdir(dir_sim_env)
model = deepcopy(point_foot_quadruped)
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
isfile(path_jac) && rm(path_jac)
ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
instantiate_residual!(sim, path_res, path_jac)
