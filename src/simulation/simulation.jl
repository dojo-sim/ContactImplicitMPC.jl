mutable struct Simulation
	model::ContactModel
	env::Environment{<:World,<:FrictionCone}
	con::ContactMethods
	res::ResidualMethods
	rz::Any
	rθ::Any
end

function Simulation(model::ContactModel, env::Environment)
	Simulation(model, env, ContactMethods(), ResidualMethods(), zeros(0, 0), zeros(0, 0))
end

function get_simulation(model::String, env::String, sim_name::String;
		model_name = model,
		gen_base = true,
		gen_dyn = true,
		approx = false)

	#TODO: assert model exists

	dir_model = joinpath(pwd(), "src/dynamics/", model)
	dir_sim   = joinpath(pwd(), "src/simulation", model, sim_name)

	dir_base = joinpath(dir_model, "dynamics/base.jld2")
	dir_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
	dir_res = joinpath(dir_sim, "residual.jld2")
	dir_jac = joinpath(dir_sim, "jacobians.jld2")

	model = deepcopy(eval(Symbol(model)))
	env = deepcopy(eval(Symbol(env)))
	sim = Simulation(model, env)

	instantiate_base!(sim.model, dir_base)

	instantiate_dynamics!(sim.model, dir_dyn,
		derivs = approx)

	instantiate_residual!(sim,
		dir_res, dir_jac,
		jacobians = (approx ? :approx : :full))

	return sim
end
