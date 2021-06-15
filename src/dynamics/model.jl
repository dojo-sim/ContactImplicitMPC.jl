abstract type ContactModel end

struct Dimensions
    q::Int         # configuration
    u::Int         # control
	w::Int         # disturbance
	c::Int         # contact points
end

# https://github.com/HarvardAgileRoboticsLab/drake/blob/75b260c9eb250d08ffbbf3daa80758e4fe558d7f/drake/matlab/solvers/trajectoryOptimization/VariationalTrajectoryOptimization.m
function lagrangian_derivatives(model::ContactModel, q, v)
	D1L = -1.0 * C_fast(model, q, v)
    D2L = M_fast(model, q) * v
	return D1L, D2L
end

# function dynamics(model::ContactModel, h, q0, q1, u1, w1, λ1, q2)
# 	# evalutate at midpoint
# 	qm1 = 0.5 * (q0 + q1)
#     vm1 = (q1 - q0) / h[1]
#     qm2 = 0.5 * (q1 + q2)
#     vm2 = (q2 - q1) / h[1]
#
# 	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
# 	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)
#
# 	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
# 		+ transpose(B_fast(model, qm2)) * u1
# 		+ transpose(A_fast(model, qm2)) * w1
# 		+ transpose(J_fast(model, q2)) * λ1
# 		- h[1] * model.joint_friction .* vm2)
# end

function dynamics(model::ContactModel, h, q0, q1, u1, w1, Λ1, q2)
	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	# return 0.0
	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		# + transpose(J_fast(model, q2)) * λ1
		+ Λ1
		- h[1] * model.joint_friction .* vm2)

end

function dynamics(model::ContactModel, env::Environment, h, q0, q1, u1, w1, λ1, q2)
	Λ1 = transpose(J_func(model, env, q2)) * λ1 #@@@@ maybe need to use J_fast
	dynamics(model, h, q0, q1, u1, w1, Λ1, q2)
end

mutable struct BaseMethods
	L::Any
	M::Any
	B::Any
	A::Any
	# J::Any
	C::Any
	k::Any
end

function BaseMethods()
	function f()
		error("Not Implemented: use instantiate_base!")
		return nothing
	end
	return BaseMethods(fill(f, 6)...)
end

mutable struct DynamicsMethods
	d::Any
	dy::Any
	dθ::Any
	dq0::Any
	dq1::Any
	du1::Any
	dw1::Any
	dγ1::Any
	db1::Any
	dΛ1::Any
	∂q2::Any
end

function DynamicsMethods()
	function f()
		error("Not Implemented: use instantiate_dynamics!")
		return nothing
	end
	return DynamicsMethods(fill(f, 11)...)
end

"""
	get_model(name::String, surf::String)
	Helper function that provides a model where fast functions have been instantiated.
"""
function get_model(name::String;
	model_name::String = name, approx = false)

	#TODO: assert model exists
	path = joinpath(@__DIR__, name)
	model = eval(Symbol(model_name * (surf != "flat" ? "_" * surf : "")))
	instantiate_base!(model, joinpath(path, dynamics, "base.jld2"))
	instantiate_dynamics!(model, joinpath(path, dynamics, "dynamics.jld2"),
		derivs = approx)

	return model
end

function get_gait(name::String, gait::String)
	#TODO: assert model exists
	path = joinpath(@__DIR__, name)
	gait_path = joinpath(path, "gaits/" * gait * ".jld2")

	res = JLD2.jldopen(gait_path) # z̄ x̄ ū h̄ q u γ b

	return res["q"], res["u"], res["γ"], res["b"], mean(res["h̄"])
end

function model_name(model::ContactModel)
    name = Symbol(string(typeof(model).name)[10:end-1])
    return name
end
