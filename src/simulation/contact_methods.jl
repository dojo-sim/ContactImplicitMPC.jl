mutable struct ContactMethods
	J
	dJ
	ϕ
	d
	dλ1
	cf
	dcf
	vs
	vsq2
	vsq1h
	mdvs
	mdψη
	rc
	rcz
	rcθ
end

function ContactMethods()
	function f()
		error("Not Implemented: use instantiate_contact_methods!")
		return nothing
	end
	return ContactMethods(fill(f, 15)...)
end

function contact_forces(model::ContactModel, env::Environment{<:World,LinearizedCone}, γ1, b1, q2, k)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nf = Int(nb / nc)
	ne = dim(env)
	λ1 = vcat([transpose(rotation(env, k[(i-1) * (ne - 1) .+ (1:ne)])) * [friction_mapping(env) * b1[(i-1) * nf .+ (1:nf)]; γ1[i]] for i = 1:nc]...) # TODO: make efficient
end

function contact_forces(model::ContactModel, env::Environment{<:World,NonlinearCone}, γ1, b1, q2, k)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	ne = dim(env)
	λ1 = vcat([transpose(rotation(env, k[(i-1) * (ne - 1) .+ (1:ne)])) * [b1[(i-1) * (ne-1) .+ (1:(ne-1))]; γ1[i]] for i = 1:nc]...) # TODO: make efficient
end

function velocity_stack(model::ContactModel, env::Environment{<:World,LinearizedCone}, q1, q2, k, h)
	nc = model.dim.c
	ne = dim(env)
	v = J_fast(model, q2) * (q2 - q1) / h[1]
	v_surf = [rotation(env, k[(i-1) * (ne - 1) .+ (1:ne)]) * v[(i-1) * ne .+ (1:ne)] for i = 1:nc]
	vT_stack = vcat([[v_surf[i][1:ne-1]; -v_surf[i][1:ne-1]] for i = 1:nc]...)
end

function velocity_stack(model::ContactModel, env::Environment{<:World,NonlinearCone}, q1, q2, k, h)
	nc = model.dim.c
	ne = dim(env)
	v = J_fast(model, q2) * (q2 - q1) / h[1]
	v_surf = [rotation(env, k[(i-1) * ne .+ (1:ne)]) * v[(i-1) * ne .+ (1:ne)] for i = 1:nc]
	vT_stack = vcat([v_surf[i][1:ne-1] for i = 1:nc]...)
end
