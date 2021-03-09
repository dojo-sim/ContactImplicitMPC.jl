
function easy_lin_step(lin::LinStep14, q1_ref, q2_ref,
	model::ContactDynamicsModel; u2_ref=zeros(SVector{model.dim.u,Float64}),
	ρ0=1e0, outer_iter=5, inner_iter=100, tol=1e-8, step_print::Bool=false, z_init=3e-2)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	nz = nq + 4nγ + 4nb

	function rz(z::AbstractVector)
		θ = [q1_ref; q2_ref; u2_ref]
		return triv_lin_rzθ(model, lin, z, θ, ρ)
	end

	# Jacobian
	function Rz(z::AbstractVector)
		# differentiate r
		return easy_Rz(model, lin.bil_terms, lin.bil_vars, lin.Rz0, z)
	end

    # initialize
	z = z_init * (ones(eltype(q1_ref), nz) + ones(eltype(q2_ref), nz) +  ones(eltype(u2_ref), nz))
	z[1:nq] = copy(q2_ref)

    ρ = ρ0 # barrier parameter
    flag = false

    for k = 1:outer_iter
        for i = 1:inner_iter
            # compute residual, residual Jacobian
            res = rz(z)
            if norm(res) < tol
                step_print ? println("   iter ($i) - norm: $(scn(norm(res)))") : nothing
                # return z, true
                flag = true
                # continue
                break
            end

            jac = Rz(z)
            Δ = jac \ res

            # line search the step direction
            α = 1.0
            iter = 0
            while check_variables(model, z - α * Δ) # backtrack inequalities
                α = 0.5 * α
                iter += 1
                if iter > 50
                    @error "backtracking line search fail"
                    flag = false
                    return z, false
                end
            end

            while norm(rz(z - α * Δ))^2.0 >= (1.0 - 0.001 * α) * norm(res)^2.0
                α = 0.5 * α
                # println("   α = $α")
                iter += 1
                if iter > 50
                    @error "line search fail"
                    flag = false
                    return z, false
                end
            end

			# update
            z = z - α * Δ
        end

        ρ = 0.1 * ρ
        step_print ? println("ρ: $(scn(ρ))") : nothing
    end

	_Rz = easy_Rz(model, lin.bil_terms, lin.bil_vars, lin.Rz0, z)
	_Rθ = lin.Rθ0
	∇ = easy_∇(model, lin.bil_terms, lin.bil_vars, lin.Rz0, lin.Rθ0, z)
    return z, ∇, flag, _Rz, _Rθ
end






function easy_lin_step!(sim, t)
    # unpack
    model = sim.model
    q = sim.q
    u = sim.u
    w = sim.w
    h = sim.h
    ip = sim.ip
    z = ip.z
    θ = ip.θ

    # initialize
    if sim.sim_opts.warmstart
        z .+= sim.sim_opts.z_warmstart * rand(ip.num_var)
        sim.ip_opts.κ_init = sim.sim_opts.κ_warmstart
    else
        z_initialize!(z, model, q[t+1])
    end
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t], h)

    # solve
    status = interior_point!(ip, opts = sim.ip_opts)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(model, z)
        sim.q[t+2] = copy(q2)
        sim.γ[t] = γ
        sim.b[t] = b

        if sim.ip_opts.diff_sol
            nq = model.dim.q
            nu = model.dim.u
            nc = model.dim.c
            nb = model.dim.b

            sim.dq2dq0[t] = view(ip_data.δz, 1:nq, 1:nq)
            sim.dq2dq1[t] = view(ip_data.δz, 1:nq, nq .+ (1:nq))
            sim.dq2du[t] = view(ip_data.δz, 1:nq, 2 * nq .+ (1:nu))
            sim.dγdq0[t] = view(ip_data.δz, nq .+ (1:nc), 1:nq)
            sim.dγdq1[t] = view(ip_data.δz, nq .+ (1:nc), nq .+ (1:nq))
            sim.dγdu[t] = view(ip_data.δz, nq .+ (1:nc), 2 * nq .+ (1:nu))
            sim.dbdq0[t] = view(ip_data.δz, nq + nc .+ (1:nb), 1:nq)
            sim.dbdq1[t] = view(ip_data.δz, nq + nc .+ (1:nb), nq .+ (1:nq))
            sim.dbdu[t] = view(ip_data.δz, nq + nc .+ (1:nb), 2 * nq .+ (1:nu))
        end
    end
    return status
end





@with_kw struct SimulatorOptions{T}
    warmstart::Bool = true
    z_warmstart::T = 0.001
    κ_warmstart::T = 0.001
end

struct Simulator{S,nq,nu,nc,nb,nw}
    model

    T::Int
    h::S

    q::Vector{SArray{Tuple{nq},S,1,nq}}
    u::Vector{SArray{Tuple{nu},S,1,nu}}
    γ::Vector{SArray{Tuple{nc},S,1,nc}}
    b::Vector{SArray{Tuple{nb},S,1,nb}}
    w::Vector{SArray{Tuple{nw},S,1,nw}}

    dq2dq0::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2dq1::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2du::Vector{SizedArray{Tuple{nq,nu},S,2,2}}
    dγdq0::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdq1::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdu::Vector{SizedArray{Tuple{nc,nu},S,2,2}}
    dbdq0::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdq1::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdu::Vector{SizedArray{Tuple{nb,nu},S,2,2}}

    ip::InteriorPoint{S}
    ip_opts::InteriorPointOptions{S}

    sim_opts::SimulatorOptions{S}
end

function simulator(model, q0::SVector, q1::SVector, h::S, T::Int;
    u = [@SVector zeros(model.dim.u) for t = 1:T],
    w = [@SVector zeros(model.dim.w) for t = 1:T],
    ip_opts = InteriorPointOptions{S}(),
    r! = r!, rz! = rz!, rθ! = rθ!,
    rz = spzeros(num_var(model), num_var(model)),
    rθ = spzeros(num_var(model), num_data(model)),
    sim_opts = SimulatorOptions{S}()) where S

    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b

    q = [q0, q1, [@SVector zeros(nq) for t = 1:T]...]
    γ = [@SVector zeros(nc) for t = 1:T]
    b = [@SVector zeros(nb) for t = 1:T]

    dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:T]
    dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:T]
    dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
    dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
    dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:T]

    ip = interior_point(
        num_var(model),
        num_data(model),
        inequality_indices(model),
        r! = r!, rz! = rz!, rθ! = rθ!,
        rz = rz,
        rθ = rθ)

    Simulator(
        model,
        T, h,
        q,
        u,
        γ,
        b,
        w,
        dq2dq0,
        dq2dq1,
        dq2du,
        dγdq0,
        dγdq1,
        dγdu,
        dbdq0,
        dbdq1,
        dbdu,
        ip,
        ip_opts,
        sim_opts)
end






T = Float64
model = get_model("quadruped")
model = get_model("particle")
nq = model.dim.q
nc = model.dim.c
nu = model.dim.u
nw = model.dim.w
nz = num_var(model)
nθ = num_data(model)

ip_opts = InteriorPointOptions(κ_init=1e-3, κ_tol=1e-2)
ip = interior_point(
	num_var(model),
	num_data(model),
	inequality_indices(model),
	r! = model.res.r,
	rz! = model.res.rz,
	rθ! = model.res.rθ,
	rz = model.spa.rz_sp,
	rθ = model.spa.rθ_sp) # not correct

q0 = zeros(SizedVector{nq,T})
q1 = zeros(SizedVector{nq,T})
u1 = rand(SizedVector{nu,T})
w1 = rand(SizedVector{nw,T})
h = 0.01
θ_initialize!(ip.θ, model, q0, q1, u1, w1, h)
z_initialize!(ip.z, model, q1)

@btime status = interior_point!(ip, opts = ip_opts)

z0 = deepcopy(ip.z)
θ0 = deepcopy(ip.θ)
κ0 = deepcopy(ip.κ)
lin = LinStep14(model, z0, θ0, κ0[1])

# residual
function r!(r, z, θ, κ)
    @warn "approx"
	r_approx!(lin, r, z, θ, κ)
end

# residual Jacobian wrt z
function rz!(rz, z, θ, κ)
	@warn "approx"
	rz_approx!(lin, rz, z, θ, κ)
end

# residual Jacobian wrt θ
function rθ!(rθ, z, θ, κ)
	@warn "approx"
	rθ_approx!(lin, rθ, z, θ, κ)
	nothing
end

ip_approx = interior_point(
	num_var(model),
	num_data(model),
	inequality_indices(model),
	r! = r!, rz! = rz!, rθ! = rθ!,
	rz = model.spa.rz_sp,
	rθ = model.spa.rz_sp) # not correct

θ_initialize!(ip_approx.θ, model, q0, q1, u1, w1, h)
z_initialize!(ip_approx.z, model, q1)

status = interior_point!(ip_approx, opts = ip_opts)
@test norm(ip_approx.z - z0, Inf) < 1e-8
@test norm(ip_approx.θ - θ0, Inf) < 1e-8
@test norm(ip_approx.r, Inf) < 1e-5

function lin_step(lin::LinStep14, )


function inner(b; lin=100, pol=200)
    return lin+pol*b
end

function outer(a, b; kwargs...)
    @show kwargs
    return a + inner(b; kwargs...)
    return nothing
end

outer(10,10; lin=100, pol=300)



mutable struct ControlTraj12{T,nq,nu,nw,nc,nb}
	N::Int                       # horizon length
	lin::Vector{LinStep14{T}}    # linearization point length=N
	q::Vector{SizedVector{nq,T}} # trajectory of q's   length=N+2
	u::Vector{SizedVector{nu,T}} # trajectory of u's   length=N
	w::Vector{SizedVector{nw,T}} # trajectory of w's   length=N
	γ::Vector{SizedVector{nc,T}} # trajectory of γ's   length=N
	b::Vector{SizedVector{nb,T}} # trajectory of b's   length=N
	d::Vector{SizedVector{nq,T}} # dynamics violation  length=N
	δz::SparseMatrixCSC{T,Int}   # solution gradients  length=N
end


function ControlTraj12(lin::Vector{LinStep14{T}}, q::Vector{SizedVector{nq,T}},
	u::Vector{SizedVector{nu,T}}, w::Vector{SizedVector{nw,T}},
	γ::Vector{SizedVector{nc,T}}, b::Vector{SizedVector{nb,T}}) where {T,nq,nu,nw,nc,nb}
	N = length(q)-2
	@assert length(u) == length(w) == length(γ) == length(b) == N
	nz = length(lin[1].z0)
	nθ = length(lin[1].θ0)

	d = [zeros(SizedVector{nq,T}) for k=1:T]
	δz = [spzeros(nz,nθ) for k=1:T]
	return ControlTraj12{T,nq,nu,nw,nc,nb}(lin,q,u,w,γ,b,d,δz)
end

N = 10
lin = [ for k = 1:N]
ControlTraj12(lin, q, u, w, γ, b)
