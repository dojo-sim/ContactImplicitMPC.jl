
include("newton.jl")



# vis = Visualizer()
# open(vis)

opts.r_tol = 1e-6
core1 = Newton(H, h, model)
linearization!(model, ref_traj0, impl0)
@time newton_solve!(model, core1, impl0, cost0, ref_traj0, opts, initial_offset=false)
@time newton_solve!(model, core1, impl0, cost0, ref_traj0, opts, warm_start=true, initial_offset=true)
















T = Float64
κ = 1e-4
# model = get_model("quadruped")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b

# time
h = h̄
H = length(u)
# H = 15

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
impl0 = ImplicitTraj(H, model)

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

sim0 = simulator(model, q0, q1, h, H,
	p = open_loop_policy([SVector{model.dim.u}(h * u[i]) for i=1:H], h),
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
    sim_opts = SimulatorOptions(warmstart = true))

simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
traj0 = deepcopy(sim0.traj)

norm(q[1]   - traj0.q[1])/norm(traj0.q[1])
norm(q[2]   - traj0.q[2])/norm(traj0.q[2])
norm(q[3]   - traj0.q[3])/norm(traj0.q[3])
norm(h*u[1] - traj0.u[1])/norm(traj0.u[1])
norm(h*γ[1] - traj0.γ[1])/norm(traj0.γ[1])
norm(h*b[1] - traj0.b[1])/norm(traj0.b[1])
norm(h*b[2] - traj0.b[2])/norm(traj0.b[2])
norm(h*b[3] - traj0.b[3])/norm(traj0.b[3])


nz = num_var(model)
r0 = zeros(nz)
model.res.r(r0, traj0.z[1], traj0.θ[1], κ)
@test norm(r0) < 1e-8

# vis = Visualizer()
# open(vis)
visualize!(vis, model, traj0.q)

linearization!(model, ref_traj0, impl0)
@time implicit_dynamics!(model, ref_traj0, impl0, κ=κ)
@test mean([norm(d) for d in impl0.d]) < 1e-6
mean([norm(d) for d in impl0.d])
plot([norm(d,Inf) for d in impl0.d])







δz_ = [impl0.δq0[1] impl0.δq1[1] impl0.δu1[1]]
@test norm(impl0.δz[1][1:size(δz_)[1], 1:size(δz_)[2]] - δz_) == 0.0


cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-2*SizedVector{nq}([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,])), H),
    Qu=fill(Diagonal(3e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-6*ones(SizedVector{nc})), H),
    Qb=fill(Diagonal(1e-6*ones(SizedVector{nb})), H),
    )
opts = NewtonOptions()

core0 = Newton(H, h, model)
core0.r
traj1 = deepcopy(traj0)
for t = 1:H
    traj1.q[t+2] .+= ones(nq)
    traj1.u[t] .+= ones(nu)
    traj1.w[t] .+= ones(nw)
    traj1.γ[t] .+= ones(nc)
    traj1.b[t] .+= ones(nb)
end
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, traj0, ref_traj0, opts)
norm(core0.r.r) == 0.0
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, traj1, ref_traj0, opts)
norm(core0.r.r) == 0.0
#
# off = 0
# core0.r.r[off .+ (1:nq)]
# off += nq
# core0.r.r[off .+ (1:nu)]
# off += nu
# core0.r.r[off .+ (1:nc)]
# off += nc
# core0.r.r[off .+ (1:nb)]
# off += nb
# core0.r.r[off:end]
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, ref_traj0, ref_traj0, opts)
# @allocated residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, opts)
# @code_warntype residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, opts)
# @benchmark residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, opts)

@time jacobian!(model, core0, core0.j, impl0, cost0, opts)
@time jacobian!(model, core0, core0.j, impl0, cost0, opts)
# @allocated jacobian!(model, core0, core0.j, impl0, cost0, opts)
# @code_warntype jacobian!(model, core0, core0.j, impl0, cost0, opts)
# @benchmark jacobian!(model, core0, core0.j, impl0, cost0, opts)

# core1 = Newton(H, h, model)
# @time newton_solve!(model, core1, impl0, cost0, ref_traj0, opts)
visualize!(vis, model, ref_traj0.q, Δt=5*h)

scn(norm(core0.r.r, 1)/length(core0.r.r))
