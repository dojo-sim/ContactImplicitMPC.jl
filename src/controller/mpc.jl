function mpc_init(model,
    q1_init, q2_init, q1_ref, q2_ref, u_ref, q_ref, γ_ref, b_ref;
    Nsample::Int=3, Nmpc::Int=20, z_init=3e-2)

    N = length(u_ref)+1
    nq = model.dim.q
    nu = model.dim.u
    nγ = model.dim.γ
    nb = model.dim.b
    nr = 2nq+nu+2nγ+4nb # outer problem residual size

    q_truth = [q2_init-(q2_init-q1_init)/Nsample, q2_init]
    u_truth = []
    dist_truth = []
    γ_truth = []
    q1_init_rot = deepcopy(q1_init)
    q2_init_rot = deepcopy(q2_init)
    q1_ref_rot = deepcopy(q1_ref)
    q2_ref_rot = deepcopy(q2_ref)
    u_rot = deepcopy(u_ref)
    q_rot = deepcopy([q_ref[1:end-2]; [mpc_stride_q(model,q1_ref), mpc_stride_q(model,q2_ref)]])

    γ_rot = deepcopy(γ_ref)
    b_rot = deepcopy(b_ref)
    traj = zeros((Nmpc-1)*nr)
    set_var_traj!(Nmpc, model, traj, deepcopy(u_rot), :u)
    set_var_traj!(Nmpc, model, traj, deepcopy(q_rot), :q)
    set_var_traj!(Nmpc, model, traj, deepcopy(γ_rot), :γ)
    set_var_traj!(Nmpc, model, traj, deepcopy(b_rot), :b)

    lint0 = lintraj(q1_ref, q2_ref, deepcopy([q_ref[1:end-2]; [q1_ref, q2_ref]]),
        u_ref, model, ρ0=1e-3, outer_iter=1, z_init=z_init)

    compute_times = zeros(0)
    return Nsample, Nmpc,
        q_truth, u_truth, dist_truth, γ_truth,
        q1_init_rot, q2_init_rot, q1_ref_rot, q2_ref_rot,
        traj, u_rot, q_rot, γ_rot, b_rot, lint0, compute_times
end


function mpc_step(model, model_oa, N, Nsample, Nmpc,
    q_truth, u_truth, dist_truth, γ_truth,
    q1_init_rot, q2_init_rot, q1_ref_rot, q2_ref_rot,
    traj, u_rot, q_rot, γ_rot, b_rot, lint0, compute_times;
    model_oa_rd=deepcopy(model_oa),
    mpc_iter::Int=20,
    dist_amp::T=0.045,
    dist_list=nothing,
    ρ0::T=1e-3,
    res_tol::T=1e-5,
    outer_iter::Int=1,
    solver_outer_iter::Int=1,
    solver_inner_iter::Int=10,
    Qu=fill(Diagonal(5e-3*[0.2, 30.0]), N-1),
    Qq=fill(Diagonal(2e-2*ones(probsize.nq)), N-1),
    Qγ=fill(Diagonal(1e-100*ones(probsize.nγ)), N-1),
    Qb=fill(Diagonal(1e-100*ones(2probsize.nb)), N-1),
    sim_outer_iter::Int=3,
    sim_ρ0::T=1e-3,
    z_init=3e-2,
    control_mode::Symbol=:mpc,
    ) where T

    nq = model.dim.q
    nu = model.dim.u
    nγ = model.dim.γ
    nb = model.dim.b
    nr = 2nq+nu+2nγ+4nb # outer problem residual size
    probsize = ProblemSize(Nmpc,nq,nu,nc)

    model_hf = deepcopy(model)
    model_hf.dt /= Nsample
    model_oa_hf = deepcopy(model_oa)
    model_oa_hf.dt /= Nsample
    model_oa_rd_hf = deepcopy(model_oa_rd)
    model_oa_rd_hf.dt /= Nsample

    for l = 1:mpc_iter
        # Visualization
        plt = plot(xlims=(-2, 15))
        scatter!(-ones(1), q1_init_rot[2:2], color=:red, markersize=8.0, label="y")
        scatter!(-ones(1), q1_init_rot[4:4], color=:blue, markersize=8.0, label="l1")
        scatter!(-ones(nq), q1_init_rot)
        scatter!(zeros(nq), q2_init_rot)
        plot!(reshape(vcat(q_rot...), (nq,N-1))', label="q");# display(plt)

        # Warm Starting
        utr_wrm = deepcopy(get_var_traj(Nmpc,model,traj,:u))
        qtr_wrm = deepcopy(get_var_traj(Nmpc,model,traj,:q))
        γtr_wrm = deepcopy(get_var_traj(Nmpc,model,traj,:γ))
        btr_wrm = deepcopy(get_var_traj(Nmpc,model,traj,:b))
        rotate_seq!(utr_wrm)
        rotate_seq!(qtr_wrm)
        rotate_seq!(γtr_wrm)
        rotate_seq!(btr_wrm)

        # TrajOpt Solve
        compute_time = @elapsed traj = easy_lin_trajopt(probsize, model, q1_init_rot, q2_init_rot, q1_ref_rot, q2_ref_rot, lint0[1:Nmpc-1],
                                    ρ0=ρ0,
                                    res_tol=res_tol,
                                    outer_iter=outer_iter,
                                    solver_outer_iter=solver_outer_iter,
                                    solver_inner_iter=solver_inner_iter,
                                    utr_ref=deepcopy(u_rot)[1:Nmpc-1],
                                    qtr_ref=deepcopy(q_rot)[1:Nmpc-1],
                                    γtr_ref=deepcopy(γ_rot)[1:Nmpc-1],
                                    btr_ref=deepcopy(b_rot)[1:Nmpc-1],
                                    utr_wrm=utr_wrm[1:Nmpc-1],
                                    qtr_wrm=qtr_wrm[1:Nmpc-1],
                                    γtr_wrm=γtr_wrm[1:Nmpc-1],
                                    btr_wrm=btr_wrm[1:Nmpc-1],
                                    Qu=Qu[1:Nmpc-1],
                                    Qq=Qq[1:Nmpc-1],
                                    Qγ=Qγ[1:Nmpc-1],
                                    Qb=Qb[1:Nmpc-1],
                                    u_amp=0.0,
                                    live_plot=false,
                                    z_init=z_init,
                                    )
        push!(compute_times, compute_time)

        # Add disturbance and simulate
        u_control = get_var(Nmpc,model,traj,:u,1)
        # Open loop control
        control_mode == :openloop ? u_control = u_rot[1] : nothing
        if dist_list == nothing
            disturbance = dist_amp*rand(model_oa.dim.u-model.dim.u) # forces applied on the body along the X and Y axes.
        else
            if l > length(dist_list)
                disturbance = zeros(model_oa.dim.u - model.dim.u)
            else
                disturbance = dist_list[l]
            end
        end
        q_ver, γ_ver, = simulate(q_truth[end-1], q_truth[end], Nsample, model_oa_rd_hf,
            u2=fill([u_control; disturbance]/Nsample, Nsample),
            ρ0=sim_ρ0, outer_iter=sim_outer_iter, z_init=z_init)

        push!(q_truth, q_ver...)
        push!(u_truth, u_control/Nsample)
        push!(dist_truth, disturbance/Nsample)
        push!(γ_truth, γ_ver...)

        # Visualization
        plot!(reshape(vcat(get_var_traj(Nmpc,model,traj,:q)...), (nq,Nmpc-1))',
            label=nothing, linewidth=3.0, linestyle=:dot)#; display(plt)
        scatter!([i/Nsample for i=1:Nsample], reshape(vcat(q_ver...), (nq,Nsample))',
            label=nothing, linewidth=3.0)#; display(plt)
        scatter!([1-1/Nsample, 1.0], reshape(vcat(q_truth[end-1:end]...), (nq,2))',
            label=nothing, markersize=5.0); display(plt)

        # Updates
        q1_init_rot = q2_init_rot
        q2_init_rot = q_truth[end]

        q1_ref_rot = q2_ref_rot
        q2_ref_rot = q_rot[1]
        rotate_seq!(u_rot)
        rotate_seq!(q_rot)
        rotate_seq!(γ_rot)
        rotate_seq!(b_rot)
        rotate_seq!(lint0)
        # typeof(model) <: WalkerNew ? q_rot[end] = mpc_stride_q(model, q_rot[end]) : nothing

        q_rot[end] = mpc_stride_q(model, SVector{nq,T}(q_rot[end]))
        if typeof(model) <: HopperAlt11 || typeof(model) <: Hopper3DAlt11 # update the height
            q_seq = view(q_truth, (length(q_truth)-Nsample+1:length(q_truth)))
            γ_seq = [γ[1] for γ in γ_truth[end-Nsample+1:end]]
            γ_max, i_max = findmax(γ_seq)
            @show scn.(γ_seq, digits=0)
            if γ_max > 0.1
                @show "***************************************in contact"
                if typeof(model) <: HopperAlt11
                    model.alt = q_seq[i_max][2] - q_seq[i_max][4] * cos(q_seq[i_max][3])
                elseif typeof(model) <: Hopper3DAlt11
                    model.alt = ϕ_no_alt(model, q_seq[i_max])[1]
                end
            end
            @show model.alt
        elseif typeof(model) <: WalkerAlt12 # update the height
            q_seq = view(q_truth, (length(q_truth)-Nsample+1:length(q_truth)))
            γ1_seq = [γ[1] for γ in γ_truth[end-Nsample+1:end]]
            γ2_seq = [γ[2] for γ in γ_truth[end-Nsample+1:end]]
            γ1_max, i1_max = findmax(γ1_seq)
            γ2_max, i2_max = findmax(γ2_seq)
            if γ1_max > 0.5
                @show "************************* 111111111111111 in contact"
                model.alt1 = ϕ_no_alt(model, q_seq[i1_max])[1]
            end
            if γ2_max > 0.5
                @show "************************* 222222222222222 in contact"
                model.alt2 = ϕ_no_alt(model, q_seq[i2_max])[2]
            end
            @show model.alt1
            @show model.alt2
        elseif typeof(model) <: QuadrupedAlt15 || typeof(model) <: WalkerAlt # update the height
            threshold = 1.0
            typeof(model) <: QuadrupedAlt15 ? threshold = 0.15 : nothing
            typeof(model) <: WalkerAlt ? threshold = 0.15 : nothing
            q_seq = view(q_truth, (length(q_truth)-Nsample+1:length(q_truth)))
            γ1_seq = [γ[1] for γ in γ_truth[end-Nsample+1:end]]
            γ2_seq = [γ[2] for γ in γ_truth[end-Nsample+1:end]]
            γ3_seq = [γ[3] for γ in γ_truth[end-Nsample+1:end]]
            γ4_seq = [γ[4] for γ in γ_truth[end-Nsample+1:end]]
            γ1_max, i1_max = findmax(γ1_seq)
            γ2_max, i2_max = findmax(γ2_seq)
            γ3_max, i3_max = findmax(γ3_seq)
            γ4_max, i4_max = findmax(γ4_seq)
            if γ1_max > threshold
                @show "************************* 111111111111111 in contact"
                model.alt1 = ϕ_no_alt(model, q_seq[i1_max])[1]
            end
            if γ2_max > threshold
                @show "************************* 222222222222222 in contact"
                model.alt2 = ϕ_no_alt(model, q_seq[i2_max])[2]
            end
            if γ3_max > threshold
                @show "************************* 333333333333333 in contact"
                model.alt3 = ϕ_no_alt(model, q_seq[i3_max])[3]
            end
            if γ4_max > threshold
                @show "************************* 444444444444444 in contact"
                model.alt4 = ϕ_no_alt(model, q_seq[i4_max])[4]
            end
            @show model.alt1
            @show model.alt2
            @show model.alt3
            @show model.alt4
        end
    end

    return Nsample, Nmpc,
        q_truth, u_truth, dist_truth, γ_truth,
        q1_init_rot, q2_init_rot, q1_ref_rot, q2_ref_rot,
        traj, u_rot, q_rot, γ_rot, b_rot, lint0, compute_times
end


function mpc_init(model,
    q1_init, q2_init, q1_ref, q2_ref, u_ref, q_ref, γ_ref, b_ref;
    Nsample::Int=3, Nmpc::Int=20, z_init=3e-2)

    N = length(u_ref)+1
    nq = model.dim.q
    nu = model.dim.u
    nγ = model.dim.γ
    nb = model.dim.b
    nr = 2nq+nu+2nγ+4nb # outer problem residual size

    q_truth = [q2_init-(q2_init-q1_init)/Nsample, q2_init]
    u_truth = []
    dist_truth = []
    γ_truth = []
    q1_init_rot = deepcopy(q1_init)
    q2_init_rot = deepcopy(q2_init)
    q1_ref_rot = deepcopy(q1_ref)
    q2_ref_rot = deepcopy(q2_ref)
    u_rot = deepcopy(u_ref)
    q_rot = deepcopy([q_ref[1:end-2]; [mpc_stride_q(model,q1_ref), mpc_stride_q(model,q2_ref)]])

    γ_rot = deepcopy(γ_ref)
    b_rot = deepcopy(b_ref)
    traj = zeros((Nmpc-1)*nr)
    set_var_traj!(Nmpc, model, traj, deepcopy(u_rot), :u)
    set_var_traj!(Nmpc, model, traj, deepcopy(q_rot), :q)
    set_var_traj!(Nmpc, model, traj, deepcopy(γ_rot), :γ)
    set_var_traj!(Nmpc, model, traj, deepcopy(b_rot), :b)

    lint0 = lintraj(q1_ref, q2_ref, deepcopy([q_ref[1:end-2]; [q1_ref, q2_ref]]),
        u_ref, model, ρ0=1e-3, outer_iter=1, z_init=z_init)

    compute_times = zeros(0)
    return Nsample, Nmpc,
        q_truth, u_truth, dist_truth, γ_truth,
        q1_init_rot, q2_init_rot, q1_ref_rot, q2_ref_rot,
        traj, u_rot, q_rot, γ_rot, b_rot, lint0, compute_times
end




# MPC options
@with_kw mutable struct MPC11Options{T}
    H_sample::Int=3
    H_mpc::Int=10
end

mutable struct MPC11{T,nq,nu,nw,nc,nb}
    q0_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q1_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q_sim::Vector{SizedArray{Tuple{nq},T,1,1,Vector{T}}} # history of simulator configurations
    u_sim::Vector{SizedArray{Tuple{nu},T,1,1,Vector{T}}} # history of simulator controls
    w_sim::Vector{SizedArray{Tuple{nw},T,1,1,Vector{T}}} # history of simulator disturbances
    γ_sim::Vector{SizedArray{Tuple{nc},T,1,1,Vector{T}}} # history of simulator impact forces
    b_sim::Vector{SizedArray{Tuple{nb},T,1,1,Vector{T}}} # history of simulator friction forces
    m_opts::MPC11Options
    ref_traj::ContactTraj
    impl::ImplicitTraj
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}} # change in q required to loop the ref trajectories
end


T = Float64
κ = 1e-4
model = get_model("quadruped")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b

# time
h = h̄
H = length(u)

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])


m_opts = MPC11Options{T}()

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

sim0 = simulator(model, q0, q1, h, H,
    u = [SVector{model.dim.u}(h * u[i]) for i=1:H],
    r! = model.res.r,
    rz! = model.res.rz,
    rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
    sim_opts = SimulatorOptions(warmstart = true))

simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
traj0 = deepcopy(sim0.traj)


function dummy_mpc(ref_traj::ContactTraj, q_stride::AbstractVector{T}, m_opts::MPC11Options)
    for l = 1:m_opts.H_mpc
        plt = plot(legend=false)
        plot!(hcat(ref_traj.q...)')
        # plot!([q[1] for q in ref_traj.q])
        rot_n_stride!(ref_traj, q_stride)
        display(plt)
    end
    return nothing
end

function rot_n_stride!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    rotate!(traj)
    mpc_stride!(traj, q_stride)
    return nothing
end

function rotate!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    @show "rot"
    q1 = copy(traj.q[1])
    u1 = copy(traj.u[1])
    w1 = copy(traj.w[1])
    γ1 = copy(traj.γ[1])
    b1 = copy(traj.b[1])
    z1 = copy(traj.z[1])
    θ1 = copy(traj.θ[1])
    for t = 1:H+1
        traj.q[t] .= traj.q[t+1]
    end
    for t = 1:H-1
        traj.u[t] .= traj.u[t+1]
        traj.w[t] .= traj.w[t+1]
        traj.γ[t] .= traj.γ[t+1]
        traj.b[t] .= traj.b[t+1]
        traj.z[t] .= traj.z[t+1]
        traj.θ[t] .= traj.θ[t+1]
    end
    traj.q[end] .= q1
    traj.u[end] .= u1
    traj.w[end] .= w1
    traj.γ[end] .= γ1
    traj.b[end] .= b1
    traj.z[end] .= z1
    traj.θ[end] .= θ1
    return nothing
end

"""
    Update the last two cofiguaraton to be equal to the fisrt two up to an constant offset.
"""
function mpc_stride!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    traj.q[end-1] .= deepcopy(traj.q[1]+q_stride)
    traj.q[end  ] .= deepcopy(traj.q[2]+q_stride)
    return nothing
end

function get_stride(model::Quadruped, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end

ref_traj0 = deepcopy(sim0.traj)
q_stride = get_stride(model, ref_traj0)

m_opts.H_mpc = 399
dummy_mpc(ref_traj0, q_stride, m_opts)
