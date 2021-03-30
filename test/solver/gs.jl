@testset "Solver: Classical Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    gs_data = ContactControl.CGSData(A, n)

    ContactControl.gs!(gs_data, A)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    A1 = deepcopy(A)
    A1.nzval .+= rand(nnz(A))
    ContactControl.gs!(gs_data, A1)
    @test norm(A1 - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = A1 \ b
    @test norm(c - x, Inf) < 1e-10
end

@testset "Solver: Modified Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    gs_data = ContactControl.MGSData(A, n)

    ContactControl.gs!(gs_data, A)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    A1 = deepcopy(A)
    A1.nzval .+= rand(nnz(A))
    ContactControl.gs!(gs_data, A1)
    @test norm(A1 - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = A1 \ b
    @test norm(c - x, Inf) < 1e-10
end


@testset "Solver: Dense Modified Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    Ad = Matrix(A)
    gs_data = ContactControl.DMGSData(Ad, n)

    ContactControl.gs!(gs_data, Ad)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10
    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = Ad \ b
    @test norm(c - x, Inf) < 1e-10

    n = 43
    m = 33
    A = rand(n,n)
    gs_data = ContactControl.DMGSData(Ad, n)
    gs!(gs_data, A)
    X = rand(n,m)
    B = rand(n,m)
    qr_matrix_solve!(gs_data, X, B)
    @test norm(A*X - B) < 1e-10
end




model = get_model("quadruped")
T = Float64
κ = 1e-4
nq = model.dim.q
nb = model.dim.b
nc = model.dim.c

ny = nb + 2nc
nx = nz - 2ny
nw = nx + ny
nz = num_var(model)
nθ = num_data(model)

ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj)
traj = deepcopy(ref_traj)
impl = ImplicitTraj(ref_traj, model, κ=κ)
implicit_dynamics!(impl, model, traj; κ=κ)

z_initialize!(impl.ip[1].z, model, traj.q[1+2]) #TODO: try alt. schemes
impl.ip[1].θ .= traj.θ[1]
status = interior_point!(impl.ip[1])

ip0 = deepcopy(impl.ip[1])
z_initialize!(ip0.z, model, traj.q[1+2])
ip0.θ .= traj.θ[1]

interior_point!(ip0)
z0  = deepcopy(traj.z[1])
z0  = rand(nz)
θ0  = deepcopy(traj.θ[1])
δz0 = zeros(nz)
r0  = zeros(nz)
rz0 = zeros(nz,nz)
rθ0 = zeros(nz,nθ)
r_cache0  = deepcopy(ip0.r_cache)
rz_cache0 = deepcopy(ip0.rz_cache)
rθ_cache0 = deepcopy(ip0.rθ_cache)

ip0.methods.r!(r0, z0, θ0, κ, r_cache0)
ip0.methods.rz!(rz0, z0, θ0, rz_cache0)
ip0.methods.rθ!(rθ0, z0, θ0, rθ_cache0)

r0
gs_data = DMGSData(rz0, nz)
gs!(gs_data, rz0)
qr_solve!(gs_data, δz0, r0)
@test norm(rz0 \ r0 - δz0, Inf) < 1e-10


off = 0
ibil = [ibil1; ibil3; ibil2]
ibil1 = Vector(nq + nc .+ (1:nc))
ibil2 = Vector(nq + 3nc + nb .+ (1:nc))
ibil3 = Vector(nq + 4nc + nb .+ (1:nb))
ilin = setdiff(1:nz, ibil)
idyn = Vector(1:nq)
irst = setdiff(ilin, idyn)
ix = off .+ (1:nq); off += nq
iy1 = off .+ (1:ny); off += ny
iy2 = off .+ [Vector(nb .+ (1:nc)); Vector(1:nb); Vector(nb+nc .+ (1:nc))]; off += ny


plot(Gray.((1e10*abs.(rz0[:,:]))))
plot(Gray.((1e10*abs.(rz0[ilin,:]))))
plot(Gray.((1e10*abs.(rz0[ibil,:]))))
plot(Gray.((1e10*abs.(rz0[ibil,ix]))))
plot(Gray.((1e10*abs.(rz0[ibil,iy1]))))
plot(Gray.((1e10*abs.(rz0[ibil,iy2]))))

Ax  = rz0[ilin,ix]
Ay1 = rz0[ilin,iy1]
Ay2 = rz0[ilin,iy2]
Dx  = rz0[idyn, ix]
Dy1 = rz0[idyn, iy1]
Rx  = rz0[irst, ix]
Ry1 = rz0[irst, iy1]
Ry2 = rz0[irst, iy2]

Y1 = Diagonal(diag(rz0[ibil,iy2]))
Y2 = Diagonal(diag(rz0[ibil,iy1]))
rdyn = r0[idyn]
rrst = r0[irst]
r1 = r0[ilin]
r2 = r0[ibil]

δz1 = zeros(nz)
δz1[[ix;iy1]] .= [Ax Ay1-Ay2*(Y1\Y2)]\(r1 - Ay2*(Y1\r2))
δz1[iy2] = Y1 \ (r2 - Y2*δz1[iy1])


norm(δz0 - δz1)
norm(δz0[ix] - δz1[ix])
norm(δz0[iy1] - δz1[iy1])
norm(δz0[iy2] - δz1[iy2])

Y2*δz0[iy1] + Y1*δz0[iy2] - r0[ibil]
Y2*δz1[iy1] + Y1*δz1[iy2] - r0[ibil]

plot(Gray.((1e10*abs.(rz0[idyn,[ix; iy1; iy2]]))))
plot(Gray.((1e10*abs.(rz0[idyn,iy1]))))
plot(Gray.((1e10*abs.(rz0[idyn,iy2]))))
plot(Gray.((1e10*abs.(rz0[idyn,iy2]))))
plot(Gray.((1e10*abs.(rz0[irst,[ix;iy1;iy2]]))))
plot(Gray.((1e10*abs.(rz0[irst,ix]))))
plot(Gray.((1e10*abs.(rz0[irst,iy1]))))
plot(Gray.((1e10*abs.(rz0[irst,iy2]))))
plot(Gray.((1e10*abs.(rz0[[ilin; ibil],[ix; iy1; iy2]]))))

@test norm(rz0 \ r0 - δz1, Inf) < 1e-10
@test norm(δz0 - δz1, Inf) < 1e-10

rz0[irst,iy2]


function schur_solve(A,B,C,D,Ai,u,v)
    S = D - C*Ai*B
    Si = inv(S)
    x = (Ai + Ai*B*Si*C*Ai) * u - Ai*B*Si * v
    y = - Si*C*Ai * u + Si * v
    @assert norm([A B; C D]*[x;y] - [u;v], Inf) < 1e-10
    return x, y
end

δx_2, δy1_2 = schur_solve(Dx, Dy1, Rx, Ry1 - Ry2*Y1\Y2, inv(Dx), rdyn, rrst - Ry2*Y1\r2)
δy2_2 = Y1 \ (r2 - Y2*δy1_2)
δz2 = zeros(nz)
δz2[ix]  .= δx_2
δz2[iy1] .= δy1_2
δz2[iy2] .= δy2_2
@test norm(δz0 - δz2, Inf) < 1e-10

mutable struct StructuredSolver13{T}
    rz::Matrix{T}
    Dx::SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Array{Int64,1}},false}
    Dy1::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    Rx::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    Ry1::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    Ry2::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    Y1::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    Y2::SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false}
    # Residual
    r::Vector{T}
    r1::SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false}
    r2::SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false}
    rdyn::SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false}
    rrst::SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false}
end

function StructuredSolver13(model::ContactDynamicsModel, rz::AbstractMatrix{T}) where {T}
    nq = model.dim.q
    nc = model.dim.c
    nb = model.dim.b
    nz = num_var(model)

    # Terms
    off = 0
    ibil1 = Vector(nq + nc .+ (1:nc))
    ibil2 = Vector(nq + 3nc + nb .+ (1:nc))
    ibil3 = Vector(nq + 4nc + nb .+ (1:nb))
    ibil = [ibil1; ibil3; ibil2]
    ilin = setdiff(1:nz, ibil)
    idyn = Vector(1:nq)
    irst = setdiff(ilin, idyn)
    # Vars
    ix = off .+ Vector(1:nq); off += nq
    iy1 = off .+ Vector(1:ny); off += ny
    iy2 = off .+ [Vector(nb .+ (1:nc)); Vector(1:nb); Vector(nb+nc .+ (1:nc))]; off += ny

    # Matrix views
    Dx  = view(rz, idyn, ix)
    Dy1 = view(rz, idyn, iy1)
    Rx  = view(rz, irst, ix)
    Ry1 = view(rz, irst, iy1)
    Ry2 = view(rz, irst, iy2)
    Y1  = view(rz, ibil, iy2)
    Y2  = view(rz, ibil, iy1)
    Y2  = view(rz, ibil, iy1)

    # Vector views
    r = zeros(nz)
    r1 = view(r, ilin)
    r2 = view(r, ibil)
    rdyn = view(r, idyn)
    rrst = view(r, irst)
    @show typeof(r)
    @show typeof(r1)
    @show typeof(r2)
    @show typeof(rdyn)
    @show typeof(rrst)

    @show typeof(Dx)
    @show typeof(Dy1)
    @show typeof(Rx)
    @show typeof(Ry1)
    @show typeof(Ry2)
    @show typeof(Y1)
    @show typeof(Y2)
    @show typeof(Y2)
    return StructuredSolver13{T}(rz,Dx,Dy1,Rx,Ry1,Ry2,Y1,Y2,r,r1,rdyn,rrst,r2)
end


solver = StructuredSolver13(model, rz0)



typeof(view(r0, idyn))
typeof(view(rz0, idyn, ix))






a = 10
a = 10
a = 10
a = 10
a = 10

# function nlog(x; h=10.0)
#     if x < 0
#         return -log.(h, -x)
#     elseif x == 0
#         return 0.0
#     elseif x > 0
#         return log.(h, x)
#     end
#     return nothing
# end
#
# function plot_mat(A)
#     Ap = A
#     mi = minimum(Ap)
#     Ap = Ap .- mi .+ 1e-5
#     Ap = log.(10, Ap)
#     ma = maximum(Ap)
#     Ap ./= ma
#     plt = plot()
#     plot!(Gray.((Ap)))
#     @show Ap
#     display(plt)
#     return nothing
# end
#
# plot_mat(rz0[1:10,1:10])
#
# rz0[1:10, 1:10] - rz0[1:10, 1:10]'
#
# rz0[1:10,1:10]
# heatmap(Matrix{T}(nlog.(rz0[1:10,1:10])), color = :greys)
# heatmap(Matrix{T}(nlog.(rz0)), color = :greys)





model
n = 43-16
n = 16
Random.seed!(100)
A = sprand(n, n, 0.16)
while rank(A) < n
    A = sprand(n, n, 0.16)
end
Ad = Matrix(A)
gs_data = ContactControl.DMGSData(Ad, n)

@benchmark ContactControl.gs!(gs_data, Ad)
@test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10
b = rand(SizedVector{n,T})
c = zeros(SizedVector{n,T})
ContactControl.qr_solve!(gs_data, c, b)
x = Ad \ b
@test norm(c - x, Inf) < 1e-10

n = 43
m = 33
A = rand(n,n)
gs_data = ContactControl.DMGSData(Ad, n)
gs!(gs_data, A)
X = rand(n,m)
B = rand(n,m)
qr_matrix_solve!(gs_data, X, B)
@test norm(A*X - B) < 1e-10
