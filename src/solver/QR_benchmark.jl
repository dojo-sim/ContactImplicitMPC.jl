using StaticArrays
using Random

Random.seed!(10)

# data
T = Float64
n = 25
A = rand(SMatrix{n,n,T,n^2})
M = A * A'
b = rand(SVector{n,T})

# Basic inverse
x0 = inv(M) * b
@benchmark $x0 = inv($M) * $b

# Basic \
x1 = M \ b
@benchmark $x1 = $M \ $b
@show norm(x1 - x0)

# StaticArrays QR
qrM = StaticArrays.qr(M)
@benchmark $qrM = StaticArrays.qr($M)
x2 = qrM.R \ (qrM.Q' * b)
@benchmark $x2 = $qrM.R \ ($qrM.Q' * $b)
@show norm(x2 - x0)

# Custom QR
dat = SDMGSData(M)
factorize!(dat, M)
@benchmark factorize!(dat, M)
qr_solve!(dat, b)
@benchmark qr_solve!(dat, b)
x3 = dat.xs
@show norm(x3 - x0)
