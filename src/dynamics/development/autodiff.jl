using AbstractTrees

# Define functions
function f(q0)
    t0 = sum(q0)^2
    t1 = exp.(q0)
    t2 = sin.(q0)
    return [sum(log(t0) * sum(t0 .* t1 .* t2))]
end

function g(q1, Q0)
    nq = size(Q0)[1]
    q2 = reshape(Q0 .* q1, (nq^2))
    return q2
end

function h(q0, Q0)
    g(f(q0), Q0)
end

# Sizes
nq = 5

# Variables
@variables q0[1:nq]
@variables q1m[1:1]
@variables Q0[1:nq,1:nq]

################################################################################
# Evaluation
################################################################################

# Generate expressions
q1 = f(q0);
q2m = g(q1m, Q0);
q2 = h(q0, Q0);

expr = Dict()
expr[:f] = build_function(q1, q0)[2]
expr[:g] = build_function(q2m, q1m, Q0)[2]
expr[:h] = build_function(q2, q0, Q0)[2]

fs = eval(expr[:f])
gs = eval(expr[:g])
hs = eval(expr[:h])

# Test functions
Random.seed!(100)
q0_ = rand(nq)
Q0_ = rand(nq, nq)
q1m_ = zeros(1)
q2_ = zeros(nq, nq)
q2m_ = zeros(nq, nq)

fs(q1m_, q0_)
gs(q2m_, q1m_, Q0_)
hs(q2_, q0_, Q0_)

@test norm(q2_ - q2m_) < 1e-10

tf = @belapsed fs(q1m_, q0_)
tg = @belapsed gs(q2m_, q1m_, Q0_)
th = @belapsed hs(q2_, q0_, Q0_)
thm = tf + tg
th / thm


################################################################################
# Jacobian
################################################################################

# Generate expressions
q1 = f(q0);
q2m = g(q1m, Q0);
q2 = h(q0, Q0);

∇q0_q1 = Symbolics.jacobian(q1, q0, simplify=true)
∇q1_q2 = Symbolics.jacobian(q2m, q1m, simplify=true)
∇q0_q2 = Symbolics.jacobian(q2, q0, simplify=true)
∇q0_q2m = ∇q1_q2 * ∇q0_q1

expr = Dict()
expr[:f] = build_function(q1, q0)[2]
expr[:g] = build_function(q2m, q1m, Q0)[2]
expr[:h] = build_function(q2, q0, Q0)[2]
expr[:∇f] = build_function(∇q0_q1, q0)[2]
expr[:∇g] = build_function(∇q1_q2, q1m, Q0)[2]
expr[:∇h] = build_function(∇q0_q2, q0, Q0)[2]

fs = eval(expr[:f])
gs = eval(expr[:g])
hs = eval(expr[:h])
∇fs = eval(expr[:∇f])
∇gs = eval(expr[:∇g])
∇hs = eval(expr[:∇h])

# Test functions
Random.seed!(100)
q0_ = rand(nq)
Q0_ = rand(nq, nq)
q1m_ = zeros(1)
q2_ = zeros(nq^2)
q2m_ = zeros(nq^2)
∇q0_q1_ = zeros(1,nq)
∇q1_q2_ = zeros(nq^2,1)
∇q0_q2_ = zeros(nq^2,nq)
∇q0_q2m_ = zeros(nq^2,nq)

fs(q1m_, q0_)
gs(q2m_, q1m_, Q0_)
hs(q2_, q0_, Q0_)
∇fs(∇q0_q1_, q0_)
∇gs(∇q1_q2_, q1m_, Q0_)
∇hs(∇q0_q2_, q0_, Q0_)
∇q0_q2m_ = ∇q1_q2_ * ∇q0_q1_

@test norm(q2_ - q2m_) < 1e-10
@test norm(∇q0_q2_ - ∇q0_q2m_) < 1e-10

tf = @belapsed fs(q1m_, q0_)
tg = @belapsed gs(q2m_, q1m_, Q0_)
th = @belapsed hs(q2_, q0_, Q0_)
thm = tf + tg
th / thm

t∇f = @belapsed ∇fs(∇q0_q1_, q0_)
t∇g = @belapsed ∇gs(∇q1_q2_, q1m_, Q0_)
t∇h = @belapsed ∇hs(∇q0_q2_, q0_, Q0_)
tmult = @belapsed ∇q0_q2m_ = ∇q1_q2_ * ∇q0_q1_
t∇hm = t∇f + t∇g + tmult
t∇h / t∇hm
