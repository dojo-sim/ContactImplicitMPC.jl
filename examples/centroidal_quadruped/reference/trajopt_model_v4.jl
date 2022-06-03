using DirectTrajectoryOptimization 
const DTO = DirectTrajectoryOptimization

function centroidal_quadruped_dyn(model, env, h, y, x, u, w) 

    # dimensions
    nq = model.nq
    nu = model.nu 

    # configurations
    
    q1⁻ = x[1:nq] 
    q2⁻ = x[nq .+ (1:nq)]
    q2⁺ = y[1:nq]
    q3⁺ = y[nq .+ (1:nq)]

    # control 
    u_control = u[1:nu] 
    γ = u[nu .+ (1:8)] 
    β = u[nu + 8 .+ (1:32)] 
    
    E = [1.0 0.0 -1.0 0.0; 
         0.0 1.0 0.0 -1.0] # friction mapping 
    J = J_func(model, env, q2⁺)
    λ = transpose(J) * [
                        [E * β[0  .+ (1:4)]; γ[1]];
                        [E * β[4  .+ (1:4)]; γ[2]];
                        [E * β[8  .+ (1:4)]; γ[3]];
                        [E * β[12 .+ (1:4)]; γ[4]];
                        [-γ[5]; E * β[16  .+ (1:4)]];
                        [-γ[6]; E * β[20  .+ (1:4)]];
                        [-γ[7]; E * β[24  .+ (1:4)]];
                        [-γ[8]; E * β[28 .+ (1:4)]]
                       ]
    [
     q2⁺ - q2⁻;
     dynamics(model, h, q1⁻, q2⁺, u_control, zeros(model.nw), λ, q3⁺)
    ]
end

# function centroidal_quadruped_dyn1(model, env, h, y, x, u, w)
#     nx = 2 * model.nq
#     [
#      centroidal_quadruped_dyn(model, env, h, y, x, u, w);
#      y[nx .+ (1:93)] - u;
#     #  y[nx + 93 .+ (1:nx)] - x[1:nx];
#     ]
# end

function centroidal_quadruped_dynt(model, env, h, y, x, u, w)
    nx = 2 * model.nq
    [
     centroidal_quadruped_dyn(model, env, h, y, x, u, w);
     y[nx .+ (1:93)] - u;
    #  y[nx + 93 .+ (1:nx)] - x[nx + 93 .+ (1:nx)];
    ]
end

function contact_constraints_inequality_1(model, env, h, x, u, w) 
    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    u_control = u[1:nu] 
    γ = u[nu .+ (1:8)] 
    β = u[nu + 8 .+ (1:32)] 
    ψ = u[nu + 8 + 32 .+ (1:8)] 
    η = u[nu + 8 + 32 + 8 .+ (1:32)] 
    sα = u[nu + 8 + 32 + 8 + 32 .+ (1:1)]

    ϕ = ϕ_func(model, env, q3)
  
    μ = model.μ_world
    fc = μ .* γ[1:8] - [
            sum(β[0 .+ (1:4)]); 
            sum(β[4 .+ (1:4)]); 
            sum(β[8 .+ (1:4)]); 
            sum(β[12 .+ (1:4)]);
            sum(β[16 .+ (1:4)]); 
            sum(β[20 .+ (1:4)]); 
            sum(β[24 .+ (1:4)]); 
            sum(β[28 .+ (1:4)]);
            ]

    [
     -ϕ; 
     -fc;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end

function contact_constraints_inequality_t(model, env, h, x, u, w) 
    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    u_control = u[1:nu] 
    γ = u[nu .+ (1:8)] 
    β = u[nu + 8 .+ (1:32)] 
    ψ = u[nu + 8 + 32 .+ (1:8)] 
    η = u[nu + 8 + 32 + 8 .+ (1:32)] 
    sα = u[nu + 8 + 32 + 8 + 32 .+ (1:1)]

    ϕ = ϕ_func(model, env, q3)[1:8]
    γ⁻ = x[nx + nu .+ (1:8)] 
    sα⁻ = x[nx + nu + 8 + 32 + 8 + 32 .+ (1:1)]
    
    μ = model.μ_world
    fc = μ .* γ[1:8] - [
            sum(β[0 .+ (1:4)]); 
            sum(β[4 .+ (1:4)]); 
            sum(β[8 .+ (1:4)]); 
            sum(β[12 .+ (1:4)]);
            sum(β[16 .+ (1:4)]); 
            sum(β[20 .+ (1:4)]); 
            sum(β[24 .+ (1:4)]); 
            sum(β[28 .+ (1:4)]);
            ]

    [
     -ϕ; 
     -fc;
     γ⁻ .* ϕ .- sα⁻;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end


function contact_constraints_inequality_T(model, env, h, x, u, w) 
    nq = model.nq
    nu = model.nu
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    ϕ = ϕ_func(model, env, q3)[1:8]
    γ⁻ = x[nx + nu .+ (1:8)] 
    sα⁻ = x[nx + nu + 8 + 32 + 8 + 32 .+ (1:1)]
   
    [
     -ϕ; 
     γ⁻ .* ϕ .- sα⁻;
    ]
end

function contact_constraints_equality(model, env, h, x, u, w) 
    nq = model.nq
    nu = model.nu 

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    γ = u[nu .+ (1:8)] 
    β = u[nu + 8 .+ (1:32)] 
    ψ = u[nu + 8 + 32 .+ (1:8)] 
    η = u[nu + 8 + 32 + 8 .+ (1:32)] 
    sα = u[nu + 8 + 32 + 8 + 32 .+ (1:1)]
   
    E = [1.0 0.0 -1.0 0.0; 
         0.0 1.0 0.0 -1.0]
    v = (q3 - q2) ./ h[1]
    vT = [
            vcat([E' * v[6 + (i-1) * 3 .+ (1:2)] for i = 1:4]...);
            vcat([E' * v[6 + (i-1) * 3 .+ (2:3)] for i = 1:4]...);
    ]
    
    ψ_stack = vcat([ψi * ones(4) for ψi in ψ]...)
    
    [
     η - vT - ψ_stack;
    ]
end