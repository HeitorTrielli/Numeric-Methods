using Plots, BenchmarkTools
include("Lista 1.jl"); # Para usar as funções definidas na lista 1
Threads.nthreads()
# Definindo as variáveis
beta, mu, alpha, delta, rho, sigma = 0.987, 2, 1/3, 0.012, 0.95, 0.007;

# Definindo o capital de steady state
kss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));

# Definindo os grids
grid_z = exp.(tauchen(7)[2]); # Valores que exp(z_t) pode tomar 
prob_z = tauchen(7)[1]; # Matriz de transição de z
grid_k = LinRange(0.75*kss, 1.25*kss, 500); # Grid de capital
v0 = zeros(length(grid_k), length(grid_z))

# A função de utilidade
utility = function(c; mu = 2)
   return (float(c)^(1 - mu) - 1)/(1 - mu)
end


# Brute Force
value_function_brute = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 10e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0

    error = 1
    count = 0
    while error > tol
        value = copy((1 - inert)*v_next + inert*value)

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len  
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state] = maximum(val)
                k_line[capital, state] = k_possible[argmax(val)]
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        count = count + 1
    end # while

    return value, k_line, c
end # function

brute_force = @btime value_function_brute();
brute_force_time = "25 seconds";

# Exploiting monotonicity
value_function_monotone = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 10e-4, 
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)
    value = v_0

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)

    error = 1
    count = 0
    while error > tol
        value = copy((1 - inert)*v_next + inert*value)

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len
                
                if count == 0
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k.> 0]    # the values of asset for which the consumption is positive
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]    #the value function at each of the assets above              
                else 
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state]]    # the values of asset for which the consumption is positive
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state], :]
                end
                
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state] = maximum(val)
                k_line[capital, state] = k_possible[argmax(val)]
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        count = count + 1
    end # while

    return value, k_line, c
end # function

monotone = @btime value_function_monotone();
monotone_time = "7.5 seconds";


# Accelerator
value_function_accelerator = function(; z = grid_z, p_z = prob_z, k = grid_k, tol = 10e-2,
    beta = 0.987, mu =2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)

    value, c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len) 

    error = 1
    count = 0

    while error > tol
        value = copy((1 - inert)*v_next + inert*value)

            for state in 1:z_len    
                for capital in 1:k_len
                    if  count%10 == 0 || count < 100
                        k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
                        v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]    
                        val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital]) .- k_possible .+ beta*v_possible*p_z[state, :]
                        v_next[capital, state] = maximum(val)
                        k_line[capital, state] = k_possible[argmax(val)]
                        c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
                    else                     
                        v_next[capital, state] = utility(c[capital, state]) + beta*value[capital,:]'*p_z[state, :]
                    end

                end # for k
            end # for z
 
        count = count + 1
        error = maximum(abs.(value - v_next))
        print(error, " ")


    end # while

    return value, k_line, c

end # function

accelerator = @time value_function_accelerator()


plot(accelerator[1])

# Multi grid
value_function_mg = function(g1, g2, g3; z = grid_z, p_z = prob_z, tol = 10e-4,
    beta = 0.987, mu =2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007)
    
    grid_k1 = LinRange(0.75*kss, 1.25*kss, g1)
    grid_k2 = LinRange(0.75*kss, 1.25*kss, g2)
    grid_k3 = LinRange(0.75*kss, 1.25*kss, g3)
    k1_len = length(grid_k1)
    k2_len = length(grid_k2)
    k3_len = length(grid_k3)
    z_len = length(grid_z)

    value, c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len) 

    error = 1
    count = 0
    while error > tol
        value = copy(v_next)

        for state in 1:z_len    
            for capital in 1:k_len
                if count % 10 == 0 || count < 20
                    k_possible = k[z[state]*(k[capital]^alpha) .- k .+ (1 - delta)*k[capital] .> 0]    # the values of asset for which the consumption is positive
                    v_possible = value[z[state]*(k[capital]^alpha) .- k .+ (1 - delta)*k[capital] .> 0, :]   #the value function at each of the assets above              
                    val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                    v_next[capital, state] = maximum(val)
                    k_line[capital, state] = k_possible[argmax(val)]
                    c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
                else
                    v_next = utility.(c) + beta*p_z*value
                end # if
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        count = count + 1
    end # while

    return value, k_line, c
end # function


