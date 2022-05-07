using Plots, BenchmarkTools, Distributions, Profile # Pacotes que estou usando

Threads.nthreads()

tauchen = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95, m = 3)

    theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
    theta_min = - theta_max + 2 * mu # definindo o menor valor do grid
    grid = LinRange(theta_min, theta_max, grid_len) # Cria um vetor de n pontos entre theta_max e theta_min em que a distancia entre os pontos sequencias é igual

    d = Normal(mu, 1) # d vira normal(mu,1), que será usado para computar a PDF dos erros na hora de achar as probabilidades de transição
    delta = (maximum(grid) - minimum(grid)) / (length(grid) - 1) # distância dos pontos subsequentes do grid

    vec_1 = cdf(d,((minimum(grid) .- rho * grid .+ delta / 2) / sigma)) # vetor das probabilidades de ir para o menor valor do grid, dado cada estado anterior do grid; cdf(d, x) retorna a cdf da distribuição d no valor x
    vec_n = 1 .- cdf(d,((maximum(grid) .- rho * grid .- delta / 2) / sigma)) # análogo para o maior valor do grid
    grid_interno = grid[2:(length(grid) - 1)] # valores não extremos do grid

    pij = function(j, i = grid) # função que vai computar o vetor de probabilidades de ir para o estado (não extremo) j dado cada estado anterior do grid
        cdf(d,((j .+ delta/2 .- rho * i) / sigma)) - cdf(d,((j .- delta / 2 .- rho * i) / sigma))                             
    end

    mat_interna = reduce(hcat, pij.(grid_interno))  # aplica pij em cada ponto do grid interno; reduce: transforma o vetor de vetores em uma matriz
    
    probs = [vec_1 mat_interna vec_n] # gerando a matriz de transição

    return probs, grid
        
end;
# Definindo as variáveis
beta, mu, alpha, delta, rho, sigma = 0.987, 2, 1/3, 0.012, 0.95, 0.007;

# Definindo o capital de steady state
kss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));

# Definindo os grids
grid_z = exp.(tauchen(7)[2]); # Valores que exp(z_t) pode tomar 
prob_z = tauchen(7)[1]; # Matriz de transição de z
grid_k = LinRange(0.75*kss, 1.25*kss, 500); # Grid de capital
v0 = zeros(length(grid_k), length(grid_z));

# A função de utilidade
utility = function(c; mu = 2)
   return (float(c)^(1 - mu) - 1)/(1 - mu)
end;

# Brute Force
value_function_brute = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
    count = 0
    error = 1
    while error > tol

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
        value = copy((1 - inert)*v_next + inert*value)

        count = count + 1
    end # while

    return value, k_line, c
end; # function

brute_force = @time value_function_brute()
brute_force_time = "21 seconds";

# Exploiting monotonicity
value_function_monotone = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 1e-4, 
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)
    value = v_0

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)

    error = 1

    for state in 1:z_len            
        for capital in 1:k_len
            k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k.> 0]    # the values of asset for which the consumption is positive
            v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]    #the value function at each of the assets above              
            val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
            v_next[capital, state] = maximum(val)
            k_line[capital, state] = k_possible[argmax(val)]
            c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
        end # for k
    end # for z

    value = copy((1 - inert)*v_next + inert*value)

    while error > tol

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state]]    # the values of asset for which the consumption is positive
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state], :]               
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state] = maximum(val)
                k_line[capital, state] = k_possible[argmax(val)]
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)

    end # while

    return value, k_line, c
end # function

monotone = @time value_function_monotone();
monotone_time = "6.3 seconds";


# Concavity
lastpos = function(val, k_possible) # Função que vai retornar o ponto maximo da função e o seu respectivo k'
    v = 0
    k = 0
    for i in 2:length(val)
        if val[i] < val[i-1]
            v = val[i-1]
            k = k_possible[i-1]
            break
        else
            v = val[i]
            k = k_possible[i]
        end
    end
    return v, k
end;

value_function_concave = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    k_len = length(grid_k)
    z_len = length(grid_z)

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
    count = 0
    error = 1
    while error > tol

        Threads.@threads for state in 1:z_len
            for capital in 1:k_len  
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state], k_line[capital, state] = lastpos(val, k_possible)
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
        count = count + 1
    end # while

    return value, k_line, c
end; # function
concave = @time value_function_concave()


# Concavity + monotonicity
value_function = function(;v_0 = v0, z = grid_z, p_z = prob_z, k = grid_k, tol = 1e-4, 
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)
    value = v_0

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)

    error = 1

    for state in 1:z_len            
        for capital in 1:k_len
            k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k.> 0]    # the values of asset for which the consumption is positive
            v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]    #the value function at each of the assets above              
            val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
            v_next[capital, state], k_line[capital, state] = lastpos(val, k_possible)
            c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
        end # for k
    end # for z

    value = copy((1 - inert)*v_next + inert*value)

    while error > tol

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state]]    # the values of asset for which the consumption is positive
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state], :]               
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state], k_line[capital, state] = lastpos(val, k_possible)
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)

    end # while

    return value, k_line, c
end # function


vf = @time value_function()
monotone = @time value_function_monotone()

# Accelerator
value_function_accelerator = function(; z = grid_z, p_z = prob_z, k = grid_k, tol = 1e-4,
    beta = 0.987, mu =2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)
    
    k_len = length(grid_k)
    z_len = length(grid_z)

    value, c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len) 

    error = 1
    count = 0
    
    Threads.@threads for state in 1:z_len    
        for capital in 1:k_len
            k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k.> 0]    # the values of asset for which the consumption is positive
            v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]    #the value function at each of the assets above              
            val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
            v_next[capital, state] = maximum(val)
            k_line[capital, state] = k_possible[argmax(val)]
            c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
        end # for k
    end # for z

    value = copy((1 - inert)*v_next + inert*value)

    while error > tol

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len
                if count <= 10 || count%10 == 0
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state]]    # the values of asset for which the consumption is positive
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k_line[capital, state], :]               
                    val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                    v_next[capital, state] = maximum(val)
                    k_line[capital, state] = k_possible[argmax(val)]
                    c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
                else
                    v_next[capital, state] = utility(c[capital,state]) + beta*value[findall(k_line[capital,state] .== k)[1],:]'*p_z[state, :]
                end
            end # for k
        end # for z
 
    
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
        count = count + 1
    end

    return value, k_line, c

end # function

accelerator = @time value_function_accelerator();

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


