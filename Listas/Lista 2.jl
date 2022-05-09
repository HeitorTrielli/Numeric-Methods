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

v0 = zeros(500, 7);

# A função de utilidade
utility = function(c; mu = 2)
   return (float(c)^(1 - mu) - 1)/(1 - mu)
end;

# Brute Force
value_function_brute = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital
    


    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
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

    end # while

    return value, k_line, c
end; # function

brute_force = @time value_function_brute();
brute_force_time = "21 seconds";

# Exploiting monotonicity
value_function_monotone = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0

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

value_function_concave = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital

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
concave = @time value_function_concave();


# Concavity + monotonicity
value_function = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital

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
        print(error, " ")
    end # while

    return value, k_line, c
end # function


vf = @time value_function();

# Accelerator
value_function_accelerator = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0


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

accelerator = @time value_function_accelerator()


## Accelerator sem usar monotonicidade
value_function_accelerator_brute = function(;v_0 = v0, k_len = 500, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = LinRange(0.75*kss, 1.25*kss, k_len); # Grid de capital

    # Value function, policy function on c, policy function on k and variable for iteration
    c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0


    error = 1
    count = 0

    value = copy((1 - inert)*v_next + inert*value)

    while error > tol

        Threads.@threads for state in 1:z_len    
            for capital in 1:k_len
                if count <= 30 || count%2 == 0
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    # the values of asset for which the consumption is positive
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
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
        print(error, " ")
    end

    return value, k_line, c

end # function

accel_brute = @time value_function_accelerator_brute()

plot(accel_brute[3])

# Multi grid
# Função que faz a interpolação da função valor
lin_interpol = function(minsize, maxsize, value; z_len = 7)    
    
    step = Int(maxsize/minsize)
    v_1 = zeros(maxsize, z_len)
    Threads.@threads for i in 1:(maxsize - step)
        v_1[i, :] = value[Int(ceil(i/step)),:] + (i-1)%step*(value[Int(ceil((i+step)/step)),:] - value[Int(ceil(i/step)),:])/step
    end
    Threads.@threads for i in (maxsize - step + 1):maxsize
        v_1[i, :] = value[minsize, :] - (maxsize-i)%step*(value[Int(ceil((maxsize)/step)),:] - value[Int(ceil((maxsize-step)/step)),:])/step

    end
    return v_1
end

value_function_mg = function(g1, g2, g3; v_0 = v0, z_len = 7, tol = 1e-4,
    beta = 0.987, mu = 2, alpha = 1 / 3, delta = 0.012, rho = 0.95, sigma = 0.007, inert = 0)

    value1 = value_function_brute(v_0 = zeros(g1, z_len), k_len = g1)[1]

    v1 = interpol(g1, g2, value1)

    value2 = value_function_brute(v_0 = v1, k_len = g2)[1]

    v2 = interpol(g2, g3, value2)

    value, k_line, c = value_function_brute(v_0 = v2, k_len = g3)

    return value, k_line, c
end # function

multigrid = @time value_function_mg(100, 500, 1000)

