using BenchmarkTools, Plots
include("Lista 1.jl"); # Para usar as funções definidas na lista 1

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
end


teste = [1 2 3 4 5]
@benchmark @. utility(grid_k - 1)
@benchmark utility.(grid_k .- 1)

z = grid_z
p_z = prob_z
k = grid_k
tol = 10e-4    
k_len = length(grid_k)
z_len = length(grid_z)
value, c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len) 

count = 0


# Concavity
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

brute = @time value_function_brute()

state = 5
capital = 500

k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
lagval = replace(lag(val), missing => 0.0)
index = findlast(val - lagval .> 0)
v_next[capital, state] = val[index]
k_line[capital, state] = k_possible[index]
c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]

@btime Threads.@threads for state in 1:z_len
            for capital in 1:k_len  
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
                val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
                v_next[capital, state] = maximum(val)
                k_line[capital, state] = k_possible[argmax(val)]
                c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
            end # for k
        end # for z

@btime for i in 1:1
    Threads.@threads for state in 1:z_len
        for capital in 1:k_len  
            k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0]    
            v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]               
            val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
            v_next[capital, state] = maximum(val)
            k_line[capital, state] = k_possible[argmax(val)]
        end # for k
    end # for z
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat

end



lastpos = function(val, k_possible)
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
end

teste, test = lastpos(val, k_possible)

@btime lastpos(val, k_possible)
@btime maximum(val)


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

        count += 1
    end # while

    return value, k_line, c
end; # function

k_line

zmat = repeat(z,1,k_len)'
kmat = repeat(k,1,z_len)

zmat.*kmat - k_line + (1-delta)*kmat

zmat
kmat


count = 0

