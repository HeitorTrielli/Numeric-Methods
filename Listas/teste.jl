include("Lista 1.jl"); # Para usar as funções definidas na lista 1

# Definindo as variáveis
beta, mu, alpha, delta, rho, sigma = 0.987, 2, 1/3, 0.012, 0.95, 0.007;

# Definindo o capital de steady state
kss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));

# Definindo os grids
grid_z = exp.(tauchen(7)[2]); # Valores que exp(z_t) pode tomar 
prob_z = tauchen(7)[1]; # Matriz de transição de z
grid_k = LinRange(0.75*kss, 1.25*kss, 500); # Grid de capital

# A função de utilidade
utility = function(c; mu = 2)
   return (float(c)^(1 - mu) - 1)/(1 - mu)
end

error = 1
z = repeat(grid_z, 1, k_len)'
p_z = prob_z
k = grid_k
tol = 10e-4    
k_len = length(grid_k)
z_len = length(grid_z)
kmat = repeat(grid_k, 1, z_len)

teste = [1 2 3; 4 5 6; 7 8 9]
teste.*(teste .> 1)

# Value function, policy function on c, policy function on k and variable for iteration

value, c, k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len), 
                            zeros(k_len, z_len), zeros(k_len, z_len) 
count = 0
state = 1
capital = 1            
k_possible = kmat.*(z*(k[capital]^alpha) .- k .+ (1 - delta)*k[capital] .> 0)    # the values of asset for which the consumption is positive


v_possible = value[z[state]*(k[capital]^alpha) .- k .+ (1 - delta)*k[capital] .> 0, :]    #the value function at each of the assets above              
                            
if  count%10 == 0 || count < 20 
    val = utility.(z[state]*(k[capital]^alpha) .- k_possible .+ (1 - delta)*k[capital]) .+ beta*v_possible*p_z[state, :]
    v_next[capital, state] = maximum(val)
    k_line[capital, state] = k_possible[argmax(val)]
    c[capital, state] = z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]
else
    v_next[capital, state] = utility(c[capital, state]) + beta*v_possible[state,:]'*p_z[state, :]
end

utility(c[capital, state])

v_possible[state, :]
