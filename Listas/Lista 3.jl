using Plots, Distributions # Pacotes que estou usando

#start
    tauchen = function(grid_len::Int64; mu::Float64 = 0.0, sigma::Float64 = 0.007, rho::Float64 = 0.95, m::Float64 = 3.0)

        theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
        theta_min = - theta_max + 2 * mu # definindo o menor valor do grid
        grid = Array(LinRange(theta_min, theta_max, grid_len)) # Cria um vetor de n pontos entre theta_max e theta_min em que a distancia entre os pontos sequencias é igual

        d = Normal(mu, 1) # d vira normal(mu,1), que será usado para computar a PDF dos erros na hora de achar as probabilidades de transição
        delta = (maximum(grid) - minimum(grid)) / (length(grid) - 1) # distância dos pontos subsequentes do grid

        vec_1 = cdf(d,((minimum(grid) .- rho * grid .+ delta / 2) / sigma)) # vetor das probabilidades de ir para o menor valor do grid, dado cada estado anterior do grid; cdf(d, x) retorna a cdf da distribuição d no valor x
        vec_n = 1 .- cdf(d,((maximum(grid) .- rho * grid .- delta / 2) / sigma)) # análogo para o maior valor do grid
        grid_interno = grid[2:(length(grid) - 1)] # valores não extremos do grid

        pij = function(j, i = grid) # função que vai computar o vetor de probabilidades de ir para o estado (não extremo) j dado cada estado anterior do grid
            cdf(d,((j + delta/2 .- rho * i) / sigma)) - cdf(d,((j - delta / 2 .- rho * i) / sigma))                             
        end

        mat_interna = reduce(hcat, pij.(grid_interno))  # aplica pij em cada ponto do grid interno; reduce: transforma o vetor de vetores em uma matriz
        
        probs = [vec_1 mat_interna vec_n] # gerando a matriz de transição

        return probs::Array{Float64,2}, Array(grid)::Array{Float64,1}
    end;
    # Definindo as variáveis
    beta, mu, alpha, delta, rho, sigma = 0.987, 2.0, 1/3, 0.012, 0.95, 0.007;

    # Definindo o capital de steady state

    k_ss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));
    z_len = 7;
    k_len = 500;
    grid_z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    grid_k = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len));
    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len)); # Grid de capital

    utility = function(c::Float64; mu::Float64 = 2.0)
        return (c^(1 - mu) - 1)/(1 - mu)
    end;

    u_line = function(c, mu = 2.0)
        c^(-mu)
    end
    
    u_inv = function(c, mu = 2.0)
        c^(-1/mu)
    end
#

# Normalizando o K para poder jogar nos polinomios de Chebyschev
knorm = (k .- k_ss)/((k.-k_ss)[k_len])


# Função que retorna o polinômio de Chebyschev de grau j avaliado no ponto x
cheb = function(x, j)
    fun = cos(j*acos(x))
    return fun
end


# Função que acha as raízes do polinomio de Chebyschev de grau j
chebroots = function(j)
    roots = zeros(j)
    for i in 1:j
        roots[i] = -cos((2*i - 1)/(2*j)*pi) 
    end
    return roots
end

gamma = ones(7,5)/5
gamma[1,:] = [0.1, 0.2, 0.3, 0.2, 0.2]

# Função que monta o consumo em função de gamma, k e d polinômios de chebyschev
cons = function(gamma, k, z, d)
    cons = zeros(d)
    for i in 1:d
        cons[i] = cheb(k, (i-1))
    end  
    return gamma[z,:]'cons
end

cons.((gamma,), 0.5, 1:7, 5)

grid_z*k^alpha - cons.((gamma,), k, 1:7, d) .+ (1 - delta)*k


# Função que retorna os resíduos dado a matriz de gammas, o capital, o estado e o numero de polinomios de chebyschev que estamos usando
resid = function(gamma, k, z, d; p_z = p_z, beta = beta, z_grid = grid_z, alpha = alpha, delta = delta)
    # (gamma,) serve para tratar a matriz como uma tupla e poder aplicar o broeadcast sem dar erro (não me pergunte, eu também não entendo)
    k1 = z_grid*k^alpha - cons.((gamma,), k, 1:7, d) .+ (1 - delta)*k
    cline = zeros(7)
    for i in 1:7
        cline[i] = cons((gamma,), k1[i], i, d)
    end

    u_line(cons(gamma, k, z, d)) - beta*((u_line.(cline) .* (z_grid * (k ^ (alpha - 1) * alpha .+ (1 - delta))))'p_z[z, :])
    return cline
end

resid(gamma, 0.5, 3, 5)


end


f(x,y) = sum(x) + y

f([1,2,3], 4) #-> 10
f([1,2,3], 5) #-> 11

f.([1,2,3], [4,5])    # -> error

f.(([1,2,3],), [4,5]) # -> [10, 11]
                      # this works because ([1,2,3],) is now considered
                      # to be scalar and is broadcasted