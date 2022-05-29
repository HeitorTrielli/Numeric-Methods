using Plots, Distributions, Optim, NLsolve, BenchmarkTools, ProfileView # Pacotes que estou usando

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

#####################################
############# Questão 1 #############
#####################################

# Normalizando o K para ficar entre [-1, 1]
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

# Função que monta o consumo em função de gamma, k e d polinômios de chebyschev
cons = function(gamma, k, z, d)
    cons = zeros(d)
    for i in 1:d
        cons[i] = cheb(k, (i-1))
    end  
    return gamma[z,:]'cons
end

# Tentativa de fazer a função de residuos sem usar for
    # Função que retorna os resíduos dado a matriz de gammas, o capital, o estado e o numero de polinomios de chebyschev que estamos usando
    #resid = function(gamma; knorm = knorm, k = grid_k, p_z = p_z, beta = beta, z_grid = grid_z, alpha = alpha, delta = delta, k_len = k_len)
    #    d = Int(length(gamma)/7)
    #    root = chebroots(d)
    #    # Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de chebyschev de grau d
    #    # Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
    #    k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
    #    k0level = k0*(k .- k_ss)[500].+k_ss
    #    
    #    # Aplica a função cons (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
    #    # (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
    #    c0 = reduce(hcat, [cons.((gamma,), k0, 1:7, d) for k0 in k0])
    #    k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
    #    k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[500]) for j in 1:d for i in 1:7 ], 7, d)
    #
    #    k1norm = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:7], 7, d))[i, j]) for j in 1:d for i in 1:7], 7, d))
    #    
    #    c1 = [cons.((gamma,), k1norm[i, j], 1:7, d) for j in 1:d for i in 1:7 ]
    #    c1 = reshape(c1, 7, d)
    #    resid = reshape([u_line(c0[i]) - beta*(u_line.(c1[i]).*z_grid*(k1level[i]^(alpha-1)*alpha + 1 - delta))'p_z[i % 7 != 0 ? i % 7 : 7, :] for i in 1:(d*7)], 7, d)
    #    return resid
    # end


# Função que retorna os resíduos
R = function(gamma; knorm = knorm, k = grid_k, p_z = p_z, beta = beta, z_grid = grid_z, alpha = alpha, delta = delta, k_len = k_len)
    d = Int(length(gamma)/7) # Número de polinomios que vou usar
    root = chebroots(d) # Raízes do polinomio de chebyschev de grau d
    k0 = root
    k0level = k0*(k .- k_ss)[500].+k_ss
    c0 = zeros(7,d)
    k1level = zeros(7,d)
    k1norm = zeros(7,d)
    resid = zeros(7,d)
    c1 = zeros(7)
    for state in 1:7
        for w in 1:d
            c0[state, w] = cons(gamma, k0[w], state, d)
            k1level[state, w] = z[state]*(k0level[w]^alpha) + (1 - delta)*k0level[w] - c0[state,w]
            k1norm[state, w] = (k1level[state, w] - k_ss)/((k .- k_ss)[500])
            c1 = cons.((gamma,), k1norm[state,w], 1:7, d) 
            resid[state, w] = u_line(c0[state, w]) - beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])
        end
    end
    return resid
end


cheby = function(d)
    gamma = ones(7,2)
    for i in 2:d
        sol = nlsolve(R, gamma)
        gamma = sol.zero
        if i != d
            gamma = [gamma zeros(7)]        
        end
    end
    return gamma
end


EEE = function(d;) 
    gamma = cheby(d)
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
    z_grid = grid_z
    d = Int(length(gamma)/7)
    k0 = knorm
    c0 = (reduce(hcat, [cons.((gamma,), k0, i, d) for i in 1:7]))'
    k1level = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')'
    k1norm = (k1level .- k_ss)/((k .- k_ss)[500])
    resid = zeros(7,500)
    for state in 1:7
        for w in 1:500
            c1 = cons.((gamma,), k1norm[state,w], 1:7, d) 
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end

    return resid'
end

plot(k, EEE(6))


#####################################
############# Questão 2 #############
#####################################
