using Plots, Distributions, Optim, NLsolve, BenchmarkTools, ProfileView, Distributed# Pacotes que estou usando

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
conscheb = function(gamma, k, z; z_len = z_len::Int64)
    d = Int.(length(gamma)/z_len)
    cons = zeros(d)
    for i in 1:d
        cons[i] = cheb(k, (i-1))
    end  
    return gamma[z,:]'cons
end

# Tentativa de fazer a função de residuos via list comprehension. Tá tudo certo até computar os resíduos
    # Função que retorna os resíduos dado a matriz de gammas, o capital, o estado e o numero de polinomios de chebyschev que estamos usando
    #resid = function(gamma; knorm = knorm, k = grid_k, p_z = p_z, beta = beta, z_grid = grid_z, alpha = alpha, delta = delta, k_len = k_len, z_len = z_len)
    #    d = Int(length(gamma)/z_len)
    #    root = chebroots(d)
    #    # Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de chebyschev de grau d
    #    # Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
    #    k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
    #    k0level = k0*(k .- k_ss)[k_len].+k_ss
    #    
    #    # Aplica a função conscheb (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
    #    # (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
    #    c0 = reduce(hcat, [conscheb.((gamma,), k0, 1:z_len) for k0 in k0])
    #    k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
    #    k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[k_len]) for j in 1:d for i in 1:z_len ], z_len, d)
    #
    #    k1norm = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:z_len], z_len, d))[i, j]) for j in 1:d for i in 1:z_len], z_len, d))
    #    
    #    c1 = [conscheb.((gamma,), k1norm[i, j], 1:z_len) for j in 1:d for i in 1:z_len ]
    #    c1 = reshape(c1, z_len, d)
    #    resid = reshape([u_line(c0[i]) - beta*(u_line.(c1[i]).*z_grid*(k1level[i]^(alpha-1)*alpha + 1 - delta))'p_z[i % z_len != 0 ? i % z_len : z_len, :] for i in 1:(d*z_len)], z_len, d)
    #    return resid
    # end


# Função que retorna os resíduos
Rcheb = function(gamma; knorm = knorm::Array{Float64,1}, k = grid_k::Array{Float64,1}, p_z = p_z::Array{Float64,2},
    beta = beta::Float64, z_grid = grid_z::Array{Float64,1}, alpha = alpha::Float64,
    delta = delta::Float64, k_len = k_len::Int64, z_len = z_len::Int64)
    d = Int(length(gamma)/z_len) # Número de polinomios que vou usar
    root = chebroots(d) # Raízes do polinomio de chebyschev de grau d
    k0 = root
    k0level = k0*(k .- k_ss)[k_len].+k_ss
    c0 = zeros(z_len,d)
    k1level = zeros(z_len,d)
    k1norm = zeros(z_len,d)
    resid = zeros(z_len,d)
    c1 = zeros(z_len)
    @sync @distributed for state in 1:z_len
        for w in 1:d
            c0[state, w] = conscheb(gamma, k0[w], state)
            k1level[state, w] = z[state]*(k0level[w]^alpha) + (1 - delta)*k0level[w] - c0[state,w]
            k1norm[state, w] = (k1level[state, w] - k_ss)/((k .- k_ss)[k_len])
            c1 = conscheb.((gamma,), k1norm[state,w], 1:z_len) 
            resid[state, w] = u_line(c0[state, w]) - beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])
        end
    end
    return resid
end

cheby_gamma = function(d::Int64; z_len = z_len::Int64)
    gamma::Array{Float64,2} = ones(z_len,2)
    for i in 2:d
        sol = nlsolve(Rcheb, gamma)
        gamma = sol.zero
        if i != d
            gamma = [gamma zeros(z_len)]        
        end
    end
    return gamma
end


EEEcheb = function(d; z_len = z_len::Int64) 
    gamma = cheby_gamma(d)
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
    z_grid = grid_z
    k0 = knorm
    c0 = (reduce(hcat, [conscheb.((gamma,), k0, i) for i in 1:z_len]))'
    k1level = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')'
    k1norm = (k1level .- k_ss)/((k .- k_ss)[k_len])
    resid = zeros(z_len,k_len)
    @sync @distributed for state in 1:z_len
        for w in 1:k_len
            c1 = conscheb.((gamma,), k1norm[state,w], 1:z_len) 
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end

    return Array(resid')::Array{Float64,2}, Array(c0')::Array{Float64,2}
end


sol_cheb = @time EEEcheb(6)

plot(k, EEEcheb(6)[1])




#####################################
############# Questão 2 #############
#####################################



############# Colocação #############
A = Array(reshape(LinRange(2,4, 7*11), 7, 11))
d = Int(length(A)/7)
# função phi; i é o indice da função, k é o capital que vamos aplicar a função, n é o numero de funções que vamos usar 
index = Int.(floor.(LinRange(1, length(grid_k), d)))
k_col = grid_k[index]

phi = function(i, k; k_col = k_col)
    if i == 1
        if k_col[i] <= k && k <= k_col[i+1]
            func = (k_col[i+1] - k)/(k_col[i+1] - k_col[i])
        else
             k < k_col[i] ? func = 1.0 : func = 0.0
        end

    elseif i == length(k_col)
        if k_col[i-1] <= k && k <= k_col[i]
            func = (k - k_col[i-1])/(k_col[i] - k_col[i-1])
        else
            k > k_col[i] ? func = 1.0 : func = 0.0      
        end

    else 
        if k_col[i-1] <= k && k <= k_col[i]
            func = (k - k_col[i-1])/(k_col[i] - k_col[i-1])
        elseif k_col[i] <= k && k <= k_col[i+1]
            func = (k_col[i+1] - k)/(k_col[i+1] - k_col[i])
        else 
            func = 0.0
        end
    end
    return func
end

conscol = function(A, k, z; z_len = z_len::Int64)
    d::Int64 = Int(length(A)/z_len)
    cons = zeros(d)
    for i in 1:d
        cons[i] = phi(i, k)
    end  
    return A[z,:]'cons
end


Rcol = function(A; k = grid_k::Array{Float64,1}, p_z = p_z::Array{Float64,2}, beta = beta::Float64,
    z_grid = grid_z::Array{Float64,1}, alpha = alpha::Float64, delta = delta::Float64,
    k_len = k_len::Int64, z_len = z_len::Int64)
    
    d = Int(length(A)/z_len) # Número de polinomios que vou usar
    index = Int.(floor.(LinRange(1, length(k), d)))
    k0 = k[index]
    c0 = zeros(z_len, d)
    k1 = zeros(z_len, d)
    resid = zeros(z_len,d)
    c1 = zeros(z_len)
    @sync @distributed for state in 1:z_len
        for w in 1:d
            c0[state, w] = conscol(A, k0[w], state)
            k1[state, w] = z[state]*(k0[w]^alpha) + (1 - delta)*k0[w] - c0[state,w]

            c1 = conscol.((A,), k1[state,w], 1:7) 

            resid[state, w] = u_line(c0[state, w]) - beta*(((u_line.(c1)).*((z_grid*(k1[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])
        end
    end
    return resid
end


EEE2 = function(A; k0 = grid_k::Array{Float64,1}, z_grid = grid_z::Array{Float64,1}, z_len = z_len::Int64) 

    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
    c0 = (reduce(hcat, [conscol.((A,), k0, i) for i in 1:z_len]))'
    k1 = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')'
    resid = zeros(z_len,k_len)
    @sync @distributed for state in 1:z_len
        for w in 1:k_len
            c1 = conscol.((A,), k1[state,w], 1:z_len) 
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end

    return resid'::Array{Float64,1}, c0'::Array{Float64,1}
end

@time nlsolve(Rcol, A)

sol_col = ProfileView.@profview nlsolve(Rcol, A)
A_col = sol_col.zero
EEE_col = EEE2(A_col)

plot(k, EEE_col[1])



############# Galerkin #############
phi_low = function(i::Int64, k::Float64; k_col = k_col::Array{Float64, 1})
    if i == 1
        k < k_col[i] ? func = 1.0 : func = 0.0

    elseif i == length(k_col)
        if k_col[i-1] <= k && k <= k_col[i]
            func = (k - k_col[i-1])/(k_col[i] - k_col[i-1])
        else 
            func = 0.0
        end

    else 
        if k_col[i-1] <= k && k <= k_col[i]
            func = (k - k_col[i-1])/(k_col[i] - k_col[i-1])
        else 
            func = 0.0
        end
    end
    return func
end



phi_high = function(i::Int64, k::Float64; k_col = k_col::Array{Float64, 1})

    if i == 1
        if k_col[i] <= k && k <= k_col[i+1]
            func = (k_col[i+1] - k)/(k_col[i+1] - k_col[i])
        else
            func = 0.0
        end

    elseif i == length(k_col)
            k > k_col[i] ? func = 1.0 : func = 0.0          

    else 
        if k_col[i] <= k && k <= k_col[i+1]
            func = (k_col[i+1] - k)/(k_col[i+1] - k_col[i])
        else 
            func = 0.0
        end
    end
    return func
end


Rcol2 = function(A::Array{Float64, 2}, k::Float64, z::Int64; 
    p_z = p_z::Array{Float64, 2}, beta = beta::Float64, z_grid = grid_z::Array{Float64, 1},
    alpha = alpha::Float64, delta = delta::Float64, k_len = k_len::Int64, z_len = z_len::Int64)

    c0 = conscol(A, k, z)
    k1 = grid_z[z]*(k^alpha) + (1 - delta)*k - c0
    c1 = conscol.((A,), k1, 1:z_len) 
    resid::Float64 = u_line(c0) - beta*(((u_line.(c1)).*((z_grid*(k1^(alpha-1)*alpha)) .+ (1-delta)))'p_z[z,:])
    return resid
end


gausscheb = function(A; nr = 15, z = grid_z::Array{Float64, 1}, z_len = z_len::Int64, k_col = k_col::Array{Float64, 1})
    d = Int(length(A)/z_len)
    roots = chebroots(nr)
    f_low = zeros(nr)
    f_high = zeros(nr)
    integral = zeros(z_len, d)
    

    @sync @distributed for state in 1:length(z)
        for i in 1:d
            if i == 1
                b = k_col[i+1]
                a = k_col[i]
                for n in 1:nr
                    k = (roots[n] + 1)*(b - a)/2 + a
                    f_high[n] = Rcol2(A, k, state)*phi_high(i, k)*(sqrt((1 - roots[n]^2)))
                end
                integral[state, i] = pi*(b - a)/(2*nr)*sum(f_high)
            elseif i == d
                b = k_col[i]
                a = k_col[i-1]
                for n in 1:nr
                    k = (roots[n] + 1)*(b - a)/2 + a
                    f_low[n] = Rcol2(A, k, state)*phi_low(i, k)*(sqrt((1 - roots[n]^2)))
                end
                integral[state, i] = pi*(b - a)/(2*nr)*sum(f_low)
            else
                b_low = k_col[i]
                a_low = k_col[i-1]
                b_high = k_col[i+1]
                a_high = b_low
                for n in 1:nr
                    k_low = (roots[n] + 1)*(b_low - a_low)/2 + a_low
                    k_high = (roots[n] + 1)*(b_high - a_high)/2 + a_high
                    f_low[n] = Rcol2(A, k_low, state)*phi_low(i, k_low)*(sqrt((1 - roots[n]^2)))
                    f_high[n] = Rcol2(A, k_high, state)*phi_high(i, k_high)*(sqrt((1 - roots[n]^2)))
                end
                integral[state, i] = pi*(b_high - a_high)/(2*nr)*sum(f_high) + pi*(b_low - a_low)/(2*nr)*sum(f_low) 
            end
        end
    end
    return integral
end


sol_galerkin = @time nlsolve(gausscheb, A)

A_galerkin = test.zero
plot(EEE2(A_galerkin)[1])