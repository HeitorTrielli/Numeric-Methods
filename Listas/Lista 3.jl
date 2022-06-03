using Plots, Distributions, NLsolve, BenchmarkTools, Distributed # Pacotes que estou usando
Threads.nthreads() # Quantas threads estamos usando
# Usando código das listas passadas
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

        return (prob = probs::Array{Float64,2}, grid = Array(grid)::Array{Float64,1})
    end;
    # Definindo as variáveis
    beta, mu, alpha, delta, rho, sigma = 0.987, 2.0, 1/3, 0.012, 0.95, 0.007;

    # Definindo o capital de steady state
    k_ss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));
    z_len = 7;
    k_len = 500;
    grid_z = exp.(tauchen(z_len).grid); # Valores que exp(z_t) pode tomar 
    grid_k = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len));
    z = copy(grid_z); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len).prob; # Matriz de transição de z
    k = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len)); # Grid de capital
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
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

# Os gráficos e os tempos estão todos chamados no fim do código #

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

# Função que monta o consumo em função de gamma, k e d polinômios de Chebyschev. d é definido pelo tamanho de gamma
conscheb = function(gamma, k, z; z_len = z_len::Int64)
    d = Int.(length(gamma)/z_len)
    cons = zeros(d)
    for i in 1:d
        cons[i] = cheb(k, (i-1))
    end  
    return gamma[z,:]'cons
end

#= Tentativa de fazer a função de residuos via list comprehension. Tá tudo certo até computar os resíduos
    Função que retorna os resíduos dado a matriz de gammas, o capital, o estado e o numero de polinomios de Chebyschev que estamos usando
    resid = function(gamma; knorm = knorm, k = grid_k, p_z = p_z, beta = beta, z_grid = grid_z, alpha = alpha, delta = delta, k_len = k_len, z_len = z_len)
        d = Int(length(gamma)/z_len)
        root = chebroots(d)
        # Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de Chebyschev de grau d
        # Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
        k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
        k0level = k0*(k .- k_ss)[k_len].+k_ss
        
        # Aplica a função conscheb (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
        # (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
        c0 = reduce(hcat, [conscheb.((gamma,), k0, 1:z_len) for k0 in k0])
        k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
        k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[k_len]) for j in 1:d for i in 1:z_len ], z_len, d)
    
        k1norm = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:z_len], z_len, d))[i, j]) for j in 1:d for i in 1:z_len], z_len, d))
        
        c1 = [conscheb.((gamma,), k1norm[i, j], 1:z_len) for j in 1:d for i in 1:z_len ]
        c1 = reshape(c1, z_len, d)
        resid = reshape([u_line(c0[i]) - beta*(u_line.(c1[i]).*z_grid*(k1level[i]^(alpha-1)*alpha + 1 - delta))'p_z[i % z_len != 0 ? i % z_len : z_len, :] for i in 1:(d*z_len)], z_len, d)
        return resid
     end
=#

# Função que retorna os resíduos dado gamma. Já leva em conta k e z. 
Rcheb = function(gamma; knorm = knorm::Array{Float64,1}, k = grid_k::Array{Float64,1}, p_z = p_z::Array{Float64,2},
    beta = beta::Float64, z_grid = grid_z::Array{Float64,1}, alpha = alpha::Float64,
    delta = delta::Float64, k_len = k_len::Int64, z_len = z_len::Int64)
    d = Int(length(gamma)/z_len) # Número de polinomios que vou usar
    k0 = chebroots(d) # Vamos fazer R ser 0 nas raízes do polinômio de Chebyschev de grau d
    k0level = k0*(k .- k_ss)[k_len].+k_ss # Trazendo as raízes pro grid de capital
    c0, k1level, k1norm, resid, c1, = zeros(z_len,d), zeros(z_len,d), zeros(z_len,d), zeros(z_len,d), zeros(z_len,d) # Pré-alocando memória
    #= Primeiro achamso c_chapeu em k0, depois usamos esse valor para achar k1
    depois normalizamos k1 para poder jogar na função c_chapeu em cada estado da natureza
    e depois calculamos os resíduos. @sync @distributed é para paralelizar o for=#
    @sync @distributed for state in 1:z_len
        for w in 1:d
            c0[state, w] = conscheb(gamma, k0[w], state) 
            k1level[state, w] = z[state]*(k0level[w]^alpha) + (1 - delta)*k0level[w] - c0[state,w] 
            k1norm[state, w] = (k1level[state, w] - k_ss)/((k .- k_ss)[k_len]) 
            c1 = conscheb.((gamma,), k1norm[state,w], 1:z_len) 
            resid[state, w] = u_line(c0[state, w]) - beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:]) 
        end # for state
    end # for w
    return resid
end

# Função que acha o gamma ótimo dado o numero de polinômios de Chebyschev que queremos usar, no caso d.
cheby_gamma = function(d::Int64; z_len = z_len::Int64)

    gamma::Array{Float64,2} = ones(z_len,2)
    for i in 2:d
        sol = nlsolve(Rcheb, gamma) # nlsolve acha o zero de um sistema de equações não linear
        gamma = sol.zero
        if i != d # Se ainda não temos uma matriz 7xd então continua aumentando o gamma que vai ser chutado, dado o gamma ótimo da ultima iteração
            gamma = [gamma zeros(z_len)]        
        end
    end
    return gamma
end

# Função que acha os erros de euler (e ja acha o consumo também) dado o numero de polinômios de Chebyschev que a gente quer
solve_cheb = function(d; z_len = z_len::Int64, knorm = knorm::Array{Float64, 1}, k_len = k_len::Int64, z_grid = grid_z::Array{Float64, 1}) 
    gamma = cheby_gamma(d)
    zmat = repeat(z,1,k_len)' #matriz com z repetidos nas linhas
    kmat = repeat(k,1,z_len)  #matriz com k repetido nas colunas
    k0 = knorm
    c0 = (reduce(hcat, [conscheb.((gamma,), k0, i) for i in 1:z_len]))' # Aplica conscheb(gamma) em cada k normalizado, em cada estado
    k1level = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')' # Acha todos os k1, dado a matriz de c0
    k1norm = (k1level .- k_ss)/((k .- k_ss)[k_len]) # Normaliza k1
    resid = zeros(z_len,k_len) #Pré alocando
    @sync @distributed for state in 1:z_len
        for w in 1:k_len
            c1 = conscheb.((gamma,), k1norm[state,w], 1:z_len) # Aplica conscheb em k1 em todos os estados para achar c1
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1level[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end
    return (euler = Array(resid')::Array{Float64,2}, c = Array(c0')::Array{Float64,2}, pol = Array(k1level')::Array{Float64,2})
end


#####################################
############# Questão 2 #############
#####################################

d = 11 # Numero de pontos limite que vamos usar
cons_level = zmat.*(kmat.^alpha) - delta*kmat # Sei que o consumo não deve fugir muito dessa matriz, então vou usar ela para guiar o chute inicial
low = minimum(cons_level)
high = maximum(cons_level)
A = Array(reshape(LinRange(low,high, z_len*d), z_len, d)) # Chute inicial para os proximos métodos
index = Int.(floor.(LinRange(1, length(grid_k), d))) # Achando os índices dos capitais que vão servir como pontos limite
k_col = grid_k[index] # Escolhendo esses capitais
nr = 15 # Numero de raízes de Chebyschev que vão 

############# Colocação #############
# função phi; i é o indice da função, k é o capital que vamos aplicar a função, n é o numero de funções que vamos usar 
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

# função que retorna sum(A(z,i)*phi(i, k))
conscol = function(A, k, z; z_len = z_len::Int64, d = d::Int64)
    cons = zeros(d)
    for i in 1:d
        cons[i] = phi(i, k)
    end  
    return A[z,:]'cons
end

# Acha os resíduos. Análogo a versão de Chebyschev, mas adaptando para o novo consumo, e sem a parte de normalizar o capital
Rcol = function(A; p_z = p_z::Array{Float64,2}, beta = beta::Float64,
    z_grid = grid_z::Array{Float64,1}, alpha = alpha::Float64, delta = delta::Float64,
    k_len = k_len::Int64, z_len = z_len::Int64, d = d::Int64, k_col = k_col::Array{Float64,1})

    k0 = k_col # Vamos calcular os resíduos nos pontos limite para poder zerar
    c0, k1, resid = zeros(z_len, d), zeros(z_len, d), zeros(z_len, d) # Pré-alocando
    c1 = zeros(z_len) # Pré-alocando
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

# Função que 1.) acha o A ótimo, dado o chute inicial A e 2.) calcula e retorna c0, a policy e os erros de euler
solve_col = function(A; k0 = grid_k::Array{Float64,1}, z_grid = grid_z::Array{Float64,1}, z_len = z_len::Int64) 
    solution = nlsolve(Rcol, A)
    A_star = solution.zero
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
    c0 = (reduce(hcat, [conscol.((A_star,), k0, i) for i in 1:z_len]))'
    k1 = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')'
    resid = zeros(z_len,k_len)
    @sync @distributed for state in 1:z_len
        for w in 1:k_len
            c1 = conscol.((A_star,), k1[state,w], 1:z_len) 
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end

    return (euler = Array(resid')::Array{Float64,2}, c = Array(c0')::Array{Float64,2}, pol = Array(k1')::Array{Float64, 2})
end


############# Galerkin #############
# Análogo a phi, mas apenas levando em conta a parte de k in [k_{i-1}, k_i]
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

# Análogo a phi, mas apenas levando em conta a parte de k in [k_i, k_{i+1}]
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

#= Análogo as funções de resíduos anteriores, mas dessa vez retorna um escalar.
As outras aplicam essa função, mas em todo k e todo z e retornam uma matriz =#
Rgal = function(A::Array{Float64, 2}, k::Float64, z::Int64; 
    p_z = p_z::Array{Float64, 2}, beta = beta::Float64, z_grid = grid_z::Array{Float64, 1},
    alpha = alpha::Float64, delta = delta::Float64, k_len = k_len::Int64, z_len = z_len::Int64)

    c0 = conscol(A, k, z)
    k1 = grid_z[z]*(k^alpha) + (1 - delta)*k - c0
    c1 = conscol.((A,), k1, 1:z_len) 
    resid::Float64 = u_line(c0) - beta*(((u_line.(c1)).*((z_grid*(k1^(alpha-1)*alpha)) .+ (1-delta)))'p_z[z,:])
    return resid
end

# Aplica a integral de galerkin usando a quadratura de Gauss-Chebyschev
gausscheb = function(A; nr = nr, z = grid_z::Array{Float64, 1}, z_len = z_len::Int64, k_col = k_col::Array{Float64, 1}, d = d::Int64)
    roots = chebroots(nr) # Raízes do polinomio de chebyschev de grau nr, a ser usado na quadratura
    f_low, f_high = zeros(nr), zeros(nr) #pré-alocando
    integral = zeros(z_len, d) #pré-alocando

    @sync @distributed for state in 1:length(z)
        for i in 1:d
            if i == 1
                b = k_col[i+1] # valor superior da integral
                a = k_col[i] #valor inferior da integral
                for n in 1:nr
                    k = (roots[n] + 1)*(b - a)/2 + a # trazendo a raíz para o nivel do capital
                    f_high[n] = Rgal(A, k, state)*phi_high(i, k)*(sqrt((1 - roots[n]^2))) # calculando os termos do somatório
                end
                integral[state, i] = pi*(b - a)/(2*nr)*sum(f_high) # calculando a integral
            elseif i == d # análogo ao anterior, mas considerando apenas os pontos abaixo de k[d]
                b = k_col[i]
                a = k_col[i-1]
                for n in 1:nr
                    k = (roots[n] + 1)*(b - a)/2 + a
                    f_low[n] = Rgal(A, k, state)*phi_low(i, k)*(sqrt((1 - roots[n]^2)))
                end
                integral[state, i] = pi*(b - a)/(2*nr)*sum(f_low)
            else # análogo aos anteriores
                b_low = k_col[i]
                a_low = k_col[i-1]
                b_high = k_col[i+1]
                a_high = b_low
                for n in 1:nr
                    k_low = (roots[n] + 1)*(b_low - a_low)/2 + a_low
                    k_high = (roots[n] + 1)*(b_high - a_high)/2 + a_high
                    f_low[n] = Rgal(A, k_low, state)*phi_low(i, k_low)*(sqrt((1 - roots[n]^2)))
                    f_high[n] = Rgal(A, k_high, state)*phi_high(i, k_high)*(sqrt((1 - roots[n]^2)))
                end
                integral[state, i] = pi*(b_high - a_high)/(2*nr)*sum(f_high) + pi*(b_low - a_low)/(2*nr)*sum(f_low) 
            end
        end
    end
    return integral
end

# Função igual solve_gal, apenas mudando o A, que antes vinha de forçar Rcol = 0, agora vem de forçar gausscheb = 0
solve_gal = function(A; k0 = grid_k::Array{Float64,1}, z_grid = grid_z::Array{Float64,1}, z_len = z_len::Int64) 
    solution = nlsolve(gausscheb, A)
    A_star = solution.zero
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)  
    c0 = (reduce(hcat, [conscol.((A_star,), k0, i) for i in 1:z_len]))'
    k1 = (zmat.*(kmat.^alpha) + (1-delta)*kmat - c0')'
    resid = zeros(z_len,k_len)
    @sync @distributed for state in 1:z_len
        for w in 1:k_len
            c1 = conscol.((A_star,), k1[state,w], 1:z_len) 
            resid[state, w] = log10(abs(1 - (u_inv(beta*(((u_line.(c1)).*((z_grid*(k1[state, w]^(alpha-1)*alpha)) .+ (1-delta)))'p_z[state,:])))/c0[state,w]))
        end
    end

    return (euler = Array(resid')::Array{Float64,2}, c = Array(c0')::Array{Float64,2}, pol = Array(k1')::Array{Float64, 2})
end


######################################
############# Resultados #############
######################################
#= Código para rodar o benchmark. Demora, por isso botei em comentário
ben_cheb = @benchmark solve_cheb(6) seconds = 60
ben_col = @benchmark solve_col(A) seconds = 60
ben_gal = @benchmark solve_gal(A) seconds = 300 =#

sol_cheb = @time solve_cheb(6)
sol_col = @time solve_col(A)
sol_gal = @time solve_gal(A)

#= Função que plotta os 1.) Erros de Euler; 2.) Função Consumo; 3.) Função política; 4.) Função política nos primeiros 100 capitais.
Caso você chame a função na variavel "cheb" ela retorna os itens acima para o metodo de colocação com polinomios de Chebyschev, caso
você chame a função em "col" ela te retorna os itens acima para o metodo da colocação elementos finitos com as funções phi, 
e qualquer outro input retorna os itens para elementos finitos com o metodo de Galerkin nas funções phi. =#
plot_cke = function(method::String; cheb = sol_cheb, col = sol_col, gal = sol_gal)
    labels = ["State 1" "State 2" "State 3" "State 4" "State 5" "State 6" "State 7"]
    if method == "cheb"
        euler = plot(k, cheb.euler, label = labels, title = "Erros de Euler \n Chebyschev", legend = :bottomleft, xlabel = "Capital")
        c = plot(k, cheb.c, label = labels, title = "Função Consumo \n Chebyschev", legend = :topleft, xlabel = "Capital")
        pol = plot(k, cheb.pol, label = labels, title = "Função Política \n Chebyschev", legend = :topleft, xlabel = "Capital")
        pol_zoom = plot(k[1:100], cheb.pol[1:100, :], label = labels, title = "Função Política Zoom \n Chebyschev", legend = :topleft, xlabel = "Capital")
    elseif method == "col"
        euler = plot(k, col.euler, label = labels, title = "Erros de Euler \n Collocation", legend = :bottomleft, xlabel = "Capital")
        c = plot(k, col.c, label = labels, title = "Função Consumo \n Collocation", legend = :topleft, xlabel = "Capital")
        pol = plot(k, col.pol, label = labels, title = "Função Política \n Collocation", legend = :topleft, xlabel = "Capital")
        pol_zoom = plot(k[1:100], col.pol[1:100, :], label = labels, title = "Função Política Zoom \n Collocation", legend = :topleft, xlabel = "Capital")
    else
        euler = plot(k, gal.euler, label = labels, title = "Erros de Euler \n Galerkin", legend = :bottomleft, xlabel = "Capital")
        c = plot(k, gal.c, label = labels, title = "Função Consumo \n Galerkin", legend = :topleft, xlabel = "Capital")
        pol = plot(k, gal.pol, label = labels, title = "Função Política \n Galerkin", legend = :topleft, xlabel = "Capital")
        pol_zoom = plot(k[1:100], gal.pol[1:100, :], label = labels, title = "Função Política Zoom \n Galerkin", legend = :topleft, xlabel = "Capital")
    end
    return (euler = euler, c = c, pol = pol, pol_zoom = pol_zoom)
end


#= Função que vai plottar quantos % um método desviou do outro 1.) nos erros de euler; 2.) na função consumo, 3.) na função política.
Como input você deve escolher dois métodos dentre: "gal", "col" ou "cheb" (tem que chamar a função em strings). Para comparar Galerkin
com a colocação tanto faz quais strings você coloca. A ordem não importa em nenhum dos casos. O título do grafico mostra quem é tratado
como numerador e quem é tratado como denominador. =#
plot_compare = function(method_1::String, method_2::String;cheb = sol_cheb, col = sol_col, gal = sol_gal)
    labels = ["State 1" "State 2" "State 3" "State 4" "State 5" "State 6" "State 7"]
    if (method_1 == "cheb" && method_2 == "col") || (method_2 == "cheb" && method_1 == "col")
        euler = plot(k, (cheb.euler ./ gal.euler) .- 1, label = labels, title = "Erros de Euler \n Chebyschev / Colocation")
        c = plot(k, (cheb.c ./ gal.c) .- 1, label = labels, legend = :bottomright, title = "Função Consumo \n Chebyschev / Colocation")
        pol = plot(k, (cheb.pol ./ gal.pol) .- 1, label = labels,title = "Função Política \n Chebyschev / Colocation")

    elseif (method_1 == "cheb" && method_2 == "gal") || (method_2 == "cheb" && method_1 == "gal")
        euler = plot(k, (cheb.euler ./ col.euler) .- 1, label = labels, legend = :bottomleft, title = "Erros de Euler \n Chebyschev / Galerkin")
        c = plot(k, (cheb.c ./ col.c) .- 1, label = labels, title = "Função Consumo \n Chebyschev / Galerkin")
        pol = plot(k, (cheb.pol ./ col.pol) .- 1, label = labels, legend = :bottomright, title = "Função Política \n Chebyschev / Galerkin")

    else
        euler = plot(k, (gal.euler ./ col.euler) .- 1, label = labels, legend = :bottomright, title = "Erros de Euler \n Galerkin / Colocação")
        c = plot(k, (gal.c ./ col.c) .- 1, label = labels, title = "Função Consumo \n Galerkin / Colocação")
        pol = plot(k, (gal.pol ./ col.pol) .- 1, label = labels, legend = :bottomright, title = "Função Politica \n Galerkin / Colocação")
    
    end
    return (euler = euler, c = c, pol = pol)
end

plots_cheb = plot_cke("cheb")
plots_cheb.euler
plots_cheb.c
plots_cheb.pol
plots_cheb.pol_zoom

plots_col = plot_cke("col")
plots_col.euler
plots_col.c
plots_col.pol
plots_col.pol_zoom


plots_gal = plot_cke("gal")
plots_gal.euler
plots_gal.c
plots_gal.pol
plots_gal.pol_zoom


cheb_gal = plot_compare("cheb", "gal")
cheb_gal.euler
cheb_gal.c
cheb_gal.pol

cheb_col = plot_compare("cheb", "col")
cheb_col.euler
cheb_col.c
cheb_col.pol

gal_col = plot_compare("gal", "col")
gal_col.euler
gal_col.c
gal_col.pol