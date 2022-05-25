using Plots, BenchmarkTools, Distributions, Distributed, ProfileView, Roots, Dierckx # Pacotes que estou usando

Threads.nthreads()
# Definindo funções e variáveis que vou usar no futuro
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
    zmat = repeat(grid_z,1,k_len)';
    kmat = repeat(grid_k,1,z_len);
    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len)); # Grid de capital
    c0 = zmat.*(kmat.^alpha) - delta*kmat;
    

    ############# Utilidade #############
    utility = function(c::Float64; mu::Float64 = 2.0)
    return (c^(1 - mu) - 1)/(1 - mu)
    end;
    v0 = utility.(c0);

    u_line = function(c, mu = 2.0)
        c^(-mu)
    end
    
    u_inv = function(c, mu = 2.0)
        c^(-1/mu)
    end
# 

############ Erros de Euler #############
# Função que nos retorna os erros de Euler
EEE = function(x; k_len = 500, z_len = 7, k_grid = grid_k, probs = p_z, mu = 2.0, z = grid_z, alpha = 1/3, delta = 0.012)
    c = x[3]
    policy = x[2]
    euler = zeros(k_len, z_len)
    for state in 1:z_len
        for capital in 1:k_len
            index = findall(policy[capital, state] .== k_grid)[1] ## Índice para saber qual linha da matriz de consumo levar em conta na esperança
            c_line = c[index,:] ## Consumo dado a escolha
            E = (u_line.(c_line).*(alpha*z*(policy[capital,state]^(alpha-1)) .+ (1-delta)))'probs[state,:] # Esperança
            euler[capital, state] = log10(abs(1 - u_inv(beta*E) / c[capital,state]))  # Erro de euler
        end
    end 
    return euler
end


#########################################
############### Questão 3 ###############
#########################################

############## Brute Force ##############
value_function_brute = function(;v_0::Array{Float64} = v0, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss, k_len::Int64 = 500, z_len::Int64 = 7)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital

    # Pré-alocando os valores da função política, da atualização da função valor e a função valor
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0 
    error = 1.0
    s = 1
    while error > tol # Checa se convergiu
        @sync @distributed for state in 1:z_len #@sync @distributed faz a parelização sincronizada
            for capital in 1:k_len  
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0] # Escolhendo os capitais onde o consumo é positivo
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :] # Escolhendo as linhas onde os capitais são positivos
                # val é um vetor com a função valor aplicada em cada possivel escolha para k'
                val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :]
                s = argmax(val) # Acha o índice do máximo
                v_next[capital, state] = val[s] # Define o valor como o máximo
                k_line[capital, state] = k_possible[s] # Define o capital como o capital do índice do máximo
            end # for k
        end # for z
        error = maximum(abs.(value - v_next)) # Checa o erro
        value = copy((1 - inert)*v_next + inert*value)
    end # while

    # Aqui acho a matriz de consumo de forma matricial
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)   
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat
    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end; # function


####### Exploiting monotonicity #########
value_function_monotone = function(;v_0::Array{Float64} = v0, k_len::Int64 = 500, z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital

    # Pré-alocando os valores da função política, da atualização da função valor e a função valor
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
    error = 1.0

    s=1 # s vai guiar a partir de onde devemos olhar pela monotonicidade
    while error > tol # Checa se convergiu 
        @sync @distributed for state in 1:z_len 
            s = 1
            for capital in 1:k_len
                # Vamos olhar apenas para os k que deixam o consumo positivo e que são maiores que o k escolhido no estado de capital anterior
                # Restante análogo à força bruta
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s]] 
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s], :] 
                val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :]
                s = argmax(val)
                v_next[capital, state] = val[s]
                k_line[capital, state] = k_possible[s]
            end # for k
        end # for z
        error = maximum(abs.(value - v_next))
        value = copy(v_next)
       
    end

    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat # Análogo a força bruta

    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end # function


############### Concavity ###############
# Função que vai retornar o ponto maximo da função, o seu respectivo k' e o índice desse ponto, aproveitando a concavidade da função que gerou os pontos do vetor
lastpos = function(val::Array{Float64}, k_possible::Array{Float64}) 
    v = 0.0
    k = 0.0
    index = 0
    for i in 2:length(val)
        if val[i] < val[i-1]
            v = val[i-1]
            k = k_possible[i-1]
            index = i-1
            break
        else
            v = val[i]
            k = k_possible[i]
            index = i
        end
    end
    return v::Float64, k::Float64, index::Int64
end;

value_function_concave = function(;v_0::Array{Float64} = v0, k_len::Int64 = 500, z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital

    # Pré-alocando os valores da função política, da atualização da função valor e a função valor
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
    error = 1.0
    while error > tol # Checa se convergiu 
        @sync @distributed for state in 1:z_len
            for capital in 1:k_len  
                k_possible = k[(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k) .> 0] # Análogo a força bruta 
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :] # Análogo a força bruta      
                val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :] # Análogo a força bruta
                v_next[capital, state], k_line[capital, state] = lastpos(val, k_possible)[1:2] # Usa a função que defini acima para achar o máximo
            end # for k
        end # for z
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
    end # while
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat # Análogo a força bruta

    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end; # function


####### Concavity + monotonicity ########
# Essa função é basicamente a mistura das ultimas duas, então não vou documentar muito
value_function = function(;v_0::Array{Float64} = v0, k_len::Int64 = 500, z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital


    # Pré-alocando
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0

    error = 1.0
    s=1 # Análogo a monotona
    while error > tol
        @sync @distributed for state in 1:z_len 
            s=1
            for capital in 1:k_len
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s]] # Análogo a monotona
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s], :] # Análogo a monotona              
                val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :] # Análogo a monotona
                v_next[capital, state], k_line[capital, state], s = lastpos(val, k_possible) # Análogo à concava
            end # for k
        end # for z
        
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
    end # while

    # Análogo a força bruta
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat 

    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end # function


#########################################
############### Questão 4 ###############
#########################################

############## Accelerator ##############
value_function_accelerator = function(;v_0::Array{Float64} = v0, k_len::Int64 = 500, z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    reset::Int64 = 10, start::Int64 = 30, kss::Float64 = k_ss)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital
    # Pré alocando
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0
    error = 1.0
    count = 0
    s = 1 # Análogo a monotona
    while error > tol

        @sync @distributed for state in 1:z_len    
            s = 1
            for capital in 1:k_len
                if count <= start || count%reset == 0 # Checa se estamos nas primeiras 30 iterações ou se estamos num multiplo de 10. Caso positivo realiza o 
                    # operador de maximização, caso contrario apenas atualiza a função valor, tomando a policy como dada
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s]] # Análogo à monotona
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0 .&& k .>= k[s], :] # Análogo à monotona 
                    val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :] # Análogo a monotona
                    v_next[capital, state], k_line[capital, state], s = lastpos(val, k_possible) # Análogo a concava
                else
                    v_next[capital, state] = utility(z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]) + beta*value[findall(k_line[capital,state] .== k)[1],:]'*p_z[state, :]
                end
            end # for k
        end # for z
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
        count += 1
    end
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat
    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}

end # function


## Accelerator sem usar monotonicidade e concavidade##
# Exatamente igual a função anterior, mas sem usar os passos de monotonicidade nem concavidade
value_function_accelerator_brute = function(;v_0::Array{Float64} = v0, k_len::Int64 = 500, z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    reset::Int64 = 10, start::Int64 = 30, kss::Float64 = k_ss)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital

    # Pré alocando
    k_line, v_next = zeros(k_len, z_len), zeros(k_len, z_len)
    value = v_0

    error = 1.0
    count = 0
    s = 1
    while error > tol
        @sync @distributed for state in 1:z_len    
            for capital in 1:k_len
                if count <= start || count%reset == 0
                    k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0] 
                    v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :]
                    val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :]
                    s = argmax(val)
                    v_next[capital, state] = val[s]
                    k_line[capital, state] = k_possible[s]
                else
                    v_next[capital, state] = utility(z[state]*(k[capital]^alpha) - k_line[capital, state] + (1 - delta)* k[capital]) + beta*value[findall(k_line[capital,state] .== k)[1],:]'*p_z[state, :]
                end
            end # for k
        end # for z
    
        error = maximum(abs.(value - v_next))
        value = copy((1 - inert)*v_next + inert*value)
        count += 1
    end

    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat

    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end # function


#########################################
############### Questão 5 ###############
#########################################

############## Multi grid ###############
# Função que faz a interpolação da função valor, mas não é uma interpolação perfeita. No finalzinho da interpolação tem uma imprecisão. É praticamente irrelevante.
lin_interpol = function(minsize::Int64, maxsize::Int64, value::Array{Float64}; z_len = 7::Int64)    
    
    step = Int(maxsize/minsize)
    v_1 = zeros(maxsize, z_len)
    Threads.@threads for i in 1:(maxsize - step)
        v_1[i, :] = value[Int(ceil(i/step)),:] + (i-1)%step*(value[Int(ceil((i+step)/step)),:] - value[Int(ceil(i/step)),:])/step
    end
    Threads.@threads for i in (maxsize - step + 1):maxsize
        v_1[i, :] = value[minsize, :] - (maxsize-i)%step*(value[Int(ceil((maxsize)/step)),:] - v_1[Int(ceil((maxsize-step))),:])/step
    end
    return v_1
end


value_function_mg = function(g1::Int64, g2::Int64, g3::Int64; z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)

    # Gerando o chute inicial para um grid de capital de 100 pontos
    k = Array(LinRange(0.75*kss, 1.25*kss, g1))
    z = tauchen(z_len)[2]
    zmat = repeat(z,1,g1)'
    kmat = repeat(k,1,z_len)
    v0 = utility.(zmat.*(kmat.^alpha) - delta*kmat)
    
    # Aplica a força bruta, adaptando para 100 pontos no grid
    value1 = value_function_brute(v_0 = v0, k_len = g1, tol = tol)[1]

    # Interpola a função valor obtida na força bruta de 100 pontos
    v1 = lin_interpol(g1, g2, value1)

    # Usa essa função valor interpolada como chute inicial
    value2 = value_function_brute(v_0 = v1, k_len = g2, tol = tol)[1]

    # Interpola a função valor anterior 500x7 para uma 5000x7
    v2 = lin_interpol(g2, g3, value2)

    # Usa essa nova função valor como chute inicial para a nova iteração com o grid de 5000 pontos
    value, k_line, c = value_function_brute(v_0 = v2, k_len = g3, tol = tol)

    return value, k_line, c
end # function

#########################################
############### Questão 6 ###############
#########################################

############ Grid Endógeno ##############
value_function_egm = function(;c0::Array{Float64,2} = c0, v0 = zeros(500,7), k::Array{Float64,1} = grid_k, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss, k_len::Int64 = 500, z_len::Int64 = 7)

    z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
    p_z = tauchen(z_len)[1]; # Matriz de transição de z
    k = Array(LinRange(0.75*kss, 1.25*kss, k_len)); # Grid de capital

    # Value function, policy function on c, policy function on k and variable for iteration
    a, k_next, c, value, cline = zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len), zeros(k_len, z_len)
    erro = 1.0
    while erro > tol
        @sync @distributed for state in 1:z_len
            zstate = z[state]
            prob = p_z[state,:]
            for capital in 1:k_len
                k_exo = k[capital]
                ck = c0[capital, :]
                prob = p_z[state,:]
                # Função que traz a CPO
                func = function(k;  zstate = zstate, k_exo = k_exo  , c = ck, prob = prob, z = grid_z, delta = 0.012, alpha = 1/3, beta =0.987)
                    E = (u_line.(c).*(alpha*z*(k_exo^(alpha-1)).+(1-delta)))'*prob
                    (1-delta)*k - k_exo - u_inv(beta*E) .+ zstate*(k^alpha)
                end
                # Acha os capitais que fazem a CPO ser zero. Isto é, acha o primeiro grid endógeno.
                a[capital, state] = fzero(func, a[capital, state])
            end # for k
            spl = Spline1D(a[:,state], k) # Função que vai interpolar k'(k_end, z) para k'(k, z)
            for capital in 1:k_len  
                k_next[capital, state] = spl(k[capital]) # Realizando a interpolação 
                c[capital, state] = z[state]*(k[capital]^alpha) + (1-delta)*k[capital] - k_next[capital, state] # Achando o consumo
            end # for k
        end # for z 
        erro = maximum(abs.(c - c0))
        c0 = copy(c)
    end
    
    # Temos o c(k,z), então interpolo para ter c'(k', z')
    for state in 1:z_len
        spl = Spline1D(k, c[:, state])
        cline[:, state] = spl(k_next[:, state])
    end  
    
    erro = 1
    # Agora iterando para achar a função valor

    while erro > 0.05
        
        @sync @distributed for state in 1:z_len
            spl = Spline1D(k, v0[:, state])
            v0[:, state] = spl(k_next[:, state]) # Temos v0[k, z], queremos v0[k', z], então interpolamos para achar v0
        end
        # Atualizando a função valor
        @sync @distributed for state in 1:z_len
            for capital in 1:k_len
                value[capital, state] = utility(c[capital, state]) + beta*v0[capital, :]'*p_z[state,:]
            end
        end
        erro = maximum(abs.(value - v0))
        v0 = copy((1 - inert)*value + inert*v0)
    end
    return value::Array{Float64,2}, k_next::Array{Float64,2}, c::Array{Float64,2}, a::Array{Float64,2}, cline::Array{Float64,2}
end; # function



# Função erros de euler revisitada para o caso do grid endógeno. Ela recebe c(k, z), k'(k, z) e c'(k', z') como input e calcula os erros de euler
EEE2 = function(c::Array{Float64,2}, policy::Array{Float64,2}, c_line::Array{Float64, 2}, k_len = 500, z_len = 7, 
    k_grid = grid_k, probs = p_z, mu = 2.0, z = grid_z, alpha = 1/3, delta = 0.012)
    euler = zeros(k_len, z_len)
    for state in 1:z_len
        for capital in 1:k_len
            E = (u_line.(c_line[capital, :]).*(alpha*z*(policy[capital]^(alpha-1)) .+ (1-delta)))'probs[state,:]
            euler[capital, state] = log10(abs(1 - u_inv(beta*E) / c[capital,state]))  
        end
    end 
    return euler
end

#########################################
############### Plotando ################
#########################################

labels = ["State 1" "State 2" "State 3" "State 4" "State 5" "State 6" "State 7"]

println("Brute Force")
brute_force = value_function_brute() 
plot(brute_force[2], brute_force[1], title = "Brute Force Value", legend = :topleft, label = labels)
plot(brute_force[2], title = "Brute Force Policy", legend = :topleft, label = labels)
plot(brute_force[2], brute_force[3], title = "Brute Force Consumption", legend = :topleft, label = labels)
plot(EEE(brute_force), title = "Brute Force", label = labels)



println("Monotone")
monotone = value_function_monotone() 
plot(monotone[2], monotone[1], title = "Monotone Value", legend = :topleft, label = labels)
plot(monotone[2], title = "Monotone Policy", legend = :topleft, label = labels)
plot(monotone[2], monotone[3], title = "Monotone Consumption", legend = :topleft, label = labels)
plot(EEE(monotone), title = "Monotone", label = labels)



println("Concavity")
concave =  value_function_concave() 
plot(concave[2], concave[1], title = "Concave Value", legend = :topleft, label = labels)
plot(concave[2], title = "Concave Policy", legend = :topleft, label = labels)
plot(concave[2], concave[3], title = "Concave Consumption", legend = :topleft, label = labels)
plot(EEE(concave), title = "Concave", label = labels)



println("Monotone + Concavity")
vf =  value_function() 
plot(vf[2], vf[1], title = "Monotone + Concave Value", legend = :topleft, label = labels)
plot(vf[2], title = "Monotone + Concave Policy", legend = :topleft, label = labels)
plot(vf[2], vf[3], title = "Monotone + Concave Consumption", legend = :topleft, label = labels)
plot(EEE(vf), title = "Monotone + Concave", label = labels)



println("Accelerator")
accelerator =  value_function_accelerator() 
plot(accelerator[2], accelerator[1], title = "Accelerator Value", legend = :topleft, label = labels)
plot(accelerator[2], title = "Accelerator Policy", legend = :topleft, label = labels)
plot(accelerator[2], accelerator[3], title = "Accelerator Consumption", legend = :topleft, label = labels)
plot(EEE(accelerator), title = "Accelerator", label = labels)



println("Accelerator Brute")
accel_brute =  value_function_accelerator_brute() 
plot(accel_brute[2], accel_brute[1], title = "Accelerator Brute Value", legend = :topleft, label = labels)
plot(accel_brute[2], title = "Accelerator Brute Policy", legend = :topleft, label = labels)
plot(accel_brute[2], accel_brute[3], title = "Accelerator Brute Consumption", legend = :topleft, label = labels)
plot(EEE(accel_brute), title = "Accelerator Brute", label = labels)



println("Multigrid")
multigrid = @time value_function_mg(100, 500, 5000) 
plot(multigrid[2], multigrid[1], title = "Multigrid Value", legend = :topleft, label = labels)
plot(multigrid[2], title = "Multigrid Policy", legend = :topleft, label = labels)
plot(multigrid[2], multigrid[3], title = "Multigrid Consumption", legend = :topleft, label = labels)
grid_kmg = Array(LinRange(0.75*k_ss, 1.25*k_ss, 5000));
plot(EEE(multigrid, k_len = 5000, k_grid = grid_kmg), title = "Multigrid", label = labels)


println("EGM")
egm =  value_function_egm(c0 = brute_force[3])
plot(egm[1], title = "EGM Value", legend = :topleft, label = labels)
plot(egm[2], title = "EGM Policy", legend = :topleft, label = labels)
plot(egm[3], title = "EGM Consumption", legend = :topleft, label = labels)
plot(EEE2(egm[3], egm[2], egm[5]), title = "EGM", legend = :bottomleft, label = labels)

