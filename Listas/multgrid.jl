using Plots, BenchmarkTools, Distributions, Distributed, ProfileView, Roots, Dierckx # Pacotes que estou usando
utility = function(c::Float64; mu::Float64 = 2.0)
    return (c^(1 - mu) - 1)/(1 - mu)
    end;
    
    u_line = function(c, mu = 2.0)
        c^(-mu)
    end
    
    u_inv = function(c, mu = 2.0)
        c^(-1/mu)
    end
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
    v0 = utility.(c0);
    

    ############# Utilidade #############
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



############## Multi grid ###############
# Função que faz a interpolação da função valor
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
lastpos = function(val::Array{Float64}, k_possible::Array{Float64}) # Função que vai retornar o ponto maximo da função e o seu respectivo k'
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

############## Accelerator ##############
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
    while error > tol # Checa se convergiu
        @sync @distributed for state in 1:z_len
            for capital in 1:k_len  
                k_possible = k[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0] # Escolhendo os capitais onde o consumo é positivo
                v_possible = value[z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k .> 0, :] # Escolhendo as linhas onde os capitais são positivos
                # val é um vetor com a função valor aplicada em cada possivel escolha de k'
                val = utility.(z[state]*(k[capital]^alpha) + (1 - delta)*k[capital] .- k_possible) + beta*v_possible*p_z[state, :]
                s = argmax(val) # Acha o índice do máximo
                v_next[capital, state] = val[s] # Define o valor como o máximo
                k_line[capital, state] = k_possible[s] # Define o capital como o capital do índice do máximo
            end # for k
        end # for z
        error = maximum(abs.(value - v_next)) # Checa o erro
        value = copy((1 - inert)*v_next + inert*value)
    end # while
    zmat = repeat(z,1,k_len)'
    kmat = repeat(k,1,z_len)   
    c = zmat.*(kmat.^alpha) - k_line + (1-delta)*kmat
    return value::Array{Float64,2}, k_line::Array{Float64,2}, c::Array{Float64,2}
end; # function


value_function_mg = function(g1::Int64, g2::Int64, g3::Int64; z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)


    k = Array(LinRange(0.75*kss, 1.25*kss, g1))
    z = tauchen(z_len)[2]
    zmat = repeat(z,1,g1)'
    kmat = repeat(k,1,z_len)

    v0 = utility.(zmat.*(kmat.^alpha) - delta*kmat)
    
    value1 = value_function_brute(v_0 = v0, k_len = g1, tol = tol)[1]

    v1 = lin_interpol(g1, g2, value1)

    value2 = value_function_brute(v_0 = v1, k_len = g2, tol = tol)[1]

    v2 = lin_interpol(g2, g3, value2)

    value, k_line, c = value_function_brute(v_0 = v2, k_len = g3, tol = tol)

    return value, k_line, c
end # function


println("Multigrid")
multigrid = @benchmark value_function_mg(100, 500, 5000) seconds = 300


value_function_mg_small = function(g1::Int64, g2::Int64; z_len::Int64 = 7, tol::Float64 = 1e-4,
    beta::Float64 = 0.987, mu::Float64 = 2.0, alpha::Float64 = 1 / 3, delta::Float64 = 0.012, rho::Float64 = 0.95, sigma::Float64 = 0.007, inert::Float64 = 0.0,
    kss::Float64 = k_ss)


    k = Array(LinRange(0.75*kss, 1.25*kss, g1))
    z = tauchen(z_len)[2]
    zmat = repeat(z,1,g1)'
    kmat = repeat(k,1,z_len)

    v0 = utility.(zmat.*(kmat.^alpha) - delta*kmat)
    
    value1 = value_function_brute(v_0 = v0, k_len = g1, tol = tol)[1]

    v1 = lin_interpol(g1, g2, value1)

    value, k_line, c = value_function_brute(v_0 = v1, k_len = g2, tol = tol)

    return value, k_line, c
end # function

println("Multigrid Small")
multigrid_small = @benchmark value_function_mg_small(100, 500) seconds = 300

@time value_function_mg(100, 500, 500)