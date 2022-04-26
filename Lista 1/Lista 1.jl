using Distributions, Random, Plots, Pkg


rho = 0.95
sigma = 0.007

# Método de Tauchen:

# Função que vai computar as probabilidades de transição e o grid
    tauchen = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95, m = 3)

        theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
        theta_min = - theta_max # definindo o menor valor do grid
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
            
    end

    tauchen(3)[2]
    

    tauchen_probs, tauchen_grid = tauchen(9)
    tauchen_round = round.(tauchen_probs, digits = 3)

    # Método de Rouwenhorst

    rouwenhorst = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95)
        theta_max = (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
        theta_min = - theta_max
        grid = LinRange(theta_min, theta_max, grid_len)

        p1 = (1+rho)/2
        p = [p1 (1 - p1) ; (1 - p1) p1]
        if grid_len > 2
            for i in 3:grid_len
                z1 = [[p zeros(i-1)]; zeros(i)']
                z2 = [[zeros(i-1) p]; zeros(i)']
                z3 = [zeros(i)'; [p zeros(i-1)]]
                z4 = [zeros(i)'; [zeros(i-1) p]]

                p = p1*(z1 + z4) + (1-p1)*(z2+z3)
            end # for
        end # if
        
        transition_matrix = p./(sum(p, dims = 2))
        
        return transition_matrix, grid
    end # function

    rouwenhorst_probs, rouwenhorst_grid = rouwenhorst(9)
    rouwenhorst_round = round.(rouwenhorst_probs, digits = 3)


# Simulando o AR(1)
    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 = function(n; mu = 0, rho = 0.95, sigma = 0.007, seed = 27) 
        Random.seed!(seed)
        errors = rand(Normal(0, sigma), n) # gerando vetor de erros
        sample = zeros(0)
        append!(sample, errors[1])

        for i in 2:n
            append!(sample, mu + rho*sample[i-1] + errors[i])
        end
        
        return sample, errors
    end


    # Simulando os Markov 
    transic = function(n; mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, type = "tauchen", m = 3)
        erro = ar1(n)[2]

        if type == "tauchen"
            probs, grid = tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)

            sim = zeros(n)
            sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]]
            
            for i in 2:n
                sim[i] = grid[findmin(abs.(grid .- erro[i] .- probs[findall(x -> x == sim[i-1], grid),:]*grid))[2]]
            end #for
        
        elseif type == "rouwen"
            probs, grid = rouwenhorst(grid_len, mu = mu, rho = rho, sigma = sigma)

            sim = zeros(n)
            sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]]
            
            for i in 2:n
                sim[i] = grid[findmin(abs.(grid .- erro[i] .- probs[findall(x -> x == sim[i-1], grid),:]*grid))[2]]
            end # for

        else  
            error("Escolha um dos métodos estudados")
        end # if-elseif

        return sim
    end

plot(transic(10000, type = "tauchen", grid_len = 4))

