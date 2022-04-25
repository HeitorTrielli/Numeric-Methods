using Distributions, Random, Plots

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
        vec_N = 1 .- cdf(d,((maximum(grid) .- rho * grid .- delta / 2) / sigma)) # análogo para o maior valor do grid
        grid_interno = grid[2:(length(grid) - 1)] # valores não extremos do grid

        pij = function(j, i = grid) # função que vai computar o vetor de probabilidades de ir para o estado (não extremo) j dado cada estado anterior do grid
            cdf(d,((j .+ delta/2 .- rho * i) / sigma)) - cdf(d,((j .- delta / 2 .- rho * i) / sigma))                             
        end

        mat_interna = reduce(hcat, map(pij, grid_interno)) # map: aplica pij em cada ponto do grid interno; reduce: transforma o vetor de vetores que vem do map em uma matriz
        
        probs = [vec_1 mat_interna vec_N] # gerando a matriz de transição

        return probs, grid
            
    end
    

    probs_tauchen, grid_tauchen = tauchen(9)
    round_tauchen = map(x -> round(x, digits = 3), probs_tauchen)


    teste = [zeros(2) p]

    zero
    
    [p zero]
    
    zeros(2)
    zeros(3)'
    
    [[zeros(2) p]; zeros(3)']
    
    z2 = [[zeros(2) p]
    z4 = [zeros(3)'; [zeros(2) p]]
    z1 = [[p zeros(2)]; zeros(3)']


    # Método de Rouwenhorst

    rouwenhorst = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95)
        theta_max = (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
        theta_min = - theta_max
        grid = LinRange(theta_min, theta_max, grid_len)

        p1 = (1+rho)/2
        p = [p1 (1 - p1) ; (1 - p1) p1]
    
        for i in 3:grid_len
            z1 = [[p zeros(i-1)]; zeros(i)']
            z2 = [[zeros(i-1) p]; zeros(i)']
            z3 = [zeros(i)'; [p zeros(i-1)]]
            z4 = [zeros(i)'; [zeros(i-1) p]]

            
        end
    end


 


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

    return sample
    end




plot(ar1(10000))





