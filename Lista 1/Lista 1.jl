using Distributions, Random, Plots

rho = 0.95
sigma = 0.007
# Método de Tauchen:

    # função que vai criar o grid, dados os parâmetros
    gridify = function (n; sigma = 0.007, rho = 0.95, m = 3)
        theta_max = m * sigma / (sqrt(1 - rho^2)) # definindo o maior valor do grid
        theta_min = - theta_max # definindo o menor valor do grid

        return LinRange(theta_min, theta_max, n) # Cria um vetor de n pontos entre theta_max e theta_min em que a distancia entre os pontos sequencias é igual
    end

    # função que vai computar as probabilidades de transição dado o grid e os parâmetros
    tauchen = function (grid; sigma = 0.007, rho = 0.95, m = 3, seed = 27)
        Random.seed!(seed) 
        d = Normal() # d vira normal(0,1), que será usado para computar a PDF dos erros na hora de achar as probabilidades de transição
        delta = (maximum(grid) - minimum(grid)) / (length(grid) - 1) # distância dos pontos subsequentes do grid

        vec_1 = cdf(d,((minimum(grid) .- rho * grid .+ delta / 2) / sigma)) # vetor das probabilidades de ir para o menor valor do grid, dado cada estado anterior do grid; cdf(d, x) retorna a cdf da distribuição d no valor x
        vec_N = 1 .- cdf(d,((maximum(grid) .- rho * grid .- delta / 2) / sigma)) # análogo para o maior valor do grid
        grid_interno = grid[2:(length(grid) - 1)] # valores não extremos do grid

        pij = function(j, i = grid) # função que vai computar o vetor de probabilidades de ir para o estado (não extremo) j dado cada estado anterior do grid
            cdf(d,((j .+ delta/2 .- rho * i) / sigma)) - cdf(d,((j .- delta / 2 .- rho * i) / sigma))                             
        end

        mat_interna = reduce(hcat, map(pij, grid_interno)) # map: aplica pij em cada ponto do grid interno; reduce: transforma o vetor de vetores que vem do map em uma matriz

        CP = [vec_1 mat_interna vec_N] # combinando os vetores de transição para ter a matriz de transição

        return CP
            
    end

    probs_tauchen = tauchen(gridify(9))
    round_tauchen = map(x -> round(x, digits = 3), probs_tauchen)


# Simulando o AR(1)
    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 = function(n; rho = 0.95, sigma = 0.007, seed = 27) 
        if rho == 1
            error("rho must be other than 1")
        end
        
        Random.seed!(seed)
        errors = rand(Normal(0, sigma^2), n)
        rhos = rho .^ (0:(n - 1))
        rhos_inverse = 1 ./ rhos

        if rho > 1
            coef_mat = rhos*transpose(rhos_inverse) .* (rhos*transpose(rhos_inverse) .>= 1)      
        else 
            coef_mat = rhos*transpose(rhos_inverse) .* (rhos*transpose(rhos_inverse) .<= 1)
        end

        sample = coef_mat * errors

    return sample
    end

plot(ar1(10000))