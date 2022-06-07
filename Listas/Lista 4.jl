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
#

beta, gamma, rho, sigma = 0.96,  1.0001, 0.96, 0.001
z_len, a_len = 9, 100
grid_z = exp.(tauchen(z_len; rho = rho, sigma = sigma).grid)
p_z = tauchen(z_len, rho = rho, sigma = sigma).prob

zmat = repeat(z,1,k_len)';
kmat = repeat(k,1,z_len);
utility = function(c::Float64; mu::Float64 = gamma)
    return (c^(1 - mu) - 1)/(1 - mu)
end;
u_line = function(c, mu = 2.0)
    c^(-mu)
end;
u_inv = function(c, mu = 2.0)
    c^(-1/mu)
end;