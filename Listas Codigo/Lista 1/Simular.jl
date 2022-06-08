module Simular

include("Tauchen.jl")
include("Rouwenhorst.jl")

using Random, Distributions, Main.Tauchen, Main.Rouwenhorst

export ar1, transic

# Simulando o AR(1)
    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 = function(n; mu = 0, rho = 0.95, sigma = 0.007, seed = 999) 
       
        Random.seed!(seed) # escolhe o seed
        errors = rand(Normal(0, sigma), n) # gerando vetor de erros
        sample = zeros(n)
        sample[1] = errors[1] # O valor inicial da simulação é o primeiro erro

        # Definindo os demais valores simulados recursivamente
        for i in 2:n
            sample[i] = mu + rho*sample[i-1] + errors[i] 
        end # for
        
        return sample, errors

    end; # function


    # Simulando os métodos discretizados 
    transic = function(n; mu = 0, rho = 0.95, sigma = 0.007, seed = 999, grid_len = 9, method = "tauchen", m = 3)
        
        erro = ar1(n)[2] # Choques do AR(1)
        cdf_erro = cdf.(Normal(0, sigma), erro) # Tomando o valor da CDF da normal(0, sigma^2) nos choques
    
        if method == "tauchen" # Simula para o método de Tauchen
            probs, grid = tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            
            # Transforma a matriz de transição em sua versão CDF
            CDF = zeros(grid_len, grid_len)
            for i in 1:grid_len
                CDF[:,i] = sum(probs[:,1:i], dims = 2)
            end # for cdf
    
            sim = zeros(n)
            sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]] # O valor inicial é o ponto do grid mais perto do primeiro choque
            
            # Escolhe o theta atual como o primeiro theta em que sua CDF é maior que a CDF do choque
            for i in 2:n
                sim[i] = grid[minimum([sum(CDF[findall(x-> x == sim[i-1], grid),:] .<= cdf_erro[i])+1, grid_len])]
            end #for sim
        
        elseif method == "rouwen" # Simula para o método de Rouwenhorst. Aplicado exatamente igual à simulação de Tauchen, apenas mudando a matriz de transição
            probs, grid = rouwenhorst(grid_len, mu = mu, rho = rho, sigma = sigma)
    
            CDF = zeros(grid_len, grid_len)
            for i in 1:grid_len
                CDF[:,i] = sum(probs[:,1:i], dims = 2)
            end #for
    
            sim = zeros(n)
            sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]]
            
            for i in 2:n
                sim[i] = grid[minimum([sum(CDF[findall(x-> x == sim[i-1], grid),:] .< cdf_erro[i])+1, grid_len])]
            end # for
    
        else  
            error("Escolha um dos métodos estudados")
        end # if - elseif - else
    
        return sim
    end; # function
end

