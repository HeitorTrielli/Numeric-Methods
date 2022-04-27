using Distributions, Random, Plots, GLM, Pkg, DataFrames, PrettyTables, RegressionTables # Pacotes que eu vou usar

################
## Questão 1: ##
################
    # Método de Tauchen:
    # Função que vai computar as probabilidades de transição e o grid
    tauchen = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95, m = 3)

        theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
        theta_min = - theta_max + 2 * mu # definindo o menor valor do grid
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
            
    end;

    tauchen_probs, tauchen_grid = tauchen(9);
    tauchen_round = round.(tauchen_probs, digits = 3); # Arredondando para ficar mais legível

################
## Questão 2: ##
################
    # Método de Rouwenhorst
    # Função que vai computar as probabilidades de transição e o grid.
    rouwenhorst = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95)
        # Fazendo o grid
        theta_max = (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1) + mu
        theta_min = - theta_max + 2 * mu
        grid = LinRange(theta_min, theta_max, grid_len)

        # Computando a matriz de transição de Rouwenhorst
        p1 = (1+rho)/2 # p inicial
        p = [p1 (1 - p1) ; (1 - p1) p1] # matriz inicial 2x2
        if grid_len > 2
            # Cada z abaixo é a matriz que tem zero nas margens e o p anterior no meio. Elas estão dispostas na mesma ordem dos slides
            for i in 3:grid_len
                z1 = [[p zeros(i - 1)]; zeros(i)'] 
                z2 = [[zeros(i - 1) p]; zeros(i)']
                z3 = [zeros(i)'; [p zeros(i - 1)]]
                z4 = [zeros(i)'; [zeros(i - 1) p]]

                p = p1 * (z1 + z4) + (1-p1) * (z2 + z3) # Computando o valor de p_n
            end # for
        end # if
        
        transition_matrix = p./(sum(p, dims = 2)) # Normalizando a matriz para cada linha somar 1
        
        return transition_matrix, grid
    end; # function

    rouwenhorst_probs, rouwenhorst_grid = rouwenhorst(9);
    rouwenhorst_round = round.(rouwenhorst_probs, digits = 3); # Arredondando para ficar mais legível


################
## Questão 3: ##
################
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

    # plotando para comparar
    ar_sim = ar1(10000)[1];
    tauch_sim = transic(10000);
    rouwen_sim = transic(10000, method = "rouwen");

    plot([ar_sim, tauch_sim], label = ["AR(1)" "Tauchen"])
    plot([ar_sim, rouwen_sim], label = ["AR(1)" "Rouwenhorst"])

    # MQE
    pretty_table([mean((ar_sim - tauch_sim).^2) mean((ar_sim - rouwen_sim).^2)], header = ["Tauchen", "Rouwenhorst"])



################
## Questão 4: ##
################
    # Fazendo o dataframe para as regressões:
    df = DataFrame(Tauchen = tauch_sim, LagTauchen = lag(tauch_sim), Rouwenhorst = rouwen_sim, LagRouwen = lag(rouwen_sim));

    # Regressão do Tauchen
    tauchen_reg = lm(@formula(Tauchen ~ 0 + LagTauchen), df);

    # Regressão do Rouwenhorst
    rouwen_reg = lm(@formula(Rouwenhorst ~ 0 + LagRouwen), df);

    regtable(tauchen_reg, rouwen_reg)
    
################
## Questão 5: ##
################
    # Tauchen:
    # Grid e transição
    tauchen_grid_2 = tauchen(9, rho = 0.99)[2];
    tauchen_probs_2 = tauchen(9, rho = 0.99)[1];

    # Simulações
    ar_sim_2 = ar1(10000, rho = 0.99)[1];
    tauch_sim_2 = transic(10000, rho = 0.99);

    # Plotando
    plot([ar_sim_2 tauch_sim_2], label = ["AR(1)" "Tauchen"])

    # Rouwenhorst
    rouwen_grid_2 = rouwenhorst(9, rho = 0.99)[2];
    rouwen_probs_2 = rouwenhorst(9, rho = 0.99)[1];

    rouwen_sim_2 = transic(10000, rho = 0.99, method = "rouwen");

    plot([ar_sim_2 rouwen_sim_2], label = ["AR(1)" "Rouwenhorst"])

    pretty_table([mean((ar_sim - tauch_sim_2).^2) mean((ar_sim - rouwen_sim_2).^2)], header = ["Tauchen", "Rouwenhorst"])

    #Regressões
    df_2 = DataFrame(Tauchen = tauch_sim_2, LagTauchen = lag(tauch_sim_2), Rouwenhorst = rouwen_sim_2, LagRouwen = lag(rouwen_sim_2));

    # Regressão do Tauchen
    tauchen_reg_2 = lm(@formula(Tauchen ~ 0 + LagTauchen), df_2);

    # Regressão do Rouwenhorst
    rouwen_reg_2 = lm(@formula(Rouwenhorst ~ 0 + LagRouwen), df_2);

    regtable(tauchen_reg_2, rouwen_reg_2)