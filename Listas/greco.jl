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
gdp_pp = function(;N_k=N_k,N_z=N_z,k_grid=k_grid,z_grid=Tauch[1]) #Criando uma função para calcular o grid de y (para ele não ficar dentro de loops)
    y=zeros(N_k,N_z)
    for state in 1:N_z
        ptf=exp(z_grid[state])
        for capital in 1:N_k
            y[capital,state] = ptf*(k_grid[capital])^alfa+(1-delta)*k_grid[capital]
        end #End for de z
    end #End for de k

    return y
end #End function

Tauch = tauchen(7)
# Definindo as variáveis
beta, mu, alpha, delta, rho, sigma = 0.987, 2.0, 1/3, 0.012, 0.95, 0.007;

# Definindo o capital de steady state

k_ss = (alpha / (1 / beta - 1 + delta))^(1 / (1 - alpha));
z_len = 7;
k_len = 500;
grid_z = exp.(tauchen(z_len)[2]); # Valores que exp(z_t) pode tomar 
k_grid = Array(LinRange(0.75*k_ss, 1.25*k_ss, k_len));
value_0 = zeros(500,7)
policy_0 = zeros(500,7)


value_function_iteration=function(;Tauch=Tauch, k_grid=k_grid, tol = .0001, beta = .987, mu = 2, alfa = 1/3, delta = 0.012, N_k = 500, N_z = 7, value_0=value_0,policy_0=policy_0)
    
    policy=copy(policy_0)
    value=copy(value_0)
    next_capital=copy(value_0)
    consumption=copy(value_0)
    EEE=copy(value_0)
    pmk=zeros(N_k,N_z)

    y=gdp_pp(N_k=N_k,N_z=N_z,z_grid=Tauch[1],k_grid=k_grid) 

    for i_iterat in 1:10000
        Threads.@threads for state in 1:N_z
            prob = Tauch[2][state,:]             
            # Preenchendo a primeira coluna das funções
            capital = 1
            candidato_1 = u(y[capital,state]-k_grid[1])+beta*prob'value_0[1,:]
            max_capital=count(x->x<=y[capital,state],k_grid)
            for t in 2:max_capital
                candidato_2 = u(y[capital,state]-k_grid[t])+beta*prob'value_0[t,:]
                # Se a função valor decrescer, passamos do ponto ótimo, logo
                if candidato_2-candidato_1<0
                    value[capital,state] = candidato_1
                    policy[capital,state] = t-1
                    break
                elseif t==max_capital
                    value[capital,state] = candidato_2
                    policy[capital,state] = max_capital
                end # End if 
                candidato_1=candidato_2
            end # End for (testando os valores possiveis para a função politica)
            
            # Para as próximas colunas
            for capital in 2:N_k
                # como a função politica é crescente em k, só precisamos testar os valores a partir da função politica do k anterior
                i_0 = policy[(capital)-1,state]
                candidato_1 = u(y[capital,state]-k_grid[i_0])+beta*(prob'value_0[i_0,:])
                max_capital=count(x->x<=y[capital,state],k_grid)
                if i_0<N_k
                    for t in (i_0+1):max_capital
                        candidato_2 = u(y[capital,state]-k_grid[t])+beta*prob'value_0[t,:]
                        # Se a função valor decrescer, passamos do ponto ótimo, logo
                        if candidato_2-candidato_1<0
                            value[capital,state] = candidato_1
                            policy[capital,state] = t-1
                            break
                        elseif t==max_capital
                            value[capital,state] = candidato_2
                            policy[capital,state] = max_capital
                        end # End if 
                        candidato_1=candidato_2
                    end # End for (testando os valores possiveis para a função politica)            
                else
                    value[capital,state] = candidato_1
                    policy[capital,state] = N_k
                end #End if
            end # End for de k
        end #End for de z
        if maximum(abs.(value-value_0))<tol
            break
        else
            value_0=copy(value)
            policy_0=copy(policy)
        end #End condição parada  
          
    end #End iteração da função valor

    Threads.@threads for state in 1:N_z #Preenchendo a função k_linha
        next_capital[:,state]=k_grid[policy[:,state]]
        pmk[:,state]=alfa.*exp(Tauch[1][state]).*next_capital[:,state].^(alfa-1).+(1-delta)
    end #End for

    consumption=y-next_capital

    Threads.@threads for state in 1:N_z # Preenchendo a matriz dos erros de equação de Euler
        prob = Tauch[2][state,:]
        for capital in 1:N_k
            c_line=consumption[policy[capital,state],:]
            E = prob'u_marginal.(c_line)
            EEE[capital,state]=log10(abs(1-(u_marginal_inverse(beta*E)/consumption[capital,state])))
        end #End for de k
    end #End for de z

    return value, policy, consumption, next_capital, EEE

end #End function

@time value_function_iteration()



value_function_iteration_brute = function(;Tauch=Tauch, k_grid=k_grid, tol = .0001, beta = .987, mu = 2, alfa = 1/3, delta = 0.012, N_k = 500, N_z = 7, value_0=value_0,policy_0=policy_0)
    
    policy=copy(policy_0)
    value=copy(value_0)
    next_capital=copy(value_0)
    consumption=copy(value_0)
    EEE=copy(value_0)
    pmk=zeros(N_k,N_z)

    y=gdp_pp(N_k=N_k,N_z=N_z,z_grid=Tauch[1],k_grid=k_grid)

    for i_iterat in 1:10000
        Threads.@threads for state in 1:N_z
            prob = Tauch[2][state,:]                      
            # Para todas as colunas
            for capital in 1:N_k
                # como a função politica é crescente em k, só precisamos testar os valores a partir da função politica do k anterior
                i_0 = 1
                candidato_1 = u(y[capital,state]-k_grid[i_0])+beta*(prob'value_0[i_0,:])
                max_capital=count(x->x<=y[capital,state],k_grid)
                if i_0<N_k
                    for t in 2:max_capital
                        candidato_2 = u(y[capital,state]-k_grid[t])+beta*prob'value_0[t,:]
                        # Se a função valor decrescer, passamos do ponto ótimo, logo
                        if candidato_2-candidato_1<0
                            value[capital,state] = candidato_1
                            policy[capital,state] = t-1
                            break
                        elseif t==max_capital
                            value[capital,state] = candidato_2
                            policy[capital,state] = max_capital
                        end # End if 
                        candidato_1=candidato_2
                    end # End for (testando os valores possiveis para a função politica)            
                else
                    value[capital,state] = candidato_1
                    policy[capital,state] = N_k
                end #End if
            end # End for de k
        end #End for de z
        if maximum(abs.(value-value_0))<tol
            break
        else
            value_0=copy(value)
            policy_0=copy(policy)
        end #End condição parada  
          
    end #End iteração da função valor

    Threads.@threads for state in 1:N_z #Preenchendo a função k_linha
        next_capital[:,state]=k_grid[policy[:,state]]
        pmk[:,state]=alfa.*exp(Tauch[1][state]).*next_capital[:,state].^(alfa-1).+(1-delta)
    end #End for

    consumption=y-next_capital

    Threads.@threads for state in 1:N_z # Preenchendo a matriz dos erros de equação de Euler
        prob = Tauch[2][state,:]
        for capital in 1:N_k
            c_line=consumption[policy[capital,state],:]
            E = prob'u_marginal.(c_line)
            EEE[capital,state]=log10(abs(1-(u_marginal_inverse(beta*E)/consumption[capital,state])))
        end #End for de k
    end #End for de z

    return value, policy, consumption, next_capital, EEE

end #End function