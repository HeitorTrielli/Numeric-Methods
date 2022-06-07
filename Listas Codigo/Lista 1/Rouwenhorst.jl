module Rouwenhorst

export rouwenhorst

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
    
    return (prob = transition_matrix, grid = grid)
end; # function

end