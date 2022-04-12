# Método de Tauchen

# função que vai fazer o grid de n pontos, dados os parâmetros
gridify <- function(n, sigma = 0.007, rho = 0.95, m = 3){
    theta_max <- m*sigma/sqrt(1-rho^2)
    theta_min <- -theta_max

    grid <- seq(theta_min, theta_max, length.out = n) # dividindo [theta_min, theta_max] em n pontos
    return(grid)

}

# função que vai trazer a matriz de transição 
tauchen <- function(grid, sigma = 0.007, rho = 0.95, m = 3){
    
    delta = (max(grid)-min(grid))/(length(grid)-1) # tamanho dos saltos entre os pontos do grid

    
    vec_1 <- pnorm((min(grid)-rho*grid+delta/2)/sigma) # vetor de transição para o menor valor do grid dado cada estado anterior; pnorm(x) retorna a cdf da normal no valor x
    vec_N <- 1- pnorm((max(grid)-rho*grid-delta/2)/sigma) # análogo para o maior valor do grid
    grid_interno <- grid[2:(length(grid)-1)] # valores não extremos do grid

    # função que retorna o vetor de transição para o estado (não extremo) j dado cada estado anterior
    pij <- function(j, i = grid){ 
        pnorm((j + delta/2-rho*i)/sigma) - pnorm((j - delta/2-rho*i)/sigma)  
    }

    mat_interna <- sapply(grid_interno, pij) # aplicando pij em cada valor do grid interno para ter a matriz de transição dos valores não extremos

    return(cbind(vec_1, mat_interna, vec_N)) # combinando as transições extremas com as não extremas
}

probs_tauchen <- tauchen(gridify(9))

