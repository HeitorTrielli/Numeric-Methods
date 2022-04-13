rho <- 0.95
sigma <- 0.007

# Método de Tauchen

    # função que vai fazer o grid de n pontos, dados os parâmetros
    gridify <- function(n, sigma = 0.007, rho = 0.95, m = 3) {
        theta_max <- m * sigma / sqrt(1 - rho^2)
        theta_min <- - theta_max

        grid <- seq(theta_min, theta_max, length.out = n) # dividindo [theta_min, theta_max] em n pontos
        return(grid)
    }

    # função que vai trazer a matriz de transição
    tauchen <- function(grid, sigma = 0.007, rho = 0.95, m = 3, seed = 27) {
        delta <- (max(grid) - min(grid)) / (length(grid) - 1) # tamanho dos saltos entre os pontos do grid

        set.seed(seed)
        vec_1 <- pnorm((min(grid) - rho * grid + delta / 2) / sigma) # vetor de transição para o menor valor do grid dado cada estado anterior; pnorm(x) retorna a cdf da normal no valor x
        vec_n <- 1 - pnorm((max(grid) - rho * grid - delta / 2) / sigma) # análogo para o maior valor do grid
        grid_interno <- grid[2:(length(grid) - 1)] # valores não extremos do grid

        # função que retorna o vetor de transição para o estado (não extremo) j dado cada estado anterior
        pij <- function(j, i = grid) {
            pnorm((j + delta / 2 - rho * i) / sigma) - pnorm((j - delta / 2 - rho * i) / sigma)
        }

        mat_interna <- sapply(grid_interno, pij) # aplicando pij em cada valor do grid interno para ter a matriz de transição dos valores não extremos

        return(cbind(vec_1, mat_interna, vec_n)) # combinando as transições extremas com as não extremas
    }

    probs_tauchen <- tauchen(gridify(9))
    round_tauchen <- round(probs_tauchen, digits = 3) # arredondando em três casas decimais



# Simulando o AR(1)

    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 <- function(n, rho = 0.95, sigma = 0.007, seed = 27) {
        if(rho == 1) { stop('rho must be other than 1') }

        set.seed(seed)
        errors <- rnorm(n, 0, sigma^2) # gerando vetor de erros
        rhos <- c(rho^(0:(n - 1))) 
        rhos_inv <- 1 / rhos

        if(rho > 1) {
            coef_mat <- rho %*% t(rhos_inv) * (abs(rhos%*%t(rhos_inv) >= 1)) # matriz que atribui o peso correto (rho^t-i) para cada e_t
        } else {
            coef_mat <- rhos%*%t(rhos_inv) * (abs(rhos%*%t(rhos_inv) <= 1)) # análogo
        }
        sample <- coef_mat%*%errors  


        return(sample)
    }

    plot(ar1(10000, rho = 0.95), type = 'l')


# Simulando o Tauchen

    tauchen_sim <- function(n, rho = 0.95, sigma = 0.007, seed = 27){
        set.seed(seed)
        
    }

