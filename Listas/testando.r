library(dplyr)
library(stats)

rho <- 0.95
sigma <- 0.007
mu <- 0


# Método de Tauchen
    # função que vai trazer a matriz de transição
    tauchen <- function(grid_len, mu = 0, sigma = 0.007, rho = 0.95, m = 3) {
        if (rho >= 1){stop("Para aplicar Tauchen, a série deve ser estacionária")}
        
        theta_max <- m * sigma / sqrt(1 - rho^2) + mu
        theta_min <- - theta_max
        grid <- seq(theta_min, theta_max, length.out = grid_len) # dividindo [theta_min, theta_max] em n pontos

        delta <- (max(grid) - min(grid)) / (length(grid) - 1) # tamanho dos saltos entre os pontos do grid

        vec_1 <- pnorm((min(grid) - rho * grid + delta / 2) / sigma, mu) # vetor de transição para o menor valor do grid dado cada estado anterior; pnorm(x) retorna a cdf da normal no valor x
        vec_n <- 1 - pnorm((max(grid) - rho * grid - delta / 2) / sigma, mu) # análogo para o maior valor do grid
        grid_interno <- grid[2:(length(grid) - 1)] # valores não extremos do grid

        # função que retorna o vetor de transição para o estado (não extremo) j dado cada estado anterior
        pij <- function(j, i = grid) {
            pnorm((j + delta / 2 - rho * i) / sigma, mu) - pnorm((j - delta / 2 - rho * i) / sigma, mu)
        }

        mat_interna <- sapply(grid_interno, pij)

        return(list(probs = cbind(vec_1, mat_interna, vec_n), grid = grid)) # combinando as transições extremas com as não extremas
    }

    tauchen_grid <- tauchen(9)$grid
    tauchen_pdf <- tauchen(9)$probs
    tauchen_round <- round(tauchen_pdf, digits = 3) # arredondando em três casas decimais

    

# Método de Rouwenhorst
    rouwenhorst <- function(n, mu = 0, sigma = 0.007, rho = 0.95){
        
        # Fazendo o grid
        theta_max <- (sigma/sqrt(1-rho^2)) * sqrt(n-1)
        theta_min <- - theta_max
        grid <- seq(theta_min, theta_max, length.out = n)


        # Fazendo a matriz de transição
        p1 <- (1+rho)/2
        p <- matrix(c(p1, 1-p1, 1-p1, p1), nrow = 2)

        for (i in 3:n){
            z1 <- cbind(rbind(p, rep(0,(i - 1))), rep(0, i))
            z2 <- rbind(cbind(rep(0, (i - 1)), p), rep(0,i))
            z3 <- cbind(rbind(rep(0, (i - 1)), p), rep(0,i))
            z4 <- cbind(rep(0, i), rbind(rep(0,(i - 1)), p))


            p <- p1*(z1 + z4) + (1-p1)*(z2 + z3)
        }
        return(list(probs = p/rowSums(p), grid = grid))
    }

    rouwen_grid <- rouwenhorst(9)$grid
    rouwen_pdf <- rouwenhorst(9)$probs
    rouwen_round <- round(rouwen_pdf, digits = 3)


# Simulando o AR(1)

    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27) {
        
        set.seed(seed)
        errors <- rnorm(n, 0, sigma) # gerando vetor de erros
        sample <- c(errors[1])

        for (i in 2:n){
            sample[i] = mu + rho*sample[i-1] + errors[i]
        }
        return(list(sample = sample, errors = errors))
    }



plot(ar1(10000)[[1]], type = 'l')

# Simulando os Markov:
    
    transic <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, type = 'tauchen', m = 3) {

        erro <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)$errors

        if(type == 'tauchen'){
            tauch <- tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            grid <- tauch$grid
            probs <- tauch$probs

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]

            # A simulação escolhe theta_i como o valor do grid mais perto de E(theta_i|theta_{i-1}) + erro[i]
            for (i in 2:n){
                sim[i] <- grid[which(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))==min(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))))]
            }
            
            return(sim)

        } else if (type == 'rouwen') { # Análogo à simulação de Tauchen

            rouwen <- rouwenhorst(n, mu = mu, rho = rho, sigma = sigma)
            grid <- rouwen$grid
            probs <- rouwen$probs

            sim <- c()
            sim[1] <- grid[which(abs(grid - erro[1]) == min(abs(grid - erro[1])))]

            for (i in 2:n){
                sim[i] <- grid[which(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),])) == min(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))))]
            }
            return(sim)

        } else {stop('Escolha um dos métodos estudados')}
    }

transic(10000)


plot(ar1(10000)$sample, type = 'l')
lines(transic(10000, grid_len = 4), col = 'red')

n = 10000

rouwen <- rouwenhorst(n, mu = mu, rho = rho, sigma = sigma)
grid <- rouwen$grid
probs <- rouwen$probs

sim <- c()
sim[1] <- grid[which(abs(grid - erro[1]) == min(abs(grid - erro[1])))]

for (i in 2:n){
    sim[i] <- grid[which(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),])) == min(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))))]
}
return(sim)







