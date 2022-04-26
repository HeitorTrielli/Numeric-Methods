library(dplyr)
library(stats)

rho <- 0.95
sigma <- 0.007
mu <- 0
################
## Questão 1: ##
################
# Método de Tauchen
    # função que vai trazer a matriz de transição e o grid
    tauchen <- function(grid_len, mu = 0, sigma = 0.007, rho = 0.95, m = 3) {
        if (rho >= 1){stop("Para aplicar Tauchen, a série deve ser estacionária")}
        
        theta_max <- m * sigma / sqrt(1 - rho^2) + mu
        theta_min <- - theta_max
        grid <- seq(theta_min, theta_max, length.out = grid_len) # dividindo [theta_min, theta_max] em n pontos

        delta <- (max(grid) - min(grid)) / (length(grid) - 1) # tamanho dos saltos entre os pontos do grid

        vec_1 <- pnorm((min(grid) - rho * grid + delta / 2) / sigma, mu) # vetor de transição para o menor valor do grid dado cada estado anterior; pnorm(x) retorna a cdf da normal no valor x
        vec_n <- 1 - pnorm((max(grid) - rho * grid - delta / 2) / sigma, mu) # análogo para o maior valor do grid
        
        if(grid_len > 2){

            grid_interno <- grid[2:(length(grid) - 1)] # valores não extremos do grid

            # função que retorna o vetor de transição para o estado (não extremo) j dado cada estado anterior
            pij <- function(j, i = grid) {
                pnorm((j + delta / 2 - rho * i) / sigma, mu) - pnorm((j - delta / 2 - rho * i) / sigma, mu)
            }

            mat_interna <- sapply(grid_interno, pij)

            return(list(probs = cbind(vec_1, mat_interna, vec_n), grid = grid)) # combinando as transições extremas com as não extremas 
        } else {
            return(list(probs = cbind(vec_1, vec_n), grid = grid))
        } # else

    } # function

    tauchen_grid <- tauchen(9)$grid
    tauchen_probs <- tauchen(9)$probs
    tauchen_round <- round(tauchen_probs, digits = 3) # arredondando em três casas decimais

################
## Questão 2: ##
################
# Método de Rouwenhorst
    rouwenhorst <- function(grid_len, mu = 0, sigma = 0.007, rho = 0.95){
        
        # Fazendo o grid
        theta_max <- (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
        theta_min <- - theta_max
        grid <- seq(theta_min, theta_max, length.out = grid_len)


        # Fazendo a matriz de transição
        p1 <- (1+rho)/2
        p <- matrix(c(p1, 1-p1, 1-p1, p1), nrow = 2)

        if(grid_len >= 3){

            for (i in 3:grid_len){
                z1 <- cbind(rbind(p, 0), 0) # adciona um vetor de zeros em baixo de p, depois adciona a direita
                z2 <- rbind(cbind(0, p), 0) # análogo
                z3 <- cbind(rbind(0, p), 0) # análogo
                z4 <- cbind(0, rbind(0, p)) # análogo

                p <- p1*(z1 + z4) + (1-p1)*(z2 + z3)
            } #for
        } #if
        return(list(probs = p/rowSums(p), grid = grid))

    } #function 

    rouwen_grid <- rouwenhorst(9)$grid
    rouwen_probs <- rouwenhorst(9)$probs
    rouwen_round <- round(rouwen_probs, digits = 3)


################
## Questão 3: ##
################
# Simulando o AR(1)

    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27) {
        
        set.seed(seed)
        errors <- rnorm(n, 0, sigma) # gerando vetor de erros
        sample <- c(errors[1])

        for (i in 2:n){
            sample[i] = mu + rho*sample[i-1] + errors[i]
        } # for
        return(list(sample = sample, errors = errors))
    } # function 


# Simulando os Markov:
 transic <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, type = 'tauchen', m = 3) {

        erro <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)$errors
        cdf_erro <- pnorm(erro, mu, sigma)
        
        if (type == 'tauchen'){

            tauch <- tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            grid <- tauch$grid
            probs <- tauch$probs

            cdf <- probs
            for (i in 2:grid_len){
                cdf[,i] <- rowSums(probs[,1:i])
            } #for cdf

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]

            for (i in 2:n){
                sim[i] <- grid[min(sum(cdf[which(grid == sim[i-1]),] <= cdf_erro[i])+1, grid_len)]
            } # for sim
        } else if (type == 'rouwen'){
            
            rouwen <- rouwenhorst(grid_len, mu = mu, rho = rho, sigma = sigma)
            grid <- rouwen$grid
            probs <- rouwen$probs

            cdf <- probs
            for (i in 2:grid_len){
                cdf[,i] <- rowSums(probs[,1:i])
            } #for cdf

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]

            for (i in 2:n){
                sim[i] <- grid[min(sum(cdf[which(grid == sim[i-1]),] <= cdf_erro[i])+1, grid_len)]
            } # for sim
        }

    return(sim)
 }
# plotando para comparar

ar_sim <- ar1(10000)$sample
tauch_sim <- transic(10000)
rouwen_sim <- transic(10000, type = 'rouwen')

plot(ar_sim, type = 'l')
lines(tauch_sim, col = 'red')

plot(ar_sim, type = 'l')
lines(rouwen_sim, col = 'blue')


################
## Questão 4: ##
################

# Regressão do Tauchen:

lm(tauch_sim ~ lag(tauch_sim) - 1)


# Regressão do Rouwenhorst

lm(rouwen_sim ~ lag(rouwen_sim) - 1)




################
## Questão 5: ##
################

# Tauchen:

tauchen_grid_2 <- tauchen(9, rho = 0.99)$grid
tauchen_probs_2 <- tauchen(9, rho = 0.99)$probs

ar_sim_2 <- ar1(10000, rho = 0.99)$sample
tauch_sim_2 <- transic(10000, rho = 0.99)

lm(tauch_sim_2 ~ lag(tauch_sim_2))


plot(ar_sim_2, type = 'l')
lines(tauch_sim_2, col = 'red')

# Rouwenhorst
rouwen_grid_2 <- rouwenhorst(9, rho = 0.99)$grid
rouwen_probs_2 <- rouwenhorst(9, rho = 0.99)$probs

rouwen_sim_2 <- transic(10000, rho = 0.99, type = 'rouwen')

lm(rouwen_sim_2 ~ lag(rouwen_sim_2) - 1)

plot(ar_sim_2, type = 'l')
lines(transic(10000, rho = 0.99, type = 'rouwen'), col = 'blue')




