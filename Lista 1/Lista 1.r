# Heitor Trielli Zanoni Ferreira - Código - Lista 1

library(dplyr)
library(stats)
library(stargazer)


################
## Questão 1: ##
################
    # Método de Tauchen
    # função que vai trazer a matriz de transição e o grid
    tauchen <- function(grid_len, mu = 0, sigma = 0.007, rho = 0.95, m = 3) {
        if (rho >= 1){stop("Para aplicar Tauchen, a série deve ser estacionária")};
        
        theta_max <- m * sigma / sqrt(1 - rho^2) + mu
        theta_min <- - theta_max + 2 * mu
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
    tauchen_round <- round(tauchen_probs, digits = 3) # Arredondando para ficar mais legível

    stargazer(tauchen_probs)


################
## Questão 2: ##
################
    # Método de 
    # Função que vai computar as probabilidades de transição e o grid. O grid é similar ao de Tauchen
    rouwenhorst <- function(grid_len, mu = 0, sigma = 0.007, rho = 0.95){
        
        # Fazendo o grid
        theta_max <- (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1) + mu
        theta_min <- - theta_max + 2 * mu
        grid <- seq(theta_min, theta_max, length.out = grid_len)


        # Computando a matriz de transição de Rouwenhorst
        p1 <- (1+rho)/2
        p <- matrix(c(p1, 1-p1, 1-p1, p1), nrow = 2)

        if(grid_len >= 3){
            # Cada z abaixo é a matriz que tem zero nas margens e o p anterior no meio. Elas estão dispostas na mesma ordem dos slides
            for (i in 3:grid_len){
                z1 <- cbind(rbind(p, 0), 0) 
                z2 <- rbind(cbind(0, p), 0) 
                z3 <- cbind(rbind(0, p), 0) 
                z4 <- cbind(0, rbind(0, p)) 

                p <- p1*(z1 + z4) + (1-p1)*(z2 + z3) # Computando o valor de p_n
            } #for
        } #if

        transition_matrix <- p/rowSums(p) # Normalizando a matriz para cada linha somar 1

        return(list(probs = transition_matrix, grid = grid))

    } #function 

    rouwen_grid <- rouwenhorst(9)$grid
    rouwen_probs <- rouwenhorst(9)$probs
    rouwen_round <- round(rouwen_probs, digits = 3) # Arredondando para ficar mais legível


################
## Questão 3: ##
################
    # Simulando o AR(1)
    # Cria n valores de um AR(1) utilizando a formula de que y_t = sum_i theta^(t-i) e_i, i = 1, ..., t, assumindo que y_1 = e_1
    ar1 <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 999) {
        
        set.seed(seed) # escolhe o seed
        errors <- rnorm(n, 0, sigma) # gerando vetor de erros
        sample <- c(errors[1]) # O valor inicial da simulação é o primeiro erro

        # Definindo os demais valores simulados recursivamente
        for (i in 2:n){
            sample[i] = mu + rho*sample[i-1] + errors[i]
        } # for
    
        return(list(sample = sample, errors = errors))
    
    } # function 


    # Simulando os métodos discretizados:
    transic <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 999, grid_len = 9, method = 'tauchen', m = 3) {

        erro <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)$errors  # Choques do AR(1)
        cdf_erro <- pnorm(erro, mu, sigma) # Tomando o valor da CDF da normal(0, sigma^2) nos choques
        
        if (method == 'tauchen'){ # Simula para o método de Tauchen

            tauch <- tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            grid <- tauch$grid
            probs <- tauch$probs

        # Transforma a matriz de transição em sua versão CDF
            cdf <- probs
            for (i in 2:grid_len){
                cdf[,i] <- rowSums(probs[,1:i])
            } #for cdf

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]# O valor inicial é o ponto do grid mais perto do primeiro choque

            # Escolhe o theta atual como o primeiro theta em que sua CDF é maior que a CDF do choque
            for (i in 2:n){
                sim[i] <- grid[min(sum(cdf[which(grid == sim[i-1]),] <= cdf_erro[i])+1, grid_len)]
            } # for sim
        } else if (method == 'rouwen'){ # Simula para o método de Rouwenhorst. Aplicado exatamente igual à simulação de Tauchen, apenas mudando a matriz de transição
            
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
                sim[i] <- grid[min(sum(cdf[which(grid == sim[i-1]),] < cdf_erro[i])+1, grid_len)]
            } # for sim
        } # if - else if

        return(sim)
    } # function 

    # Simulando
    ar_sim <- ar1(10000)$sample
    tauch_sim <- transic(10000)
    rouwen_sim <- transic(10000, method = 'rouwen')

    # Plotando para comparar
    plot(ar_sim, type = 'l')
    lines(tauch_sim, col = 'red')

    plot(ar_sim, type = 'l')
    lines(rouwen_sim, col = 'blue')

    # MQE:
    # Tauchen
    mean((ar_sim - tauch_sim)^2)

    # Rouwenhorst
    mean((ar_sim - rouwen_sim)^2)


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
    # Grid e transição
    tauchen_grid_2 <- tauchen(9, rho = 0.99)$grid
    tauchen_probs_2 <- tauchen(9, rho = 0.99)$probs

    # Simulando o AR(1) e o Tauchen
    ar_sim_2 <- ar1(10000, rho = 0.99)$sample
    tauch_sim_2 <- transic(10000, rho = 0.99)

    # Plotando
    plot(ar_sim_2, type = 'l')
    lines(tauch_sim_2, col = 'red')

    # Rouwenhorst
    # Grid e transição
    rouwen_grid_2 <- rouwenhorst(9, rho = 0.99)$grid
    rouwen_probs_2 <- rouwenhorst(9, rho = 0.99)$probs

    # Simulando Rouwenhorst
    rouwen_sim_2 <- transic(10000, rho = 0.99, method = 'rouwen')

    # Plotando
    plot(ar_sim_2, type = 'l')
    lines(rouwen_sim_2, col = 'blue')

    # Regressão de Tauchen
    lm(tauch_sim_2 ~ lag(tauch_sim_2))

    # Regressão de Rouwenhorst
    lm(rouwen_sim_2 ~ lag(rouwen_sim_2) - 1)


    # MQE:
    # Tauchen
    mean((ar_sim - tauch_sim_2)^2)

    # Rouwenhorst
    mean((ar_sim - rouwen_sim_2)^2)
