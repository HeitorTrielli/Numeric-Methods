library(stats)

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

teste <- rnorm(10000, 0, 0.007)

tauchen_grid <- function(n, mu = 0, sigma = 0.007, rho = 0.95, m = 3) {
        if (rho >= 1){stop("Para aplicar Tauchen, a série deve ser estacionária")}
        theta_max <- m * sigma / sqrt(1 - rho^2) + mu
        theta_min <- - theta_max

        grid <- seq(theta_min, theta_max, length.out = n) # dividindo [theta_min, theta_max] em n pontos
        return(grid)
    }

    # função que vai trazer a matriz de transição
    tauchen <- function(n, mu = 0, sigma = 0.007, rho = 0.95, m = 3, seed = 27, grid_len = 9) {
        grid <- tauchen_grid(grid_len)
        delta <- (max(grid) - min(grid)) / (length(grid) - 1) # tamanho dos saltos entre os pontos do grid

        set.seed(seed)
        vec_1 <- pnorm((min(grid) - rho * grid + delta / 2) / sigma, mu) # vetor de transição para o menor valor do grid dado cada estado anterior; pnorm(x) retorna a cdf da normal no valor x
        vec_n <- 1 - pnorm((max(grid) - rho * grid - delta / 2) / sigma, mu) # análogo para o maior valor do grid
        grid_interno <- grid[2:(length(grid) - 1)] # valores não extremos do grid

        # função que retorna o vetor de transição para o estado (não extremo) j dado cada estado anterior
        pij <- function(j, i = grid) {
            pnorm((j + delta / 2 - rho * i) / sigma, mu) - pnorm((j - delta / 2 - rho * i) / sigma, mu)
        }

        mat_interna <- sapply(grid_interno, pij)

        return(list(prob = cbind(vec_1, mat_interna, vec_n), grid = grid, delta = delta)) # combinando as transições extremas com as não extremas
    }












transic <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, tauchen = TRUE){       
    if(grid_len <= 2){
        stop('Escolha ao menos três pontos')
    }

    if (tauchen == TRUE){
    grid <- tauchen_grid(grid_len, mu = mu, rho = rho, sigma = sigma)
    grid_interno <- grid[2:(length(grid)-1)]
    delta <- (max(grid) - min(grid))/(length(grid)-1)
    ar <- ar1(10000, mu = mu, rho = rho, sigma = sigma, seed = seed)
    ar_sim <- ar[[1]]
    errors <- ar[[2]]


    v1 <- min(grid)*(rho*ar_sim + errors <= min(grid) + delta/2)
    vn <- max(grid)*(rho*ar_sim + errors >= max(grid) - delta/2)
    v_int <- sapply(grid_interno, function(x){x*(x <= rho*ar_sim + delta/2 + errors)*( x >= rho*ar_sim - delta/2 + errors)})  
    return(v1+rowSums(v_int)+vn)
    }
}





grid <- tauchen_grid(grid_len, mu = 0, rho = 0.95, sigma = 0.007)

transic2 <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, tauchen = TRUE){       
    if(grid_len <= 2){
        stop('Escolha ao menos três pontos')
    }

    if (tauchen == TRUE){
    grid <- tauchen_grid(grid_len, mu = mu, rho = rho, sigma = sigma)
    grid_interno <- grid[2:(length(grid)-1)]
    delta <- (max(grid) - min(grid))/(length(grid)-1)
    ar <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)
    ar_sim <- ar[[1]]

    v1 <- min(grid)*(ar_sim <= min(grid)+delta/2)
    vn <- max(grid)*(ar_sim >= max(grid)-delta/2)
    v_int <- sapply(grid_interno, function(x){x*(ar_sim > x-delta/2)*(ar_sim < x+delta/2)})

    return( v1 + rowSums(v_int) + vn )
    }
}



grid <- tauchen_grid(grid_len, mu = mu, rho = rho, sigma = sigma)
grid_interno <- grid[2:(length(grid)-1)]
delta <- (max(grid) - min(grid))/(length(grid)-1)
ar <- ar1(10000, mu = mu, rho = rho, sigma = sigma, seed = seed)
ar_sim <- ar[[1]]
errors <- ar[[2]]


v1 <- min(grid)*(rho*ar_sim + errors <= min(grid) + delta/2)

vn <- max(grid)*(rho*ar_sim + errors >= max(grid) - delta/2)

v_int <- sapply(grid_interno, function(x){x*(x <= rho*ar_sim + delta/2 + errors)*( x > rho*ar_sim - delta/2 + errors)})


unique(v1+rowSums(v_int)+vn)


n <- 10000 
probs_tauchen <- tauchen(grid_len)[[1]]
tauchen_cdf <- matrix(0, grid_len, grid_len)
errors <- ar1(10000)$errors

tauchen_cdf[,1] <- probs_tauchen[,1]
for (i in 2:grid_len){
tauchen_cdf[,i] <- rowSums(probs_tauchen[,1:i])
}

sample <- c(grid[which(abs(grid-errors[1]) == min(abs(grid-errors[1])))])

for (i in 2:n){
    sample[i] <- grid[which(abs(tauchen_cdf[which(grid == sample[i-1]),]-pnorm(errors[i],0, sigma)) == min(abs(tauchen_cdf[which(grid == sample[i-1]),]-pnorm(errors[i],0, sigma))))]
}


            tauchen <- tauchen(9)
            grid <- tauchen[[2]]
            probs <- tauchen[[1]]
            
            erro <- ar1(10000)[[2]]

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]


            for (i in 2:10000){
                sim[i] <- grid[which(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))==min(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))))]
            }

sim

ar <- ar1(10000)[[1]]


plot(ar, type = 'l')

lines(sim, col = 'red')




            erro <- ar1(10000)$errors

            tauch <- tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            grid <- tauch$grid
            probs <- tauch$probs

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]

            # A simulação escolhe theta_i como o valor do grid mais perto de E(theta_i|theta_{i-1}) + erro[i]
            for (i in 2:n){
                sim[i] <- grid[which(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))==min(abs(grid - erro[i] - weighted.mean(grid, probs[which(grid == sim[i-1]),]))))]
                print(sim[i])
            }
            




grid_len <- 5
        
        # Fazendo o grid
        theta_max <- (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
        theta_min <- - theta_max
        grid <- seq(theta_min, theta_max, length.out = grid_len)


        # Fazendo a matriz de transição
        p1 <- (1+rho)/2
        p <- matrix(c(p1, 1-p1, 1-p1, p1), nrow = 2)

        if(grid_len >= 3){

            for (i in 3:grid_len){
                z1 <- cbind(rbind(p, rep(0,(i - 1))), rep(0, i))
                z2 <- rbind(cbind(rep(0, (i - 1)), p), rep(0,i))
                z3 <- cbind(rbind(rep(0, (i - 1)), p), rep(0,i))
                z4 <- cbind(rep(0, i), rbind(rep(0,(i - 1)), p))


                p <- p1*(z1 + z4) + (1-p1)*(z2 + z3)
            }
        }
        




# Simulando o por transição


    transic <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, tauchen = TRUE){       
        
        if (tauchen == TRUE){
        grid <- tauchen_grid(grid_len, mu = mu, rho = rho, sigma = sigma)
        grid_interno <- grid[2:(length(grid)-1)]
        delta <- (max(grid) - min(grid))/(length(grid)-1)
        ar <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)
        ar_sim <- ar$sample
        errors <- ar$errors

        v1 <- min(grid)*(ar_sim*rho <= min(grid)+delta/2)
        vn <- max(grid)*(ar_sim*rho >= max(grid)-delta/2)
        v_int <- sapply(grid_interno, function(x){x*(rho*ar_sim > x-delta/2)*(rho*ar_sim < x+delta/2)})

        return( v1 + rowSums(v_int) + vn )
        }
    }


    transic_cezar <- function(n, mu = 0, rho = 0.95, sigma = 0.007, seed = 27, grid_len = 9, tauchen = TRUE){ 

        tauchen_pdf <- tauchen(grid_len)[[1]]
        tauchen_cdf <- matrix(0, grid_len, grid_len)
        errors <- ar1(n)[[2]]

        tauchen_cdf[,1] <- tauchen_pdf[,1]
        for (i in 2:grid_len){
        tauchen_cdf[,i] <- rowSums(tauchen_pdf[,1:i])
        }

        grid[which(abs(grid-errors[1])==min(grid-errors[1]))]


        return(tauchen_cdf)
    }

transic_cezar(100)


## Idéia: simular Markov normal e a cada passo ver se o erro é alto o suficiente pra mudar o valro do theta.




max(between(0.95*ar1(10000)[[1]] + ar1(10000)[[2]], grid[3], grid[4]))




grid_len <- 9
rho <- 0.95
sigma <- 0.007
mu <- 0
m <- 3
n <- 10000
seed <- 27
grid
            
            
            erro <- ar1(n, mu = mu, rho = rho, sigma = sigma, seed = seed)$errors


            set.seed(seed)
            tauch <- tauchen(grid_len, mu = mu, rho = rho, sigma = sigma, m = m)
            grid <- tauch$grid
            probs <- tauch$probs

            sim <- c()
            sim[1] <- grid[which(abs(grid-erro[1])== min(abs(grid-erro[1])))]

            # A simulação escolhe theta_i como o valor do grid mais perto de E(theta_i|theta_{i-1}) + erro[i]
            for (i in 2:n){
                test <- sample(grid, 1,  prob = probs[which(grid == sim[i-1]),], replace = TRUE)
                sim[i] <- grid[which(abs(grid - erro[i] - test) == min(abs(grid - erro[i] - test)))]
            }

plot(ar1(10000, rho = rho)[[1]], type = 'l')
lines(sim, col = 'red')

