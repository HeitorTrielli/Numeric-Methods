using Plots, Distributions, NLsolve, BenchmarkTools, Distributed, ProfileView, Roots, LinearAlgebra # Pacotes que estou usando

Threads.nthreads() # Quantas threads estamos usando

# Usando código das listas passadas
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

        return (prob = probs::Array{Float64,2}, grid = Array(grid)::Array{Float64,1})
    end;
#

beta, gamma, rho, sigma = 0.96,  1.0001, 0.9, 0.01
alpha = 1/3
z_len, a_len = 9, 201

utility = function(c::Float64; mu::Float64 = gamma)
    return (c^(1 - mu) - 1)/(1 - mu)
end;


######################################
############### item a ###############
######################################

w = exp.(tauchen(z_len, rho = rho, sigma = sigma).grid); # Valores que exp(z_t) pode tomar 
p_z = tauchen(z_len, rho = rho, sigma = sigma).prob; # Matriz de transição de z

######################################
############### item b ###############
######################################

a_min = -10.0
a_max = 20.0
grid_a = Array(LinRange(a_min, a_max, a_len))
a = copy(grid_a)
zmat = repeat(w, 1, a_len)'
amat = repeat(grid_a, 1, z_len)

r = 0.040
c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)
a = Array(LinRange(a_min, a_max, a_len))


#= r is the interest rate, value is the initial guess for the value function, a_len is the length of the asset grid, z_len is the
length of the states grid, w is the wages grid, p_z is the Tauchen transition matrix, a is the asset grid=#
value_function_brute = function(;r = r::Float64, value = v0::Matrix{Float64}, tol = 1e-4::Float64, inert = 0.0::Float64, 
    a_len = a_len::Int64, z_len = z_len::Int64, w = w::Vector{Float64}, p_z = p_z::Array{Float64, 2}, a = a::Vector{Float64})

    #Pre alocating
    ro, a_possible, c_possible = Matrix{Array}(undef, a_len, z_len), Matrix{Array}(undef, a_len, z_len), Matrix{Array}(undef, a_len, z_len)
    val, v_possible = Matrix{Array}(undef, a_len, z_len), Matrix{Array}(undef, a_len, z_len) 
    v_next, s = zeros(a_len, z_len), zeros(a_len, z_len)

    @sync @distributed for state in 1:z_len # setting what will remain unchanged
        ro[1:a_len, state] .= (w[state] + (1+r)*a[asset] .- a .> 0 for asset in 1:a_len) # Budget constraint (BC)
        a_possible[1:a_len, state] .= (a[w[state] + (1+r)*a[asset] .- a .> 0] for asset in 1:a_len) # Assets that obey BC
        c_possible[1:a_len, state] .= (w[state] + (1+r)*a[asset] .- a_possible[asset, state] for asset in 1:a_len) # Consumption that obeys BC
    end # for

    error = 1.0

    while error > tol # checking convergence
        v_possible = Matrix{Array}(undef, a_len, z_len)
        @sync @distributed for state in 1:z_len
            v_possible[1:a_len, state] .= (value[ro[asset, state], :] for asset in 1:a_len) # Value function for each asset that obeys BC
            val[1:a_len, state] .= (utility.(c_possible[asset, state]) + beta*v_possible[asset, state]*p_z[state,:] for asset in 1:a_len) # Expected value of each of those assets
        end # for

        s = argmax.(val) # Taking the maximum

        @sync @distributed for state in 1:z_len 
            v_next[1:a_len, state] .= (val[asset, state][s[asset,state]] for asset in 1:a_len) # finding the next value function
        end # for
        error = maximum(abs.(value - v_next))
        value = (1 - inert)*v_next + inert*value # updating the value function
    end # while
    a_line = a[s] # Optimal policy function

    # Finding the optimal consumption
    zmat = repeat(w, 1, a_len)'
    amat = repeat(a, 1, z_len)
    c = zmat+(1+r)*amat - a_line
    return (val = value::Array{Float64,2}, pol = a_line::Array{Float64,2}, c = c::Array{Float64,2}) 
end; # function


solved = @time value_function_brute();
pol = solved.pol

pi_t = function(policy, pi_last; a_len = a_len::Int64, z_len = z_len::Int64, p_z = p_z::Array{Float64,2})
    Q = [p_z[:, state]'.*(policy.== grid_a[asset]) for asset in 1:a_len, state in 1:z_len]
    pi_t = [sum(pi_last.*Q[asset, state]) for asset in 1:a_len, state in 1:z_len]
    return pi_t::Matrix{Float64}

end

iterate_pi = function(policy; p_0 = 1, tol = 1e-8::Float64, pol = pol::Matrix{Float64})
    erro = 1.0
    if p_0 == 1
        pi_zero = pi_t(policy, ones(a_len, z_len)/(a_len*z_len))
    else
        pi_zero = p_0
    end
    iter = 0
    while erro > tol && iter < 100
        pi_one = pi_t(policy, pi_zero)
        erro = maximum(abs.(pi_one - pi_zero))
        pi_zero = 1*pi_one
        println(iter, ": erro = ", erro)
        iter += 1
    end
    return pi_zero::Array{Float64,2}
end

statdist = @time iterate_pi(pol)
ed = sum(statdist.*pol)
plot(a, sum(statdist, dims = 2))

rce = function(r_low, r_high; v0 = v0::Matrix{Float64}, tol = 1e-4)
    if r_low > r_high
        raux = copy(r_low)
        r_low = copy(r_high)
        r_high = raux
    end
        
    pol_low = value_function_brute(r = r_low ).pol
    pol_high = value_function_brute(r = r_high).pol
    
    statdist_low = iterate_pi(pol_low)
    statdist_high = iterate_pi(pol_high)

    ed_low = sum(pol_low.*statdist_low)
    ed_high = sum(pol_high.*statdist_high)

    if ed_low*ed_high > 0
        @error "Escolha taxas de juros que gerem excessos de demanda por ativo de sinais distintos"
    end

    ed_mid = 1
    r_mid = (r_low+r_high)/2
    r_mid = (r_low+r_high)/2
    pol_mid = @time value_function_brute(r = r_mid).pol
    statdist_mid = @time iterate_pi(pol_mid)
    ed_mid = sum(statdist_mid.*pol_mid)
    println("ed = ", ed_low, "; ", ed_mid, "; ", ed_high, " Δr = ", r_high - r_low, " r = ", r_low, "; ", (r_high+r_low)/2, "; ", r_high)

    if ed_mid > 0
        r_high = copy(r_mid)
        r_mid = (r_high+r_low)/2
        ed_high = ed_mid
    else
        r_low = copy(r_mid)
        r_mid = (r_low + r_high)/2
        ed_low = ed_mid
    end    

    while abs(ed_mid) > tol && r_high - r_low > 1e-8

        r_mid = (r_low+r_high)/2
        pol_mid = @time value_function_brute(r = r_mid).pol
        statdist_mid = @time iterate_pi(pol_mid)
        ed_mid = sum(statdist_mid.*pol_mid)
        display(plot(a, statdist_mid, title = r_mid))
        if ed_mid > 0
            r_high = copy(r_mid)
            r_mid = (r_high+r_low)/2
            ed_high = ed_mid
        else
            r_low = copy(r_mid)
            r_mid = (r_low + r_high)/2
            ed_low = ed_mid
        end    
        println("ed = ", ed_low, "; ", ed_mid, "; ", ed_high, " Δr = ", r_high - r_low, " r = ", r_low, "; ", r_mid, "; ", r_high)
    end
    return r_mid
end


rstar1 = @time rce(0.039, 0.04)

rstar1
solve_star = value_function_brute(r = rstar1)
statdist_star = @time iterate_pi(solve_star.pol)
plot(a, sum(statdist_star, dims = 2))

##############################
########### item d ###########
##############################
rho = 0.97

r = 0.039
c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)

w = exp.(tauchen(z_len, rho = rho, sigma = sigma)[2]); # Valores que exp(z_t) pode tomar 
p_z = tauchen(z_len, rho = rho, sigma = sigma)[1]; # Matriz de transição de z
a = Array(LinRange(a_min, a_max, a_len))

pol2 = @time value_function_brute(r = rstar1).pol
statdist2 = @time iterate_pi(pol2)
ed = sum(statdist2.*pol2)
plot(statdist2)

rstar2 = rce(0.039, 0.04)
# Deu o mesmo que rstar 1 na pratica com relação a juros, mas a distribuição parece mudar um pouco

##############################
########### item e ###########
##############################
rho = 0.9
gamma = 5.0

r = 0.03
c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)

w = exp.(tauchen(z_len, rho = rho, sigma = sigma)[2]); # Valores que exp(z_t) pode tomar 
p_z = tauchen(z_len, rho = rho, sigma = sigma)[1]; # Matriz de transição de z
a = Array(LinRange(a_min, a_max, a_len))

pol3 = @time value_function_brute().pol
statdist3 = @time iterate_pi(pol3)
ed = sum(statdist3.*pol3)
plot(statdist3)

rstar3 = rce(0.03, 0.034)
# deu menor que rstar1
