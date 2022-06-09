using Plots, Distributions, NLsolve, BenchmarkTools, Distributed, ProfileView, Roots, LinearAlgebra # Pacotes que estou usando

Threads.nthreads() # Quantas threads estamos usando

# Usando código das listas passadas
    tauchen = function(grid_len::Int64; mu::Float64 = 0.0, sigma::Float64 = 0.007, rho::Float64 = 0.95, m::Float64 = 3.0)

        theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
        theta_min = - theta_max + 2 * mu # definindo o menor valor do grid
        grid = Array(LinRange(theta_min, theta_max, grid_len))

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

        return (prob = probs::Array{Float64,2}, grid = grid)
    end;
#


beta, gamma, rho, sigma = 0.96,  1.0001, 0.9, 0.01
alpha = 1/3
zlen, alen = 9, 151

utility = function(c::Float64; mu::Float64 = gamma)
    return (c^(1 - mu) - 1)/(1 - mu)
end;
u_line = function(c, mu = 2.0)
    c^(-mu)
end;
u_inv = function(c, mu = 2.0)
    c^(-1/mu)
end;
######################################
############### item a ###############
######################################

w = exp.(tauchen(zlen, rho = rho, sigma = sigma).grid); # Valores que exp(z_t) pode tomar 
pz = tauchen(zlen, rho = rho, sigma = sigma).prob; # Matriz de transição de z

######################################
############### item b ###############
######################################

a_min = -10.0
a_max = 30.0
grid_a = Array(LinRange(a_min, a_max, alen))
a = copy(grid_a)
zmat = repeat(w, 1, alen)'
amat = repeat(grid_a, 1, zlen)

r = 0.038
c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)


#= r is the interest rate, value is the initial guess for the value function, alen is the length of the asset grid, zlen is the
length of the states grid, w is the wages grid, p_z is the Tauchen transition matrix, a is the asset grid=#

function value_function_brute(;r = r::Float64, value = v0::Matrix{Float64}, tol = 1e-4::Float64, inert = 0.0::Float64, 
    a_len = alen::Int64, z_len = zlen::Int64, w = w::Vector{Float64}, p_z = pz::Array{Float64, 2}, a = grid_a::Vector{Float64})

    #Pre alocating
    ArrMat::Type = Matrix{Array} # for convenience
    ro::ArrMat = Matrix{Array}(undef, a_len, z_len)
    # Since everything was undef a z, similar is also convenient
    a_possible::ArrMat = similar(ro)
    c_possible::ArrMat = similar(ro)
    val::ArrMat = similar(ro)
    v_possible::ArrMat = similar(ro)
    v_next, s = zeros(a_len, z_len), zeros(a_len, z_len)

    @sync @distributed for state in 1:z_len # setting what will remain unchanged
        ro[1:a_len, state] .= (w[state] + (1+r)*a[asset] .- a .> 0 for asset in 1:a_len) # Budget constraint (BC)
        a_possible[1:a_len, state] .= (a[w[state] + (1+r)*a[asset] .- a .> 0] for asset in 1:a_len) # Assets that obey BC
        c_possible[1:a_len, state] .= (w[state] + (1+r)*a[asset] .- a_possible[asset, state] for asset in 1:a_len) # Consumption that obeys BC
    end # for

    error = 1.0
    count = 0
    while error > tol # checking convergence
        if count%10 == 0 
            v_possible = Matrix{Array}(undef, a_len, z_len)
            @sync @distributed for state in 1:z_len
                v_possible[1:a_len, state] .= (value[ro[asset, state], :] for asset in 1:a_len) # Value function for each asset that obeys BC
                val[1:a_len, state] .= (utility.(c_possible[asset, state]) + beta*v_possible[asset, state]*p_z[state,:] for asset in 1:a_len) # Expected value of each of those assets
            end # for

            s = argmax.(val) # Taking the maximum
            a_line = a[s]
            @sync @distributed for state in 1:z_len 
                v_next[1:a_len, state] .= (val[asset, state][s[asset,state]] for asset in 1:a_len) # finding the next value function
            end # for
        else
            @sync @distributed for state in 1:z_len 
                v_next[1:a_len, state] .= (utility.(w[state] + (1+r)a[asset] - a_line[asset, state]) + (beta*value[findall(a_line[asset, state] .== a),:]*pz[state,:])[1] for asset in 1:a_len) # finding the next value function
            end 
        end
        error = maximum(abs.(value - v_next))
        value = (1 - inert)*v_next + inert*value # updating the value function
        count += 1
    end # while
    a_line = a[s] # Optimal policy function
    # Finding the optimal consumption
    zmat = repeat(w, 1, a_len)'
    amat = repeat(a, 1, z_len)
    c = zmat+(1+r)*amat - a_line
    return (val = value::Array{Float64,2}, pol = a_line::Array{Float64,2}, c = c::Array{Float64,2}) 
end; # function


solved = @time value_function_brute()
pol = solved.pol


pi_t = function(policy, pi_last; a_len = alen::Int64, z_len = zlen::Int64, p_z = pz::Array{Float64,2})
    Q = [p_z[:, state]'.*(policy.== grid_a[asset]) for asset in 1:a_len, state in 1:z_len]
    pi_t = [sum(pi_last.*Q[asset, state]) for asset in 1:a_len, state in 1:z_len]
    return pi_t::Matrix{Float64}

end

iterate_pi = function(policy; p_0 = 1, tol = 1e-6::Float64, pol = pol::Matrix{Float64}, a_len = alen)
    erro = 1.0
    if p_0 == 1
        pi_zero = pi_t(policy, ones(a_len, z_len)/(a_len*z_len))
    else
        pi_zero = p_0
    end
    iter = 0
    while erro > tol && iter < 150
        pi_one = pi_t(policy, pi_zero)
        erro = maximum(abs.(pi_one - pi_zero))
        pi_zero = 1*pi_one
        iter += 1
    end
    return pi_zero::Array{Float64,2}
end

statdist = @time iterate_pi(pol)
ed = sum(statdist.*pol)

rce = function(r_low, r_high; v0 = v0::Matrix{Float64}, tol = 1e-4)
    if r_low > r_high
        raux = copy(r_low)
        r_low = copy(r_high)
        r_high = raux
    end
        
    pol_low = @time value_function_brute(r = r_low ).pol
    pol_high = @time value_function_brute(r = r_high).pol
    
    statdist_low = @time iterate_pi(pol_low)
    statdist_high = @time iterate_pi(pol_high)

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
    display(plot(a, sum(statdist_mid, dims = 2), title = r_mid))
    println("ed = ", ed_low, "; ", ed_mid, "; ", ed_high, " Δr = ", r_high - r_low, " r = ", r_low, "; ", (r_high+r_low)/2, "; ", r_high)
    println(" ")

    if ed_mid > 0
        r_high = copy(r_mid)
        r_mid = (r_high+r_low)/2
        ed_high = ed_mid
    else
        r_low = copy(r_mid)
        r_mid = (r_low + r_high)/2
        ed_low = ed_mid
    end    

    while abs(ed_mid) > tol && r_high - r_low > 1e-12

        r_mid = (r_low+r_high)/2
        pol_mid = @time value_function_brute(r = r_mid).pol
        statdist_mid = @time iterate_pi(pol_mid)
        ed_mid = sum(statdist_mid.*pol_mid)
        display(plot(a, sum(statdist_mid, dims = 2), title = r_mid))
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
        println(" ")
    end
    return r_mid
end


rstar1 = @time rce(0.036, 0.038)
rstar1
solve_star = @time value_function_brute(r = rstar1)
statdist_star = @time iterate_pi(solve_star.pol)

##############################
########### item d ###########
##############################
rho = 0.97
r = 0.038

c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)
w = exp.(tauchen(z_len, rho = rho, sigma = sigma)[2]); # Valores que exp(z_t) pode tomar 
pz = tauchen(z_len, rho = rho, sigma = sigma)[1]; # Matriz de transição de z

pol2 = @time value_function_brute().pol
statdist2 = @time iterate_pi(pol2)
ed = sum(statdist2.*pol2)

rstar2 = rce(0.036, 0.038)

solve_star2 = @time value_function_brute(r = rstar2)
statdist_star2 = @time iterate_pi(solve_star2.pol)
# Deu um juros um pouco menor!

##############################
########### item e ###########
##############################
rho = 0.9
gamma = 5.0

r = 0.024
c0 = zmat+r*amat
v0 = utility.(c0)/(1-beta)

w = exp.(tauchen(z_len, rho = rho, sigma = sigma)[2]); # Valores que exp(z_t) pode tomar 
p_z = tauchen(z_len, rho = rho, sigma = sigma)[1]; # Matriz de transição de z
a = Array(LinRange(a_min, a_max, alen))

pol3 = @time value_function_brute().pol
statdist3 = @time iterate_pi(pol3)
ed = sum(statdist3.*pol3)

rstar3 = rce(0.024, 0.026)

solve_star3 = @time value_function_brute(r = rstar3)
statdist_star3 = @time iterate_pi(solve_star3.pol)
ed = sum(statdist_star3)
# deu menor que rstar2 



######## Gráficos ########
round(rstar1, digits = 4)
round(rstar2, digits = 4)
round(rstar3, digits = 4)



plot_caed = function(method::String)
    labels = ["State 1" "State 2" "State 3" "State 4" "State 5" "State 6" "State 7" "State 8" "State 9"]
    if method == "item c"
        c = plot(a, solved.c, label = labels, title = "Consumo baseline \n r = 0.038", legend = :topleft)
        val = plot(a, solved.val, label = labels, title = "Função valor baseline \n r = 0.38", legend = :topleft)
        pol = plot(a, solved.pol, label = labels, title = "Função politica baseline \n r = 0.38", legend = :topleft)
        dist = plot(a, sum(statdist, dims = 2), label = labels, title = "Consumo baseline \n r = 0.38", legend = :topleft)
    elseif method == "item d"
        c = plot(a, solve_star.c, label = labels, title = "Consumo baseline \n r = 0.0376", legend = :topleft)
        val = plot(a, solve_star.val, label = labels, title = "Função valor baseline \n r = 0.0376", legend = :topleft)
        pol = plot(a, solve_star.pol, label = labels, title = "Função politica baseline \n r = 0.0376", legend = :topleft)
        dist = plot(a, sum(statdist_star, dims = 2), label = labels, title = "Consumo baseline \n r = 0.0376", legend = :topleft)
    elseif method == "item e"
        c = plot(a, solve_star2.c, label = labels, title = "Consumo ρ = 0.97 \n r = 0.0374", legend = :topleft)
        val = plot(a, solve_star2.val, label = labels, title = "Função valor ρ = 0.97 \n r = 0.0374", legend = :topleft)
        pol = plot(a, solve_star2.pol, label = labels, title = "Função politica ρ = 0.97 \n r = 0.0374", legend = :topleft)
        dist = plot(a, sum(statdist_star2, dims = 2), label = labels, title = "Consumo ρ = 0.97 \n r = 0.0374", legend = :topleft)
    else
        c = plot(a, solve_star3.c, label = labels, title = "Consumo γ = 5 \n r = 0.0.0247", legend = :topleft)
        val = plot(a, solve_star3.val, label = labels, title = "Função valor γ = 5 \n r = 0.0247", legend = :topleft)
        pol = plot(a, solve_star3.pol, label = labels, title = "Função politica γ = 5 \n r = 0.0247", legend = :topleft)
        dist = plot(a, sum(statdist_star3, dims = 2), label = labels, title = "Consumo γ = 5 \n r = 0.0247", legend = :topleft)
    end
    return (val = val, c = c, pol = pol, dist = dist)
end

# plotando o item c
item_c = plot_caed("item c")
item_c.c
item_c.pol
item_c.dist
item_c.val

## plotando o item d
item_d = plot_caed("item d")
item_d.c
item_d.pol
item_d.dist
item_d.val

## plotando o item e
item_e = plot_caed("item e")
item_e.c
item_e.pol
item_e.dist
item_e.val

## plotando o item f
item_f = plot_caed("item f")
item_f.c
item_f.pol
item_f.dist
item_f.val
