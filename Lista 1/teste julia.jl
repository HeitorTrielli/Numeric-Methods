using Distributions, Pkg, StatsBase

Pkg.add("StatsBase")


ar1 = function(n; mu = 0, rho = 0.95, sigma = 0.007, seed = 27) 
    Random.seed!(seed)
    errors = rand(Normal(0, sigma), n) # gerando vetor de erros
    sample = zeros(0)
    append!(sample, errors[1])

    for i in 2:n
        append!(sample, mu + rho*sample[i-1] + errors[i])
    end
    
    return sample, errors
end

rouwenhorst = function (grid_len; mu = 0, sigma = 0.007, rho = 0.95)
    theta_max = (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
    theta_min = - theta_max
    grid = LinRange(theta_min, theta_max, grid_len)

    p1 = (1+rho)/2
    p = [p1 (1 - p1) ; (1 - p1) p1]
    if grid_len > 2
        for i in 3:grid_len
            z1 = [[p zeros(i-1)]; zeros(i)']
            z2 = [[zeros(i-1) p]; zeros(i)']
            z3 = [zeros(i)'; [p zeros(i-1)]]
            z4 = [zeros(i)'; [zeros(i-1) p]]

            p = p1*(z1 + z4) + (1-p1)*(z2+z3)
        end # for
    end # if
    
    transition_matrix = p./(sum(p, dims = 2))
    
    return transition_matrix, grid
end # function


sigma, rho, grid_len = 0.007, 0.95, 9
m = 3


theta_max = (sigma/sqrt(1-rho^2)) * sqrt(grid_len-1)
theta_min = - theta_max
grid = LinRange(theta_min, theta_max, grid_len)

p1 = (1+rho)/2
p = [p1 (1 - p1) ; (1 - p1) p1]
if grid_len > 2
    for i in 3:grid_len
        z1 = [[p zeros(i-1)]; zeros(i)']
        z2 = [[zeros(i-1) p]; zeros(i)']
        z3 = [zeros(i)'; [p zeros(i-1)]]
        z4 = [zeros(i)'; [zeros(i-1) p]]

        p = p1*(z1 + z4) + (1-p1)*(z2+z3)
    end # for
end # if

transition_matrix = p./(sum(p, dims = 2))
        



round_rouwen = round.(transition_matrix, digits = 3)

teste = [1 2 3; 4 5 6; 7 8 9]

teste[findmin(abs.(teste .- 3))]


mu = 0

theta_max = m * sigma / (sqrt(1 - rho^2)) + mu # definindo o maior valor do grid
theta_min = - theta_max # definindo o menor valor do grid
grid = LinRange(theta_min, theta_max, grid_len) # Cria um vetor de n pontos entre theta_max e theta_min em que a distancia entre os pontos sequencias é igual

d = Normal(mu, 1) # d vira normal(mu,1), que será usado para computar a PDF dos erros na hora de achar as probabilidades de transição
delta = (maximum(grid) - minimum(grid)) / (length(grid) - 1) # distância dos pontos subsequentes do grid

vec_1 = cdf(d,((minimum(grid) .- rho * grid .+ delta / 2) / sigma)) # vetor das probabilidades de ir para o menor valor do grid, dado cada estado anterior do grid; cdf(d, x) retorna a cdf da distribuição d no valor x
vec_n = 1 .- cdf(d,((maximum(grid) .- rho * grid .- delta / 2) / sigma)) # análogo para o maior valor do grid
grid_interno = grid[2:(length(grid) - 1)] # valores não extremos do grid

pij = function(j, i = grid) # função que vai computar o vetor de probabilidades de ir para o estado (não extremo) j dado cada estado anterior do grid
    cdf(d,((j .+ delta/2 .- rho * i) / sigma)) - cdf(d,((j .- delta / 2 .- rho * i) / sigma))                             
end


vec_1
vec_n

mat_interna = reduce(hcat, map(pij, grid_interno)) # map: aplica pij em cada ponto do grid interno; reduce: transforma o vetor de vetores que vem do map em uma matriz

probs_tauch = [vec_1 mat_interna vec_n] # gerando a matriz de transição

return probs, grid


n=10000

erro = ar1(n)[2]

sim = zeros(n)


sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]]

for i in 2:(n-1)
    sim[i] = grid[findmin(abs.(grid .- erro[2] .- probs[findall(x -> x == sim[1], grid),:]*grid))[2]]

    print(sim[i] )
end #for


prob, grid = rouwenhorst(grid_len, mu = mu, rho = rho, sigma = sigma)

sim = zeros(n)
sim[1] = grid[findmin(abs.(grid .- erro[1]))[2]]

for i in 2:n
    sim[i] = grid[findmin(abs.(grid .- erro[i] .- prob[findall(x -> x == sim[i-1], grid),:]*grid))[2]]
end # for
