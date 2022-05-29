
d = Int(length(gamma)/7)
root = chebroots(d)
# Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de chebyschev de grau d
# Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
k0level = k0*(k .- k_ss)[500].+k_ss

# Aplica a função cons (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
# (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
c0 = reduce(hcat, [cons.((gamma,), k0, 1:7, d) for k0 in k0])
k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[500]) for j in 1:d for i in 1:7 ], 7, 2)

c1 = reinterpret(Float64, reshape([cons.((gamma,), k1[i,j], i, d) for j in 1:d for i in 1:7], 7, d))
resid = u_line.(c0)' - beta*(u_line.(c1).*z_grid.*(k0level.^(alpha-1)*alpha .+ (1 - delta))')'p_z'















z_grid = grid_z


gamma = ones(7,2)

d = Int(length(gamma)/7)
root = chebroots(d)
# Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de chebyschev de grau d
# Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
k0level = k0*(k .- k_ss)[500].+k_ss


# Aplica a função cons (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
# (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
c0 = reduce(hcat, [cons.((gamma,), k0, 1:7, d) for k0 in k0])
k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
z_grid.*(k0level.^alpha)'

k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[500]) for j in 1:d for i in 1:7 ], 7, 2)
k1 = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:7], 7, d))[i, j]) for j in 1:d for i in 1:7], 7, d))

c1 = reinterpret(Float64, reshape([cons.((gamma,), (k1[i,j]-k_ss)/((k .- k_ss)[500]), i, d) for j in 1:d for i in 1:7], 7, d))
resid = u_line.(c0)' - beta*(u_line.(c1).*z_grid.*(k0level.^(alpha-1)*alpha .+ (1 - delta))')'p_z'


d = Int(length(gamma)/7)
root = chebroots(d)
# Escolhe os k0 como os valores de k normalizado que estão mais pertos das raizes do polinomio de chebyschev de grau d
# Para cada raiz, ele acha o índice do argmin do valor absoluto do vetor de k normalizado - a raiz.
k0 = knorm[[argmin(abs.(knorm .- root)) for root in root]]
k0level = k0*(k .- k_ss)[500].+k_ss

# Aplica a função cons (avaliada direto nos 7 estados) para cada k0. Isso retorna um array de arrays, por isso o reduce(hcat, ...)
# (gamma,) serve para tratar a matriz como uma tupla e poder aplicar a função direto em todos os estados de uma vez (não me pergunte, eu também não entendo)
c0 = reduce(hcat, [cons.((gamma,), k0, 1:7, d) for k0 in k0])
k1level = z_grid.*(k0level.^alpha)' - c0 .+ ((1 - delta)*k0level)'
k1norm = reshape([(k1level[i,j]-k_ss)/((k .- k_ss)[500]) for j in 1:d for i in 1:7 ], 7, 2)
k1 = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:7], 7, d))[i, j]) for j in 1:d for i in 1:7], 7, d))

c1 = reinterpret(Float64, reshape([cons.((gamma,), k1[i,j], i, d) for j in 1:d for i in 1:7], 7, d))
resid = u_line.(c0)' - beta*(u_line.(c1).*z_grid.*(k0level.^(alpha-1)*alpha .+ (1 - delta))')'p_z'
return resid







k1 = maximum.(reshape([hcat(-1, minimum.(reshape([hcat(1, k1norm[i, j]) for j in 1:d for i in 1:7], 7, d))[i, j]) for j in 1:d for i in 1:7], 7, d))
