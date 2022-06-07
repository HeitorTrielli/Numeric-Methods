# Heitor Trielli Zanoni Ferreira - Código - Lista 1
include("Tauchen.jl")
include("Rouwenhorst.jl")
include("Simular.jl")
using Plots, GLM, Pkg, DataFrames, PrettyTables, RegressionTables  # Pacotes que eu vou usar
using .Tauchen, .Rouwenhorst, .Simular # Módulos que estou usando

################
## Questão 1: ##
################

# A função que faz Tauchen está definida no arquivo "Tauchen.jl" e eu consigo chamar ela aqui por causa do using .Tauchen, que é o módulo definido no Tauchen.jl

tauchen_probs, tauchen_grid = tauchen(9);
tauchen_round = round.(tauchen_probs, digits = 3); # Arredondando para ficar mais legível


################
## Questão 2: ##
################

rouwenhorst_probs, rouwenhorst_grid = rouwenhorst(9);
rouwenhorst_round = round.(rouwenhorst_probs, digits = 3); # Arredondando para ficar mais legível


################
## Questão 3: ##
################
 
    # Simulando 
    ar_sim = ar1(10000)[1];
    tauch_sim = transic(10000);
    rouwen_sim = transic(10000, method = "rouwen");

    # Plotando para comparar
    plot([ar_sim, tauch_sim], label = ["AR(1)" "Tauchen"])
    plot([ar_sim, rouwen_sim], label = ["AR(1)" "Rouwenhorst"])

    # MQE
    pretty_table([mean((ar_sim - tauch_sim).^2) mean((ar_sim - rouwen_sim).^2)], header = ["Tauchen", "Rouwenhorst"])


################
## Questão 4: ##
################
    # Fazendo o dataframe para as regressões:
    df = DataFrame(Tauchen = tauch_sim, LagTauchen = lag(tauch_sim), Rouwenhorst = rouwen_sim, LagRouwen = lag(rouwen_sim));

    # Regressão do Tauchen
    tauchen_reg = lm(@formula(Tauchen ~ 0 + LagTauchen), df);

    # Regressão do Rouwenhorst
    rouwen_reg = lm(@formula(Rouwenhorst ~ 0 + LagRouwen), df);

    # Tabela com as regressões
    regtable(tauchen_reg, rouwen_reg)
    
    
################
## Questão 5: ##
################
    # Tauchen:
    # Grid e transição
    tauchen_grid_2 = tauchen(9, rho = 0.99)[2];
    tauchen_probs_2 = tauchen(9, rho = 0.99)[1];

    # Simulando o AR(1) e Tauchen
    ar_sim_2 = ar1(10000, rho = 0.99)[1];
    tauch_sim_2 = transic(10000, rho = 0.99);

    # Plotando
    plot([ar_sim_2 tauch_sim_2], label = ["AR(1)" "Tauchen"])

    # Rouwenhorst
    # Grid e transição
    rouwen_grid_2 = rouwenhorst(9, rho = 0.99)[2];
    rouwen_probs_2 = rouwenhorst(9, rho = 0.99)[1];

    # Simulação
    rouwen_sim_2 = transic(10000, rho = 0.99, method = "rouwen");

    # Plotando
    plot([ar_sim_2 rouwen_sim_2], label = ["AR(1)" "Rouwenhorst"])

    # Tabela com os MQE
    pretty_table([mean((ar_sim - tauch_sim_2).^2) mean((ar_sim - rouwen_sim_2).^2)], header = ["Tauchen", "Rouwenhorst"])

    # Regressões
    df_2 = DataFrame(Tauchen = tauch_sim_2, LagTauchen = lag(tauch_sim_2), Rouwenhorst = rouwen_sim_2, LagRouwen = lag(rouwen_sim_2));

    # Regressão do Tauchen
    tauchen_reg_2 = lm(@formula(Tauchen ~ 0 + LagTauchen), df_2);

    # Regressão do Rouwenhorst
    rouwen_reg_2 = lm(@formula(Rouwenhorst ~ 0 + LagRouwen), df_2);

    # Tabela com as regressões
    regtable(tauchen_reg_2, rouwen_reg_2)
