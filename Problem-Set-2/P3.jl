using PyPlot
rc("font", size=9)
using LaTeXStrings
using Distributions

mutable struct Estimates{TF<:AbstractFloat, TI<:Int}
    nReplic::TI
    nPara::TI
    vPara::Array{TF,1}
    mρ::Array{TF,2}
    vBias::Array{TF,1}
    vSE::Array{TF,1}
    vRMSE::Array{TF,1}
end

function Estimates(; nReplic::Int = 1000,
    vPara = collect(0:0.1:1))

    nPara = length(vPara)
    mρ = zeros(nReplic, nPara)
    vBias = zeros(nPara)
    vSE = zeros(nPara)
    vRMSE = zeros(nPara)
    return Estimates(nReplic, nPara, vPara, mρ, vBias, vSE, vRMSE)
end

mutable struct Data{TF<:AbstractFloat, TI<:Int}
    N::TI
    T::TI
    mY::Array{TF}
end

function Data(N::Int, T::Int)
    mY = randn(N,T+1) # Store u_it
    return Data(N, T, mY)
end

function gen_data!(D::Data, ρ::AbstractFloat, α::AbstractFloat)
    vα = randn(D.N)
    D.mY[:,1] .= α .+ 0.5.*vα .+ D.mY[:,1] # Set Y_0
    for t = 2:D.T+1
        D.mY[:,t] .= ρ.*D.mY[:,t-1] .+ vα .+ D.mY[:,t]
    end
end

function pooling_OLS(D::Data)
    vY = D.mY[:, 2:D.T+1]
    vY = reshape(vY, length(vY), 1)
    vY_lag = D.mY[:, 1:D.T]
    vY_lag = reshape(vY_lag, length(vY_lag), 1)
    return (vY\vY_lag)[1]
end

function fixed_effects(D::Data)
    mY = D.mY .- mean(D.mY, 2)
    vY = (mY[:, 2:D.T+1])[:]
    vY_lag = (mY[:, 1:D.T])[:]
    return (vY\vY_lag)[1]
end

function first_difference(D::Data)
    mY = D.mY[:, 2:D.T+1] .- D.mY[:, 1:D.T]
    vY = (mY[:, 2:D.T])[:]
    vY_lag = (mY[:, 1:D.T-1])[:]
    return (vY\vY_lag)[1]
end

function anderson_hsiao(D::Data)
    mY = D.mY[:, 2:D.T+1] .- D.mY[:, 1:D.T]
    vY = (mY[:, 2:D.T])[:]
    vY_lag = (mY[:, 1:D.T-1])[:]
    vZ = D.mY[:, 1:D.T-1]
    return dot(vZ, vY)/dot(vZ, vY_lag)
end

function bias(vEst::Array, para::AbstractFloat)
    return mean(vEst) - para
end

function se(vEst::Array)
    Avg = mean(vEst)
    return sqrt(mean(vEst.^2) - Avg^2)
end

function simulate!(OLS::Estimates, FE::Estimates, FD::Estimates, AH::Estimates, N::Int, T::Int; α::AbstractFloat = 0.0)
    for (i, ρ) in enumerate(OLS.vPara)
        for n = 1:OLS.nReplic
            D = Data(N,T)
            gen_data!(D, ρ, α)
            OLS.mρ[n,i] = pooling_OLS(D)
            FE.mρ[n,i] = fixed_effects(D)
            FD.mρ[n,i] = first_difference(D)
            AH.mρ[n,i] = anderson_hsiao(D)
        end
        for Est in [OLS, FE, FD, AH]
            Est.vBias[i] = bias(Est.mρ[:,i], ρ)
            Est.vSE[i] =  se(Est.mρ[:,i])
            Est.vRMSE[i] = sqrt(Est.vBias[i]^2 + Est.vSE[i]^2)
        end
    end
end

function plot_stat(stats::Array, labels::Array; fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(2.2, 4))
    for (i, stat) in enumerate(stats)
        ax[:plot](0:0.1:1, stat, marker=".", linewidth = 1, label = labels[i])
    end
    xlabel(L"$\rho$")
    ax[:legend](loc = "best", frameon = true)
    tight_layout(pad = 0.1)
    savefig(string("Figure/", fname, ".pdf"))
    close(fig)
end

function plot_est(est::Array, lab::String; fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(3.2, 3))
    ax[:hist](est, bins = 20, alpha = 0.8, density = true, label = lab)
    xlims = ax[:get_xlim]()
    rang = collect(linspace(xlims[1], xlims[2], 100))
    d = Normal(mean(est), std(est))
    ax[:plot](rang, pdf.(d,rang))
    ax[:legend](loc = "best", frameon = true)
    tight_layout(pad = 0.1)
    savefig(string("Figure/", fname, ".pdf"))
    close(fig)
end

OLSa, FEa, FDa, AHa = Estimates(), Estimates(), Estimates(), Estimates()
OLSb3, FEb3, FDb3, AHb3 = Estimates(), Estimates(), Estimates(), Estimates()
OLSb9, FEb9, FDb9, AHb9 = Estimates(), Estimates(), Estimates(), Estimates()
OLSc, FEc, FDc, AHc = Estimates(), Estimates(), Estimates(), Estimates()

OLS = [OLSa, OLSb3, OLSb9, OLSc]
FE = [FEa, FEb3, FEb9, FEc]
FD = [FDa, FDb3, FDb9, FDc]
AH = [AHa, AHb3, AHb9, AHc]
T = [6, 3, 9, 18]
αs = [0.0, 0.0, 0.0, 5.0]
fnames_bias = ["3a_bias", "3b3_bias", "3b9_bias", "3c_bias"]
fnames_se = ["3a_se", "3b3_se", "3b9_se", "3c_se"]
fnames_rmse = ["3a_rmse", "3b3_rmse", "3b9_rmse", "3c_rmse"]
fnames_hist_OLS = ["3a_OLS", "3b3_OLS", "3b9_OLS", "3c_OLS"]
fnames_hist_FE = ["3a_FE", "3b3_FE", "3b9_FE", "3c_FE"]
fnames_hist_FD = ["3a_FD", "3b3_FD", "3b9_FD", "3c_FD"]
fnames_hist_AH = ["3a_AH", "3b3_AH", "3b9_AH", "3c_AH"]

for i = 4#1:length(OLS)
    srand(10)
    @time simulate!(OLS[i], FE[i], FD[i], AH[i], 100, T[i]; α = αs[i])
    plot_stat([OLS[i].vBias, FE[i].vBias, FD[i].vBias, AH[i].vBias], ["OLS", "FE", "FD", "AH"]; fname = fnames_bias[i])
    plot_stat([OLS[i].vSE, FE[i].vSE, FD[i].vSE, AH[i].vSE], ["OLS", "FE", "FD", "AH"]; fname = fnames_se[i])
    plot_stat([OLS[i].vRMSE, FE[i].vRMSE, FD[i].vRMSE, AH[i].vRMSE], ["OLS", "FE", "FD", "AH"]; fname = fnames_rmse[i])
    plot_est(OLS[i].mρ[:,8], "OLS", fname = fnames_hist_OLS[i])
    plot_est(FE[i].mρ[:,8], "FE", fname = fnames_hist_FE[i])
    plot_est(FD[i].mρ[:,8], "FD", fname = fnames_hist_FD[i])
    plot_est(AH[i].mρ[:,8], "AH", fname = fnames_hist_AH[i])
end
