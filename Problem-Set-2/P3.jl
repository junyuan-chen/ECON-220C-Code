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

function gen_data!(D::Data, ρ::AbstractFloat)
    vα = randn(D.N)
    D.mY[:,1] .= 0.5.*vα .+ D.mY[:,1] # Set Y_0
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

function simulate!(OLS::Estimates, FE::Estimates, FD::Estimates, AH::Estimates, N::Int, T::Int)
    for (i, ρ) in enumerate(OLS.vPara)
        for n = 1:OLS.nReplic
            D = Data(N,T)
            gen_data!(D, ρ)
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

###
function plot_stat(stats::Array, labels::Array; fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(2.2, 4))
    for (i, stat) in enumerate(stats)
        ax[:plot](0:0.1:1, stat, marker=".", linewidth = 1, label = labels[i])
    end
    xlabel(L"$\rho$")
    ax[:legend](loc = "upper right", frameon = true)
    tight_layout(pad = 0.1)
    savefig(string("Figure/", fname, ".pdf"))
end
###
function plot_est(est::Array, lab::String; fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(3.2, 3))
#    for (i, est) in enumerate(ests)
        ax[:hist](est, bins = 20, alpha = 0.8, density = true, label = lab)
        xlims = ax[:get_xlim]()
        @show rang = collect(linspace(xlims[1], xlims[2], 100))
        d = Normal(mean(est), std(est))
        ax[:plot](rang, pdf(d,rang))
#    end
    ax[:legend](loc = "upper right", frameon = true)
    tight_layout(pad = 0.1)
    #savefig(string("Figure/", fname, ".pdf"))
end
###
function plot_est(ests::Array, labs::Array; fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(5.2, 3))
    for (i, est) in enumerate(ests)
        ax[:hist](est, bins = 20, alpha = 0.8, label = labs[i])
    end
    ax[:legend](loc = "upper right", frameon = true)
    tight_layout(pad = 0.1)
    #savefig(string("Figure/", fname, ".pdf"))
end
#plot_est([OLSa.mρ[:,8], FEa.mρ[:,8], FDa.mρ[:,8]], ["OLS", "FE", "FD"])

###

OLSa, FEa, FDa, AHa = Estimates(), Estimates(), Estimates(), Estimates()

srand(10)
@time simulate!(OLSa, FEa, FDa, AHa, 100, 9)
plot_stat([OLSa.vBias, FEa.vBias, FDa.vBias, AHa.vBias], ["OLS", "FE", "FD", "AH"]; fname = "3a_bias")
plot_stat([OLSa.vSE, FEa.vSE, FDa.vSE, AHa.vSE], ["OLS", "FE", "FD", "AH"]; fname = "3a_se")
plot_stat([OLSa.vRMSE, FEa.vRMSE, FDa.vRMSE, AHa.vRMSE], ["OLS", "FE", "FD", "AH"]; fname = "3a_rmse")


plot_est(OLSa.mρ[:,8], "OLS")
plot_est(FEa.mρ[:,8], "FE")
plot_est(FDa.mρ[:,8], "FD")
plot_est(AHa.mρ[:,8], "AH")
