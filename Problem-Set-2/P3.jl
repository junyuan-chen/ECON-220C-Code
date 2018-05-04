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


OLSa, FEa, FDa, AHa = Estimates(), Estimates(), Estimates(), Estimates()

srand(10)
@time simulate!(OLSa, FEa, FDa, AHa, 100, 6)
