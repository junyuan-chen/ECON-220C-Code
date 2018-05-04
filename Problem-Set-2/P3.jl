mutable struct Estimates{TF<:AbstractFloat, TI<:Int}
    nReplic::TI
    ρ_OLS::Array{TF,1}
    ρ_FE::Array{TF}
    ρ_FD::Array{TF}
    ρ_AH::Array{TF}
end

function Estimates(nReplic::Int)
    ρ_OLS = zeros(nReplic)
    ρ_FE = zeros(nReplic)
    ρ_FD = zeros(nReplic)
    ρ_AH = zeros(nReplic)
    return Estimates(nReplic, ρ_OLS, ρ_FE, ρ_FD, ρ_AH)
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


function simulate(Est::Estimates, N::Int, T::Int, ρ::AbstractFloat)
    for i = 1:Est.nReplic
        D = Data(N,T)
        gen_data!(D, ρ)
        Est.ρ_OLS[i] = pooling_OLS(D)
        Est.ρ_FE[i] = fixed_effects(D)
    end
end

N = 1000
Est1 = Estimates(N)

srand(10)
simulate(Est1, 100, 6, 0.1)
