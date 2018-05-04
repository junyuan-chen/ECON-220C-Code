mutable struct Estimates{TF<:AbstractFloat, TI<:Int}
    ρ_OLS::Array{TF}
    ρ_FE::Array{TF}
    ρ_FD::Array{TF}
    ρ_AH::Array{TF}
    nReplic::TI
end

function Estimates(nReplic::Int)
    ρ_OLS::zeros(nReplic)
    ρ_FE::zeros(nReplic)
    ρ_FD::zeros(nReplic)
    ρ_AH::zeros(nReplic)
    return Estimates(ρ_OlS, ρ_FE, ρ_FD, ρ_AH, nReplic)
end

mutable struct Data{TF<:AbstractFloat, TI<:Int}
    mY::Array{TF}
    N::TI
    T::TI
end

function Data(N::Int, T::Int)
    mY = randn(N,T+1) # Store u_it
    return Data(mY, N, T)
end

function gen_data!(D::Data, ρ::AbstractFloat)
    vα = randn(D.N,1)
    D.mY[:,1] .= 0.5*vα .+ D.mY[:,1] # Set Y_0
    for t = 2:D.T+1
        D.mY[:,t] .= ρ.*D.mY[:,t-1] .+ vα .+ D.mY[:,t]
    end
end
