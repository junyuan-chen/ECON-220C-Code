mutable struct Estimates {TF<:AbstractFloat, TI<:Int}
    β_h::Array{TF}
    σ_β1::Array{TF}
    σ_β2::Array{TF}
    nReplic::TI
end

function Estimates(nReplic::Int)
    β_h = zeros(nReplic)
    σ_β1 = zeros(nReplic)
    σ_β2 = zeros(nReplic)
end

function



