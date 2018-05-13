using PyPlot
rc("font", size=9)
using LaTeXStrings

function plot_qs(Q::Function, L::AbstractFloat, R::AbstractFloat; N::Int = 200, fname::String = "no_name")
    matplotlib[:style][:use]("seaborn-whitegrid")
    fig, ax = subplots(figsize=(4, 3.5))
    X = linspace(L, R, N)
    for i = 1:2:10 #1:5:I*5
        F = x -> Q(x,i)
        ax[:plot](X, F.(X), linewidth = 1, label = "n = $i ")
        println(minimum(F.(X)))
    end
    xlabel(L"$\theta$")
    ylabel(L"$Q_n(\theta)$")
    ax[:legend](frameon = true)
    tight_layout(pad = 0.1)
    savefig(string("Figure/", fname, ".pdf"))
end

function Q(x::AbstractFloat, n::Int)
    if x >= -2 && x <= -1
        return (2.0^(1-n)-1)*x + 2.0^(2-n)-2
    elseif x > -1 && x <= 0
        return (1-2.0^(1-n))*x
    elseif x > 0 && x <= 1-2.0^(-n-1)
        return (2.0-2.0^(n+1))/(2.0^(n+1)-1)*x
    elseif x > 1-2.0^(-n-1) && x < 1
        return (2-2.0^(n+1))/(2.0^(n+1)-1)*x + 2.0^(1-n) -2
    elseif x == 1
        return 0
    else
        throw(DomainError())
    end
end

plot_qs(Q, -2.0, 1.0; fname = "P1c")
