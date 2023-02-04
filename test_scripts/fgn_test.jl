using Random, Plots, GSL, Distributions, LinearAlgebra

function autocov(δt, H)
    if δt<0
        δt = abs(δt)
    elseif δt==0
        return 1.
    else
        γₖ = 1/2*((δt+1)^(2H) + abs(δt-1)^(2H) - (2δt)^(2H))
        return γₖ
    end
end


function fgn(H, T, n, m; seed=nothing, kwargs...)
    rng = MersenneTwister(seed)
    δt = T/n
    Gn = randn(rng, (m,n))
    if H==0.5
        fGn = Gn
    else
        
    end
    return fGn
end

begin
    m = 10
    n = 1
    seed = 1234
    H = 0.3

end