using Random, Plots, GSL, Distributions, LinearAlgebra

autocov = (s,t,H) -> 1/2*(t^(2H) + s^(2H) - (t-s)^(2H))

function hosking(H, δt, n, m, gn)
end

function cholesky(H, δt, n, m, gn)
    
end

function davies_harte(H, δt, n, m, gn)
end

function fgn(H, T, n, m, method::Function; seed=nothing, kwargs...)
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
    m = 100
    n = 1
    seed = 1234
    method = hosking

end