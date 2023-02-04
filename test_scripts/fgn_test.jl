using Random, Plots, FFTW, Distributions, LinearAlgebra

function autocov(s, t, H)
    if s > t
        s, t = t, s
    end
    G = 2H
    γₕ = 1/2*(t^G + s^G - abs(t-s)^G)
    return γₕ
end

function autocov(n, H)
    G = 2H
    γₕ = 1/2*((n+1)^G + abs(n-1)^G - (2n)^G)
    return γₕ
end

function Λ(H,N)
    M = 2N - 2
    C = zeros(Float64, (1,M))
    C[1:N] .= autocov.(0:N-1, H)
    C[N+1:M] .= reverse(C[2:N-1], dims=2)
    return sqrt.(real(fft(C)))
end

function fgn(H, T, N; npath=1, seed=nothing, kwargs...)
    # fractional gaussian noise
    rng = MersenneTwister(seed)
    δt = T/n
    Gn = randn(rng, (N,npath))
    if H==0.5
        fGn = Gn
    else
        
    end
    return fGn
end

begin
    q = 10
    N = 2^q + 1
    npath = 1
    seed = 1234
    H = 0.3

end