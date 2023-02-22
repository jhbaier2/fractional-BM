using Random, FFTW, Plots

function autocov(s, t, H)
    if s > t
        s, t = t, s
    end
    G = 2H
    γₕ = 1/2*(t^G + s^G - abs(t-s)^G)
    return γₕ
end

function autocov(n, H)
    # fbm auto covariance function
    G = 2H
    γₕ = 1/2*((n+1)^G + abs(n-1)^G - 2n^G)
    return γₕ
end

function Lambda(H,N)
    M = 2N-2
    C = zeros(Float64, (1,M))
    C[1:N] .= autocov.(0:N-1, H)
    C[N+1:M] .= reverse(C[2:N-1])
    Λ = sqrt.(real(fft(C, (2,))))
    return  Λ
end

function fractional_gaussian_noise(H::Float64, N::Int, T::Real; npath=1, seed=nothing, kwargs...)
    # fractional gaussian noise
    rng = MersenneTwister(seed)
    Gn = randn(rng, (npath, 2N-2))
    δt = T/N
    if H==0.5
        fGn = Gn[:,1:N]
    else
        Λ = Lambda(H,N)
        M = size(Λ, 2)
        b = ifft(Gn, (2,))
        a = b.*Λ
        fGn = real(fft(a,2))[:, 1:Int(M/2)]
    end
    return transpose(δt .* fGn)
end

function fractional_gaussian_noise(Λ::Matrix, NT::Int, T::Real; seed=nothing)
    δt = T/N
    rng = MersenneTwister(seed)
    M = size(Λ, 2)
    Gn = randn(rng, (NT, M))
    b = ifft(Gn, (2,))
    a = b.* Λ
    fGn = real(fft(a,2))[:, 1:Int(M/2)]
    return transpose(δt .* fGn)
end

function fractional_brownian_motion(H, T, N; npath=1, seed=nothing, kwargs...)
    ### tbd
    rng = MersenneTwister(seed)
    Gn = randn(rng, (npath, N))
    if H==0.5
        fGn = Gn
    else
    end
    fBm = cumsum(fGn; dims=1)
    return fBm
end