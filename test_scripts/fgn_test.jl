using Random, Plots, FFTW, Distributions, LinearAlgebra, Test, Statistics

function autocov(n, H)
    G = 2H
    γₕ = 1/2*((n+1)^G + abs(n-1)^G - 2n^G)
    return γₕ
end

function Λ(H,N)\
    M = 2N-2
    C = zeros(Float64, (1,M))
    C[1:N] .= autocov.(0:N-1, H)
    C[N+1:M] .= reverse(C[2:N-1])
    Λ = sqrt.(real(fft(C, (2,))))
    return  Λ
end

function autocov(s, t, H)
    if s > t
        s, t = t, s
    end
    G = 2H
    γₕ = 1/2*(t^G + s^G - abs(t-s)^G)
    return γₕ
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

function fractional_guassian_noise(Λ, NT)

    rng = MersenneTwister(0) #matlab default rng
    M = size(Λ, 2)
    rands = randn(rng, (NT, M))
    #rands = [0.5377	-1.3077	-1.3499	-0.205	0.6715	1.0347	0.8884	1.4384; 1.8339	-0.4336	3.0349	-0.1241	-1.2075	0.7269	-1.1471	0.3252; -2.2588 0.3426	0.7254	1.4897	0.7172	-0.3034	-1.0689	-0.7549; 0.8622 3.5784	-0.0631	1.409	1.6302	0.2939	-0.8095	1.3703; 0.3188 2.7694	0.7147	1.4172	0.4889	-0.7873	-2.9443	-1.7115]
    b = ifft(rands, (2,))
    a = b.* Λ
    fGn = real(fft(a,2))[:, 1: Int(M/2)]
    return fGn
end

begin
    q = 10
    N = 2^q + 1
    H = 0.7
    lam = Λ(H,N)
    fgn_i = fractional_guassian_noise(lam, 5)
    fbm = cumsum(fgn_i; dims = 2)
    plot(fbm')
end

#=
-------------------TESTS--------------- -> put in test file
#test fgn
begin 
     # rands = [0.5377	-1.3077	-1.3499	-0.205	0.6715	1.0347	0.8884	1.4384; 1.8339	-0.4336	3.0349	-0.1241	-1.2075	0.7269	-1.1471	0.3252; -2.2588 0.3426	0.7254	1.4897	0.7172	-0.3034	-1.0689	-0.7549; 0.8622 3.5784	-0.0631	1.409	1.6302	0.2939	-0.8095	1.3703; 0.3188 2.7694	0.7147	1.4172	0.4889	-0.7873	-2.9443	-1.7115]
    # use rands in the fgn nfunction for variable "rands"

    q = 2
    N = 2^q + 1
    H = 0.3
    lam = Λ(H,N)
    fgm = cumsum(fgn(lam, 20)





#Test autocov 
begin
    n=10
    H=.05
    ac = autocov(n, H)
    @test ac ≈ 0.5837 atol=0.001
end

#Test  Λ 
begin
    q = 10
    N = 2^q + 1
    H = 0.3
    lam = Λ(H,N)
    @test std(lam) ≈ 0.1984 atol=0.001 
    #0.1984 is a matlab check from Shevchenko
end
=#