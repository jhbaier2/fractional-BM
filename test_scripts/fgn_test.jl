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
    Λ = sqrt.(real(fft(C)))
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

function fgn(Λ, NT)

    rng = MersenneTwister(0) #matlab default rng
    M = size(Λ, 2)
    rands = randn(NT, M)
    b = ifft(rands)
    a = broadcast(*, b, Λ)
    fGn = real(fft(a,2))[:, 1: Int(M/2)]
    return fGn  
end

begin
    q = 2
    N = 2^q + 1
    H = 0.3
    lam = Λ(H,N)
    fgn(lam, 20)
end




#=
-------------------TESTS--------------- -> put in test file
#test fgn
begin 
    rands = [[-1.3239,	0.7503,	-0.7004,	0.8401,	-1.0852,	0.0909,	-0.4381,	-1.784,]
    [-0.2794,	-0.5158,	-0.8662,	-1.2097,	1.1274,	1.0755,	-0.7616,	1.3512,]
    [0.7854,	1.3949,	-0.2259,	1.8403,	-0.7558,	-0.111,	0.2349,	0.9666,]
    [-0.9421,	-1.3855,	-1.2795,	-0.6421,	-0.5374,	-1.2935,	2.2342,	-0.1436,]
    [0.6595,	1.7095,	0.7379,	0.7423,	-0.8658,	0.6412,	0.7097,	0.2791]]
    # use rands in the fgn nfunction for variable "rands"

    q = 2
    N = 2^q + 1
    H = 0.3
    lam = Λ(H,N)
    fgn(lam, 20)

    # matlab outputs std(fGn) = 1.2170, mean = 0.0074, max = 



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