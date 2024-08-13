using Random, FFTW

function cov(H, t)  
	G = 2H
	c = (t+1)^G + abs(t-1)^G - 2t^G
	return c/2
end

function Lambda(H, N)
	M = 2N - 2
	C = zeros(M)
	C[1:N] = cov.(H, collect(0:N-1))
	C[N+1:M] = C[1:N][N-1:-1:2]
	res = real(fft(C))
	return sqrt.(res)
end

function fractional_gaussian_noise(H, n, t; nsim=1, seed=1234)
    Λ = Lambda(H, n)
    M = size(Λ)[1]
	rng = MersenneTwister(seed)
	gn = randn(rng, (M, n))
	λi = ifft(gn)
	Q = λi .* Λ
	ξ = real(fft(Q))
	δt = t/(M/2)
	idx = 1:M÷2
	return ξ[idx] .* δt^H
end;

function fractional_brownian_motion(H, n, t; F₀=0., nsim=1, seed=1234)
    fgn = fractional_gaussian_noise(H, n, t; nsim=1, seed=1234)
    fgn = append!([F₀], fgn)
    fbm = cumsum(fgn)
    return fbm
end

begin
	H = 0.1
	q = 10
	n = 2^q + 1
	t = 1
    fractional_gaussian_noise(H,n,t)
end;