# fractional-BM
Main source for fBM sim:
Georgiy Shevchenko, "Fractional Brownian Motion in a Nutshell"

We Want to find: a relationship between parameters of rough fractional stochastic volatility (see Volatility is rough - Gatheral) and Garch

How: fix rough fractional stochastic volatility parameters, use quasi maximum likelihood estimator to find garch params

Question to ask: 1. How well/poorly does it fit?
                 2. How complicated in some sort of rejection threshold?
                 3.Are there situations where rfsv can't be modelled by garch (analytic approximation for something whose vol is rfsv)?

Steps: 1. fgn paths -> rough fractional stochastic volatility paths
       2. quasi maximum likelihood with garch
       3. Further Steps? Document findings, paper etc

    
