# --------------------------------------------------------------------------------
# Update: 30 June 2017
# To estimate parameters using conditional second-order estimating equation.
# --------------------------------------------------------------------------------

source('sim.CGEE2.structured.family.r')

# With the Two-generation family data, just run function `script.addprob.str.f'
# But need to create a data frame for 'stats' including: lam, kap, beta and gam.vec; where lam is the scale parameter in Weibull distribution, kap is the shape parameter, beta is the coefficient for a variable X (binary here), and gam.vec is the vector in the second-order model for within-family association (see Zhong and Cook, Statistics in Biosciences (2017)).
# --------------------------------------------------------------------------------
# A toy example to generate family data
stats <- getstats.str.f(kap=1.2, beta=log(1.2), Ktau.pp=0.1, Ktau.ss=0.4, Ktau.ps=0.2, ad.cens=90, med.age=45)

simdata <- generate.str.f(nsize = 500, stats=stats)

# --------------------------------------------------------------------------------
# To estimate parameters based on the proposed estimating equations.
# Three options for der.method : full, simple1 and simple2; two options for cov.method: full and simple2 (WPI). See Zhong and Cook, Statistics in Biosciences (2017) for details.

est <- script.addprob.str.f(instats=stats, indata=simdata, der.method='simple1', cov.method='simple2')



