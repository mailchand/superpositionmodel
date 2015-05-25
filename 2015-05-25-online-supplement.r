#
# Implementation of the inverse Gaussian distribution
#
library(SuppDists) # install from CRAN first

#
# Multivariate Normal distribution
#
library(mvtnorm) # install from CRAN first

#
# Integrate dnorm(x) pnorm(a + bx) from -Inf to y (Owen, 1980, Eq. 10,010.1)
#
owen10_010.1 = function(a, b, y)
{
	rho = -b/sqrt(1 + b*b)
	pmvnorm(lower=c(-Inf, -Inf), upper=c(a/sqrt(1+b*b), y), 
		corr=matrix(c(1, rho, rho, 1), nrow=2))
}

#
# Integrate exp(ux) * dnorm(x | ma, sa) * pnorm(x | mb, sb) from -Inf to y
#
owen10_010.1x = function(u, ma, sa, mb, sb, y)
{
	ma1 = ma + u*sa*sa
	exp(u*ma + u*u*sa*sa/2) * owen10_010.1(a=(ma1 - mb)/sb, b=sa/sb, y=(y - ma1)/sa)
}

#
# Integrate x dnorm(x) pnorm(a + bx) from -Inf to y (Owen, 1980, Eq. 10,011.1)
#
owen10_011.1 = function(a, b, y)
{
	bb = sqrt(1 + b*b)
	b/bb * dnorm(a/bb) * pnorm(y*bb + a*b/bb) - dnorm(y)*pnorm(a + b*y)
}

#
# Integrate x exp(ux) dnorm(x | ma, sa) * pnorm(x | mb, sb) from -Inf to y
#
owen10_011.1x = function(u, ma, sa, mb, sb, y)
{
	ma1 = ma + u*sa*sa
	exp(u*ma + u*u*sa*sa/2) * 
		{sa * owen10_011.1(a=(ma1 - mb)/sb, b=sa/sb, y=(y - ma1)/sa) + 
			ma1 * owen10_010.1(a=(ma1 - mb)/sb, b=sa/sb, y=(y - ma1)/sa)}
}

#
# Numerically improved solution for exp(a) * pnorm(b) (Kiani et al., 2008)
#
exp_pnorm = function(a, b)
{
	r = exp(a) * pnorm(b)
	d = is.nan(r) & b < -5.5
	r[d] = 1/sqrt(2) * exp(a - b[d]*b[d]/2) * (0.5641882/b[d]/b[d]/b[d] - 1/b[d]/sqrt(pi))
	r
}

#
# Accuracy in unimodal and synchronous stimuli
#
acc_sync = function(d=865, c=100, mu=0.08, sigma2=10.5)
{
    pinvGauss(d, nu=c/mu, lambda=c*c/sigma2)
}

acc_sync(d=865, c=100, mu=0.08, sigma2=10.5)

#
# Accuracy in unimodal and synchronous stimuli (Eq. 8)
#
acc_async = function(d=865, c=100, mua=0.09, sigmaa2=11.3, mub=0.29, sigmab2=77.9, tau=240)
{
	# before tau
	p1 = pinvGauss(tau, nu=c/mua, lambda=c*c/sigmaa2)

	# after tau: Integral w(x, tau) * acc_sync(d-tau, c-x, mua+mub, sigmaa2+sigmab2)
	d_ = d - tau
	p2 = owen10_010.1x(u=0,
		ma=mua*tau, 
		sa=sqrt(sigmaa2*tau), 
		mb=c-(mua+mub)*d_, 
		sb=sqrt((sigmaa2+sigmab2)*d_), 
		y=c)

	p3 = exp(2*c*(mua+mub)/(sigmaa2+sigmab2)) * 
		owen10_010.1x(u=-2*(mua+mub)/(sigmaa2+sigmab2), 
		ma=mua*tau, 
		sa=sqrt(sigmaa2*tau), 
		mb=(mua+mub)*d_ + c, 
		sb=sqrt((sigmaa2+sigmab2)*d_), 
		y=c)

	p4 = exp(2*c*mua/sigmaa2) * 
		owen10_010.1x(u=0,
		ma=2*c + mua*tau, sa=sqrt(sigmaa2*tau), 
		mb=c - (mua+mub)*d_, 
		sb=sqrt((sigmaa2+sigmab2)*d_), 
		y=c)

	p5 = exp(2*c*mua/sigmaa2 + 2*c*(mua+mub)/(sigmaa2+sigmab2)) * 
		owen10_010.1x(u=-2*(mua+mub)/(sigmaa2+sigmab2), 
		ma=2*c + mua*tau, 
		sa=sqrt(sigmaa2*tau), 
		mb=c + (mua+mub)*d_, 
		sb=sqrt((sigmaa2+sigmab2)*d_), 
		y=c)

	# return value (only the first one is relevant, the others are reused for mrt_async)
	c(acc=p1 + p2 + p3 - p4 - p5, p2=p2, p3=p3, p4=p4, p5=p5)
}

acc_async(d=865, c=100, mua=0.09, sigmaa2=11.3, mub=0.29, sigmab2=77.9, tau=240)[1]

#
# Mean RT in unimodal and synchronous stimuli
#
mrt_sync = function(d=865, c=100, mu=0.09, sigma2=11.3, mum=336)
{
	# E(D | D < d) (Schwarz, 1994, Eq. 6)
	ed = c/mu * (pnorm(mu*d, c, sqrt(sigma2*d)) - exp_pnorm(2*c*mu/sigma2, (-c-mu*d)/sqrt(sigma2*d)))

	# Probability that detection occurs
	pd = acc_sync(d, c, mu, sigma2)

	# return value
	ed/pd + mum
}

mrt_sync(d=865, c=100, mu=0.09, sigma2=11.3, mum=336)

#
# Mean RT in asynchronous stimuli (Eq. 10)
#
mrt_async = function(d=865, c=100, mua=0.09, sigmaa2=11.3, mub=0.29, sigmab2=77.9, 
	mum=336, tau=240)
{
	# probability that target is detected
	p = acc_async(d, c, mua, sigmaa2, mub, sigmab2, tau)

	# H(tau) (Schwarz, 1994, Eq. 6)
	Htau = c/mua * (pnorm(mua*tau, c, sqrt(sigmaa2*tau)) - 
		exp_pnorm(2*c*mua/sigmaa2, (-c-mua*tau)/sqrt(sigmaa2*tau)))

	# tau * Integral w(x, tau) * G(d - tau)
	tauwG = tau * (p['p2'] + p['p3'] - p['p4'] - p['p5'])

	# c/muAV * (p2 - p3 - p4 + p5)
	minuends = c/(mua + mub) * (p['p2'] - p['p3'] - p['p4'] + p['p5'])

	# New terms
	d_ = d - tau
	q1 = 1/(mua + mub) * 
		owen10_011.1x(u=0,
		ma=mua*tau, sa=sqrt(sigmaa2*tau), 
		mb=c - (mua+mub)*d_, sb=sqrt((sigmaa2+sigmab2)*d_), y=c)

	q2 = 1/(mua + mub) * exp(2*c*(mua+mub)/(sigmaa2+sigmab2)) * 
		owen10_011.1x(u=-2*(mua+mub)/(sigmaa2+sigmab2), 
		ma=mua*tau, sa=sqrt(sigmaa2*tau), 
		mb=c + (mua+mub)*d_, sb=sqrt((sigmaa2+sigmab2)*d_), y=c)

	q3 = 1/(mua + mub) * exp(2*c*mua/sigmaa2) * 
		owen10_011.1x(u=0,
		ma=2*c + mua*tau, sa=sqrt(sigmaa2*tau), 
		mb=c - (mua+mub)*d_, sb=sqrt((sigmaa2+sigmab2)*d_), y=c)

	q4 = 1/(mua + mub) * exp(2*c*mua/sigmaa2 + 2*c*(mua+mub)/(sigmaa2+sigmab2)) * 
		owen10_011.1x(u=-2*(mua+mub)/(sigmaa2+sigmab2), 
		ma=2*c + mua*tau, sa=sqrt(sigmaa2*tau), 
		mb=c + (mua+mub)*d_, sb=sqrt((sigmaa2+sigmab2)*d_), y=c)

	# return value
	unname(1/p['acc'] * {Htau + tauwG + minuends - q1 + q2 + q3 - q4} + mum)
}

mrt_async(d=865, c=100, mua=0.09, sigmaa2=11.3, mub=0.29, sigmab2=77.9, 
	mum=336, tau=240)
