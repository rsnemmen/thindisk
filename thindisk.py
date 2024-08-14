"""
Thin disk spectrum
======================
"""
import numpy, pylab
import scipy.integrate
import sys


def thindisk(m,mdot,rtran,angle,mfile=None,n=50,nui=10.,nuf=17.,rout=1e5):
	"""
Computes the spectrum of a standard Shakura-Sunyaev thin accretion disk.
The units of the code are CGS.

Usage:
>>> lognu,ll=thindisk(4e6,1e-3,3.,60.,'model.dat')

Arguments:
- m : black hole mass in solar masses
- mdot : mass accretion rate in Eddington units
- rtran : inner radius of the disk in Schwarzschild radii
- angle : inclination angle in degrees
- mfile : output file
- n : desired number of points in the resulting SED, default=50
- nui,nuf : initial, final frequencies desired, default 10 and 17
- rout : outer radius, default=1e5 Rs

Output:
- lognu : array with log10(nu)
- ll : array with log10(nu*Lnu)

v1 May 2012: ported from thindisk.nb, the Mathematica worksheet.
	"""
	# Definitions
	# =============
	# Constants 
	g=6.67428e-8	# G
	k=1.38065e-16	# Boltzmann constant
	c=29979245800.
	h=6.62607e-27	# Planck constant
	sigma=5.6704e-5	# Stefan-Boltzmann constant
	solarmass=1.99e33

	mass = m*solarmass
	rs = 2.*g*mass/c**2	# Schwarzschild radius
	rin = rtran*rs	# inner radius of thin disk
	rout = rout*rs	# outer radius
	angi = numpy.radians(angle)	# inclination angle of the disk

	# Eddington luminosity and accretion rate, defined assuming an efficiency of 10%
	ledd = 1.3e38*m	# erg/s
	mdotedd = ledd/(0.1*c**2)	# g/s
	mdotp = mdot*mdotedd	# Mdot in physical units

	# Thin disk spectrum
	# ====================
	def f(r,nu):
		# Computes the integrand in equation 5.45 of Frank et al. (the
		# accretion book) as a function of radius r and frequency nu
	
		# Local flux produced by viscous dissipation
		fvisc = 3./(8.*numpy.pi)*g*mass*mdotp/r**3*( 1.-numpy.sqrt(3.*rs/r) )

		# Disk effective temperature
		temp = (fvisc/sigma)**0.25

		# Expression that will be integrated over R
		return r/(numpy.exp( h*nu/(k*temp) ) - 1.)

	# Array of log10(nu)
	lognu=numpy.linspace(nui,nuf,n)
	# Array of log10(nu*Lnu)
	ll=numpy.zeros_like(lognu)

	i=0
	for y in lognu:
		nu=10**y
	
		# x below is the integration variable (radius)
		# Let's do a trick to make the integration work for lognu>~15.5
		if y>15.5: 
			integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout/10.)[0]
		else:
			integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout)[0]

		# This is basically equation 5.45 in Frank et al.
		lumdisk = 16.*numpy.pi**2*h*numpy.cos(angi)*nu**3/(c**2)*integral
		ll[i]=numpy.log10(nu*lumdisk)
	
		i=i+1
	
	if mfile!=None: 
		numpy.savetxt(mfile,numpy.transpose((lognu,ll)))
			
	return lognu,ll
	
	


def thindisk2(lognu,m,mdot,rtran,angle,rout=1e5):
	"""
Computes the spectrum of a standard Shakura-Sunyaev thin accretion disk.
The units of the code are CGS. Given the array lognu this method computed 
the output log(nu*Lnu) computed at the input frequencies for the given 
parameters.

Usage:
>>> lognu=numpy.linspace(13,15,50)
>>> ll=thindisk2(lognu,4e6,1e-3,3.,30.)

Arguments:
- m : black hole mass in solar masses
- mdot : mass accretion rate in Eddington units
- rtran : inner radius of the disk in Schwarzschild radii
- angle : inclination angle in degrees
- rout : outer radius, default=1e5 Rs

Output:
- ll : array with log10(nu*Lnu)

v1 May 2012: ported from thindisk.
	"""
	# Definitions
	# =============
	# Constants 
	g=6.67428e-8	# G
	k=1.38065e-16	# Boltzmann constant
	c=29979245800.
	h=6.62607e-27	# Planck constant
	sigma=5.6704e-5	# Stefan-Boltzmann constant
	solarmass=1.99e33

	mass = m*solarmass
	rs = 2.*g*mass/c**2	# Schwarzschild radius
	rin = rtran*rs	# inner radius of thin disk
	rout = rout*rs	# outer radius
	angi = numpy.radians(angle)	# inclination angle of the disk

	# Eddington luminosity and accretion rate, defined assuming an efficiency of 10%
	ledd = 1.3e38*m	# erg/s
	mdotedd = ledd/(0.1*c**2)	# g/s
	mdotp = mdot*mdotedd	# Mdot in physical units

	# Thin disk spectrum
	# ====================
	def f(r,nu):
		# Computes the integrand in equation 5.45 of Frank et al. (the
		# accretion book) as a function of radius r and frequency nu
	
		# Local flux produced by viscous dissipation
		fvisc = 3./(8.*numpy.pi)*g*mass*mdotp/r**3*( 1.-numpy.sqrt(3.*rs/r) )

		# Disk effective temperature
		temp = (fvisc/sigma)**0.25

		# Expression that will be integrated over R
		return r/(numpy.exp( h*nu/(k*temp) ) - 1.)

	# The if condition just below makes sure we can use this method even
	# if the user provides only one input frequency. A for loop can't 
	# iterate over an array with only one element.
	if numpy.size(lognu)>1:
		# Array of log10(nu*Lnu)
		ll=numpy.zeros_like(lognu)
	
		i=0
		for y in lognu:
			nu=10**y
	
			# x below is the integration variable (radius)
			# Let's do a trick to make the integration work for lognu>~15.5
			if y>15.5: 
				integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout/10.)[0]
			else:
				integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout)[0]

			# This is basically equation 5.45 in Frank et al.
			lumdisk = 16.*numpy.pi**2*h*numpy.cos(angi)*nu**3/(c**2)*integral
			ll[i]=numpy.log10(nu*lumdisk)
	
			i=i+1
	else:
		nu=10**lognu
		if lognu>15.5: 
			integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout/10.)[0]
		else:
			integral=scipy.integrate.quad( lambda x: f(x,nu), rin, rout)[0]
		lumdisk = 16.*numpy.pi**2*h*numpy.cos(angi)*nu**3/(c**2)*integral
		ll=numpy.log10(nu*lumdisk)
	
	#ll=delweird(ll)
			
	return ll






def thindiskfor(lognu,mass,mdot,rtr,i):
	"""
Auxiliary method that allows us to use the fortran routine thindisk.ssd 
-- compiled using f2py, see software/fortran/thin disk f2py/ --
in the same way as nemmen.thindisk2. The problem is that the fortran routine 
takes as input arrays with 200 values (since we cannot allocate dynamic arrays). 
Formats input data appropriately.
	"""
	import thindisk	# thindisk.so

	x=numpy.zeros(200)
	x[:lognu.size]=10**lognu

	y=thindisk.ssd(x,lognu.size,mass,mdot,rtr,i)

	return y[:lognu.size]



	
	


























	
	

	
