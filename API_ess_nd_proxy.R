# Script for the evalutation of relative steady-state error in API
# NEGATIVE DISTURBANCE CASE
# PROXY DESIGN WITH DIFFERENTIAL DEGRADATION

# Load packages
library (lattice)
library(animation)

# Clean up & initialize
rm (list = ls ())
graphics.off ()

# Fixed parameters
d = log(2)/25
w1 = 5*d
w2 = d
k5 = 0;
k3 = 0.03;

# Define number of points of grid in each dimension
N = c (100, 100, 100)

# Create grid of log-spaced values for variable parameters
p1seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[1]))
p2seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[2]))
p3seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[3]))

# Preallocate array to store the essr and the ssy0
essr = array (0, N)
essrp = array (0, N)
ssy0 = array (0, N)

for (i in seq (1, N[1]))
{
	for (j in seq (1, N[2]))
	{
		for (k in seq (1, N[3]))
		{
			# Set sequence values to current parameter values
			k1 = p1seq[i]
			k2 = p2seq[j]
			k4 = p3seq[k]
			
			# Analytical steady state of y from Mathematica
			ssyw = (2*d*k1*k4)/((d+k5+w1)*(d^2-k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5+w1))))
			
			# Analytical steady state of y from Mathematica without disturbance
			ssynw = (2*d*k1*k4)/((d+k5)*(d^2-k1*k3+sqrt(d^4+(k1^2*k3^2)+2*d*k1*k3*(d+(2*k2*k4)/(d+k5)))))
			
			# Operon effect
			k4 = 0.75 * k4
			
			# Analytical steady state of p from Mathematica
			sspw = (2*d*k1*k4)/((d^2-k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5+w1)))*(d+k5+w2))
			
			# Analytical steady state of p from Mathematica withtout disturbance
			sspnw = (2*d*k1*k4)/((d+k5)*(d^2-k1*k3+sqrt(d^4+(k1^2*k3^2)+2*d*k1*k3*(d+(2*k2*k4)/(d+k5)))))
			
			essr[i,j,k] = (ssynw-ssyw)/ssynw
			ssy0[i,j,k] = ssynw
			essrp[i,j,k] = (sspnw-sspw)/sspnw
		}
	}
}

# Display the percentage of grid points with rel error less than 5%
cat ('\n\nPercentage of grid points with less than 5% relative error:\n')
cat (sum(essr<0.05)/length(essr))
cat ('\n')


# Function to create animation
anifun = function (s) 
{
	# Init graphic device
	quartz (width = 7, height = 5)
	par (mar = c (4, 4, 0, 0) + 0.2)
	# Set axes labels
	xl = expression ('k'[1] ~ 'log nM'^-1 ~ 'min'^-1)
	yl = expression ('k'[2] ~ 'log min'^-1)
	# Call empty plot
	image (log10(p1seq), log10(p2seq), essr[,,s], col = 'white', xlab = xl, ylab = yl, main = NA)
	# Add region of essr < 0.10
	.filled.contour(log10(p1seq), log10(p2seq), essrp[,,s], levels = c(-1, 0.05), col = c('khaki2', 'white'))
	# Add region of essr < 0.05
	.filled.contour(log10(p1seq), log10(p2seq), essr[,,s], levels = c(-1, 0.05), col = c('lightgreen', 'white'))
	# Calculate levels for contours of essrp
	levs = c(-1, seq (min(essrp[,,s]), max(essrp[,,s]), length.out = 20))
	levs = signif (levs, 2)
	# Add contour plots of essrp
	contour (log10(p1seq), log10(p2seq), essrp[,,s], levels = levs, labcex = 1, add=T)
}

10^(seq( log10 (min(ssy0[,,10])), log10(max(ssy0[,,10])), length.out = 10))

# Create the animation
# setwd ('~/Desktop')
# saveGIF (for (s in 1:N[3]) anifun(s), movie.name = 'API_essr_nd.gif') 
