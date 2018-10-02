# Script for the evalutation of relative steady-state error in API
# NEGATIVE DISTURBANCE CASE

# Load packages
library (lattice)
library(animation)

# Clean up & initialize
rm (list = ls ())
graphics.off ()

# Fixed parameters
d = log(2)/25
w = 5*d
k5 = 0;
k3 = 100;

# Define number of points of grid in each dimesion
N = c (100, 100, 100)

# Create grid of log-spaced values for variable parameters
p1seq = seq(0.01, 100, length.out = N[1])
p2seq = seq(0.01, 100, length.out = N[2])
p3seq = seq(0.01, 100, length.out = N[3])

# Analytical steady state of y (with disturbance) from Mathematica
# ssy = (2*d*k1*k4)/((d+w+k5)*(d^2-k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+w+k5))))

# Analytical steady state of y (WITHOUT disturbance) from Mathematica
# ssy0 = k4*(-d^2+k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5)))/(2*k3*(k2*k4+d*(d+k5)))

# Preallocate array to store the essr and the ssy0
essr = array (0, N)
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
			
			# Analytical relative steady-state error from Mathematica
			# essr[i,j,k] = (k3*(k2*k4+d*(d+k5))*((-d^2+k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5)))/(k3*(k2*k4+d*(d+k5)))-(4*d*k1)/((d+w+k5)*(d^2-k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+w+k5))))))/(-d^2+k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5)))
			
			# Analytical steady state (w=0) from Mathematica
			# ssy0[i,j,k] = k4*(-d^2+k1*k3+sqrt((d^2+k1*k3)^2+(4*d*k1*k2*k3*k4)/(d+k5)))/(2*k3*(k2*k4+d*(d+k5)))
			
			# Analytical steady state from Mathematica
			ssyw = (2*d*k1*k4)/((d+w+k5)*(d^2-k1*k3 + sqrt((d^2 + k1*k3)^2 + (4*d*k1*k2*k3*k4)/(d + w + k5))))
			# Analytical steady state from Mathematica without disturbance
			ssynw = (2*d*k1*k4)/((d+k5)*(d^2-k1*k3+sqrt(d^4+(k1^2*k3^2)+2*d*k1*k3*(d+(2*k2*k4)/(d+k5)))))
			
			essr[i,j,k] = (ssynw-ssyw)/ssynw
			ssy0[i,j,k] = ssynw
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
	xl = expression ('k'[1] ~ 'log min'^-1)
	yl = expression ('k'[2] ~ 'log nM'^-1 ~ 'min'^-1)
	# Call empty plot
	image (p1seq, p2seq, essr[,,s], col = 'white', xlab = xl, ylab = yl, main = NA)
	# Add region of essr < 0.10
	.filled.contour(p1seq, p2seq, essr[,,s], levels = c(0, 0.1), col = c('khaki2', 'white'))
	# Add region of essr < 0.05
	.filled.contour(p1seq, p2seq, essr[,,s], levels = c(0, 0.05), col = c('lightgreen', 'white'))
	# Calculate levels for contours of ssy0
	#levs = seq(min(ssy0[,,s]), max(ssy0[,,s]), length.out = 20)
	levs = c(0.125, 0.25, 0.5, 1, 2, 8)
	#levs = signif (levs, 3)
	# Add contour plots of ssy0
	contour(p1seq, p2seq, ssy0[,,s], levels = levs, labcex = 1, add=T)
}

10^(seq( log10 (min(ssy0[,,10])), log10(max(ssy0[,,10])), length.out = 10))

# Create the animation
# setwd ('~/Desktop')
# saveGIF (for (s in 1:N[3]) anifun(s), movie.name = 'API_sser_nd.gif') 
