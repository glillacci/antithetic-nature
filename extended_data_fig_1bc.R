# A universal biomolecular integral feedback controller for robust perfect adaptation
# Stephanie K. Aoki, Gabriele Lillacci, Ankit Gupta, Armin Baumschlager, David Schweingruber & Mustafa Khammash

# Simulations presented in Extended Data Figures 1b and 1c

# Load packages
library (lattice)

# Clean up & initialize
rm (list = ls ())
graphics.off ()

# Fixed parameters
da = c (log(2)/25, log(2)/120, log(2)/1440)
w = 5*da[1]
k5 = 0.5*da[1]
k3 = 0.05

# Define number of points of grid in each dimesion
N = c (100, 100, 100)

# Create grid of log-spaced values for variable parameters
mu_seq = 10^(seq (log10 (1e-0), log10 (1e4), length.out = N[1]))
theta_seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[2]))
k_seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[3]))

## The following section will loop through the 3-dimensional parameter grid for mu, theta and k
## and evaluate the steady-state error at each point

# Preallocate array to store the essr and the ssy0
essr = array (0, c(N, 3))
ssy0 = array (0, c(N, 3))

for (i in seq (1, N[1]))
{
	for (j in seq (1, N[2]))
	{
		for (k in seq (1, N[3]))
		{
			# Set sequence values to current parameter values
			k1 = mu_seq[i]
			k2 = theta_seq[j]
			k4 = k_seq[k]
			
			for (z in seq (1,3))
			{
				# Set d to the current 
				d = da[z]
				
				# Analytical steady state WITH disturbance
				ssyw = (2*d*k1*k4)/((d+w+k5)*(d^2-k1*k3 + sqrt((d^2 + k1*k3)^2 + (4*d*k1*k2*k3*k4)/(d + w + k5))))
				# Analytical steady state WIHTOUT disturbance
				ssynw = (2*d*k1*k4)/((d+k5)*(d^2-k1*k3+sqrt(d^4+(k1^2*k3^2)+2*d*k1*k3*(d+(2*k2*k4)/(d+k5)))))
				
				# Calculate and store the essr
				essr[i,j,k,z] = (ssynw-ssyw)/ssynw
				# Store the disturbance-free steady state
				ssy0[i,j,k,z] = ssynw
			}
		}
	}
}

# Display the percentage of grid points with error less than 5%
cat ('\n\nPercentage of grid points with less than 5% relative error:\n')
cat (sum(essr<0.05)/length(essr))
cat ('\n')


## Plotting function. This function will generate a 2-d plot showing a slice of the 3D parameter grid,
## i.e. a 2D plot for a fixed value of k. The argument s is the index of k_seq containing the value of k to plot.

anifun = function (s) 
{
	# Init graphic device
	quartz (width = 7, height = 5)
	par (mar = c (4, 4, 0, 0) + 0.2)
	# Set axes labels
	xl = expression ('mu' ~ 'log nM'^-1 ~ 'min'^-1)
	yl = expression ('theta' ~ 'log min'^-1)
	# Call empty plot
	image (log10(mu_seq), log10(theta_seq), essr[,,s,1], col = 'white', xlab = xl, ylab = yl, main = NA)
	# Add region of essr < 0.05 for mammalian cells
	.filled.contour(log10(mu_seq), log10(theta_seq), essr[,,s,3], levels = c(0, 0.05), col = c('khaki2', 'white'))
	# Add region of essr < 0.05 for yeast
	.filled.contour(log10(mu_seq), log10(theta_seq), essr[,,s,2], levels = c(0, 0.05), col = c('lightblue', 'white'))
	# Add region of essr < 0.05 for bacteria
	.filled.contour(log10(mu_seq), log10(theta_seq), essr[,,s,1], levels = c(0, 0.05), col = c('lightgreen', 'white'))
	# Calculate levels for contours of ssy0
	levs = 10^(seq( log10 (min(ssy0[,,s,1])), log10(max(ssy0[,,s,1])), length.out = 20))
	levs = signif (levs, 3)
	# Add contour plots of ssy0
	contour(log10(mu_seq), log10(theta_seq), ssy0[,,s,1], levels = levs, labcex = 1, add=T)
}


## Plot the figures in the paper
## This generates extended data fig. 1b
anifun(15)
## This generates extended data fig. 1c
anifun(65)