# Script for the evalutation of relative steady-state error in API
# NEGATIVE DISTURBANCE CASE WITH SATURATION

# Load packages
library (deSolve)
library (animation)

# Clean up & initialize
rm (list = ls ())
graphics.off ()

# Fixed parameters
#da = c (log(2)/25, log(2)/120, log(2)/1440)
da = c (log(2)/25)
w = 5*da[1]
k5 = 0;
k3 = 100;

# Define accuracy threshold
thr = 1e-6
# Define max final time
Tmax = 50000
# Keep track of maximum simulation error
maxerr = 0

# Define number of points of grid in each dimesion
N = c (100, 100, 1)

# Create grid of log-spaced values for variable parameters
p1seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[1]))
p2seq = 10^(seq (log10 (1e-2), log10 (1e2), length.out = N[2]))
p3seq = 10^(seq (log10 (1), log10 (1e2), length.out = N[3]))

# Preallocate array to store the essr and the ssy0
essr = array (0, c(N,3))
ssy0 = array (0, c(N,3))

# Define function for the API_nd model
APInd <- function (t, y, p) {
	Vmax = 100
	k = Vmax/(2*p[4])
	list(c(
		p[1] - p[3]*y[1]*y[2] - p[6]*y[1],
		p[2]*y[3] - p[3]*y[1]*y[2] - p[6]*y[2],
		Vmax*y[1]/(k+y[1]) - (p[5] + p[6] + p[7])*y[3]
	))
}

# Main loop
for (i in seq (1, N[1]))
{
	cat ('\n')
	cat (i)
	
	for (j in seq (1, N[2]))
	{
		for (k in seq (1, N[3]))
		{
			# Set sequence values to current parameter values
			k1 = p1seq[i]
			k2 = p2seq[j]
			k4 = p3seq[k]
			
			for (z in seq (1,length(da)))
			{
				# Set d to the current 
				d = da[z]
	
				# Solve with disturbance
				TI = 0; TF = 100
				sol = ode (y = c (0,0,0), func = APInd, times = c (TI, TF), parms = c(k1,k2,k3,k4,k5,d,w), method = 'bdf')
				while (norm (sol[2,2:4]-sol[1,2:4], type='2') > thr && TF < Tmax)
				{
					TI = TF; TF = TF + 100
					sol = ode (y = sol[2,2:4], func = APInd, times = c (TI, TF), parms = c(k1,k2,k3,k4,k5,d,w), method = 'bdf')
				}
				ss = sol[2,4]
				
				# Warn of potential quality problems
				if (TF == Tmax)
				{
					cat ('\n\nWarning: Tmax reached\n')
					cat ('Norm of difference:\n')
					cat (norm (sol[2,2:4]-sol[1,2:4], type='2'))
					cat ('\n\n')
					maxerr = max (maxerr, norm (sol[2,2:4]-sol[1,2:4], type='2'))
				}
			
				# Solve without disturbance
				TI = 0; TF = 100
				sol = ode (y = c (0,0,0), func = APInd, times = c (TI, TF), parms = c(k1,k2,k3,k4,k5,d,0), method = 'bdf')
				while (norm (sol[2,2:4]-sol[1,2:4], type='2') > thr && TF < Tmax)
				{
					TI = TF; TF = TF + 100
					sol = ode (y = sol[2,2:4], func = APInd, times = c (TI, TF), parms = c(k1,k2,k3,k4,k5,d,0), method = 'bdf')
				}
				ss0 = sol[2,4]
				
				# Warn of potential quality problems
				if (TF >= Tmax)
				{
					cat ('\n\nWarning: Tmax reached\n')
					cat ('Norm of difference:\n')
					cat (norm (sol[2,2:4]-sol[1,2:4], type='2'))
					cat ('\n\n')
					maxerr = max (maxerr, norm (sol[2,2:4]-sol[1,2:4], type='2'))
				}
							
				# Calculate essr and ssy0
				ssy0[i,j,k,z] = ss0;
				essr[i,j,k,z] = (ss0-ss)/ss0;
			}
		}
		
		cat (' ')
		cat (j)
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
	image (log10(p1seq), log10(p2seq), essr[,,s,1], col = 'white', xlab = xl, ylab = yl, main = NA)
	# Add region of essr < 0.05 for mammalian cells
	#.filled.contour(log10(p1seq), log10(p2seq), essr[,,s,3], levels = c(0, 0.05), col = c('khaki2', 'white'))
	# Add region of essr < 0.05 for yeast
	#.filled.contour(log10(p1seq), log10(p2seq), essr[,,s,2], levels = c(0, 0.05), col = c('lightblue', 'white'))
	# Add region of essr < 0.05 for bacteria
	.filled.contour(log10(p1seq), log10(p2seq), essr[,,s,1], levels = c(0, 0.05), col = c('lightgreen', 'white'))
	# Calculate levels for contours of ssy0
	levs = 10^(seq( log10 (min(ssy0[,,s,1])), log10(max(ssy0[,,s,1])), length.out = 20))
	levs = signif (levs, 3)
	# Add contour plots of ssy0
	contour(log10(p1seq), log10(p2seq), ssy0[,,s,1], levels = levs, labcex = 1, add=T)
}

# Create the animation
# setwd ('~/Desktop')
# saveGIF (for (s in 1:N[3]) anifun(s), movie.name = 'API_sser_nd.gif') 
