'''
ASTR630: Class Project
Dynamical Relaxation in Crowded Systems
Author of code: Salama Algaz
Date: May 8, 2023
Last Updated: Septemeber 30,2023

This code is part of the extention on the paper written by Adams F.C. and Laughlin, Gregory.
In the paper, simulations were done on systems of 10 giant planets randomly distributed
between 5 Au to 30 Au distance from the sun. The scattering of the planets was examined 
and the resulting eccentricity distribution was compared to observation. The positions
of each of the planets and their trajectories are plotted as well.

This code uses the HNbody integrator that was developed by Kevin P. Rauch  and Douglas P. Hamilton (2002).
It writes the options of the integration scheme into a file and then runs the hnbody code on the file.
The positions and orbital elements of the planets are then obtained and plotted. Specifically this code will
run simulations similar to that in the paper and compare it to its results.

'''

################## Packages ######################

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import random


################## Functions #####################

def run_hnbody():
	'''
	Run the HNbody.exe file created by Rauch & Hamilton (2002) with the initial
	conditions of the solar system being specified in the file inputs.hnb 
	(The file inputs.hnb is written in the write_hnbfile function)
	The HNbody code then outputs a data file for each planet containing
	its orbital elements at each time step.
	'''
	from subprocess import call
	call(["./hnbody", 'inputs.hnb'])

def write_hnbfile(N, inits, tf, ts, to):
	'''
	Arguments:
	N:	 Number of planets 			(int)
	inits:	 Inital properties of planets		(str)
	tf: 	 Final time value 			(float) [years]
	ts:      Time step	 			(float) [years]
	to:      Output time step (for plots) 		(float) [years]

	Writes an hnb file to give to the hnbody code.
	'''

	# Open file
	f = open("inputs.hnb", "w")

	# Convert time variables to strings
	tf = str(tf); ts = str(ts); to = str(to)

	# Create text (str) for the file that will be written
	text = []

	# Integration options
	text.append('# Integration Options \n')		
	text.append('Integrator = Runge-Kutta \n')
	text.append('Corrector = True \n')
	text.append('\n\n')

	# Properties of initial conditions
	text.append('# Properties of initial conditions \n')	
	text.append('AngleUnit:  deg \n')			# base unit choices: deg, rad
	text.append('LengthUnit: AU \n')			# base unit choices: m, AU (or ua), pc
	text.append('MassUnit:   Msun \n')			# base unit choices: g, Msun
	text.append('TimeUnit:   yr \n')			# base unit choices: s, h, d, yr
	text.append('StepSize:  '+ ts +' \n')			# integration time step       (units of TimeUnit)
	text.append('M  =  1 \n')				# mass of the dominant object (units of MassUnit)
	text.append('N  = '+ str(N+1) +' \n')			# number of objects (including the dominant mass)				
	text.append('\n\n')

	# Input order
	text.append('# Input Order \n')	
	text.append('InputOrder: Mass SemiMajorAxis  Eccentricity Inclination '+ \
		'LongAscendNode ArgPeriapse  MeanAnomaly \n')
	text.append('\n\n')

	# Append initial orbital elements of planets to the text
	# If its a string simply append to the text
	if type(inits) == str:
		text.append(inits + ' \n')
	# Otherwise its a list of strings so go through all elements
	else:	
		for init in inits:
			text.append(init + ' \n')
	text.append('\n\n')

	# Outuput options
	text.append('# Output Options \n')
	text.append('Tfinal = ' + tf + '\n')
	text.append('\n\n')

	text.append('OutputInterval = ' + to + '\n')
	text.append('OutputFiles    = planet%d.dat \n')
	text.append('OutputOrder    = time x1 x2 x3 semi ecc incl longasc argperi mass \n')
	text.append('OutputCoord    = Bodycentric \n')
	text.append('OutputDigits   = 16 \n')
	text.append('OutputHeader = False \n')
	text.append('\n\n')

	# Write text to file and close file
	f.writelines(text)
	f.close()

def create_inits(N, a, mlog):
	'''
	N: 	Number of planets					(int)
	a: 	list of 2 elements for limits of semi-major axes 	(list) [AU]
	mlog:	Whether the masses are distributed logarithmically 	(boolean)

	Create a list whose elements are strings.
	Each string contains the initial conditions of each planet
	
	Returns:
	inits:  list of initial conditions for each planet 		(list)
	'''

	# Create list of initial conditions
	inits = []
	# Create list of the logarithm of semi-major axis (the semi-major axis is picked randomly logarithmically)
	alogs = np.linspace(np.log10(a[0]), np.log10(a[1]), 100001)
	# Create list of possible arguments of periapsis (where the planet is at its orbit)
	w = np.linspace(0, 360, 100001)

	# Check if masses should be distrubted logarithmically or not
	if mlog:
		masses = np.linspace(-1, 1, 100001)
	# Otherwise distribute them uniformly
	else:
		masses = np.linspace(0, 0.004, 100001)

	# For each planet
	for i in range(N):
		# Pick a random semi-major axis
		ai = 10**random.choice(alogs)
		# Pick a random argument of periapsis (in degrees)
		wi = random.choice(w)

		# Pick a random mass
		if mlog:
			mi = (10**random.choice(masses))*0.001
		else:
			mi = random.choice(masses)

		# Create string of initial conditions
		string = str(mi) + '  ' + str(ai) +'  0.0    0.0     0.0  '+ str(wi) + '  0.0'
		# Add to list
		inits.append(string)

	# Return list of initial conditions
	return inits

def get_data_of_planets(N, tf, to):
	'''
	N: 	 Number of planets 		 	(int)
	tf: 	 Final time value 			(float) [years]
	to:      Output time step (for plots) 		(float) [years]
	
	This function takes in the number of planets and reads in the orbital elements of 
	each planet into its own matrix.
	The columns of the matrix are:
		time | x1 | x2 | x3 | semi | ecc | incl | longasc | argperi | mass
	where x1, x2, and x3 are the positions of the planet relative to the ecliptic plane
	The semi stands for semi-major axis, the ecc for eccentricity of orbit,
	the incl for the inclination of the orbit, the lonasc for the long ascending node
	and the argperi stands for the argument of periapsis.
	
	The number of columns of the matrix is 10 while the number of rows is determined by
	the final time value and the time step: tf/to
	
	Returns:
	list_of_datas:  list of data matrices for each planet. The data matrix contains all the  (list)
	orbital properties of the planet at each time interval.
	'''

	# List that holds the data matrices for each planet
	list_of_datas = []

	# Number of lines in data files
	num_lines = int(tf/to)

	# Go through every planet
	for i in range(N + 1):
		
		# Open data file corresponding to planet
		file = open('planet{}.dat'.format(i), 'r');
		
		# Create a temporary matrix
		temp_data = np.zeros((num_lines + 1, 10), dtype=np.ndarray)
		
		# Go through the data file and read in the numbers line by line
		for i in range(num_lines + 1):
			# Read in line and split it to a list of ten elements
			line = file.readline().split()
			# Convert the line to float then append to the matrix
			temp_data[i] = np.array(line, dtype=np.float32)
		
		# Close file
		file.close()
		# Add data matrix to the list
		list_of_datas.append(temp_data)

	# Return list of data matrices. The length of the list is equal to the number of planets
	return list_of_datas

def get_data_of_remaining_planets(list_of_datas, avg_time, tf):
	'''
	Arguments:
	list_of_datas: 	list of data matrices for each planet  		(list)
	avg_time:	Time to average over the eccentricities		(float) [years]
	tf:		Final time  					(float) [years]

	This function gets the number of the planets ejected by looking at the eccentricity of 
	the planet over the avg_time towards the end of simulation and checks to see if the eccentricity
	is larger than 1 which means the planet is unbound and has left the planetary system.
	The average is taken because the planet's eccentricity can go over 1 but return to normal values
	so an average of the eccentricity is taken to determine if the planet's eccentricity is consistantly over 1 
	
	Returns:	
	num_planets_ej:		Number of planets that have left the solar system 		(int) 	
	mean_masses:		Average mass of the sysem before and after ejection  		(list 1x2) [mj: jupiter's mass]
	major_axes:		Semi-major axes of planets in system before and after 		(list 1x2) [AU]
	eccs:			Eccentricties of planets in system at final time		(list)

	'''

	# Get the time where averaging the eccentricities begins 
	begin_time = tf - avg_time
	# Number of planets ejected from the system
	num_planets_ej = 0
	# List of the masses of the planets left
	masses_left = []
	# List of masses of all planets
	masses = []
	# List of semi-major axes (initial)
	major_axes_init = []
	# List of semi_major axis (final) 
	major_axes_left = []
	# List of eccentricities (final, initial is just zero by choice)
	eccs = []

	# Go over each body (excluding star)
	for body in list_of_datas[1:]:
		
		# Get mass of body
		masses.append(body[0, 9])
		
		# Get semi-major axis of body
		major_axes_init.append(body[0, 4])

		# Get the index where begin_time for averaging occurs
		ind, = np.where(body[:,0] == begin_time)
		
		# Get the eccentricity values from begin_time to final_time
		e = body[int(ind):,5]
		# Get the average of it
		e_avg = np.mean(e)

		# If the average eccentricity is greater than 1 it means the planet has been ejected
		if e_avg > 0.95:
			num_planets_ej += 1
		# If not then the planet is still in the system, so get its mass, semi-major axis and eccentricity
		else:
			masses_left.append(body[0,9])
			eccs.append(body[-1, 5])
			major_axes_left.append(abs(body[-1, 4]))

	# Get average of masses of planets in system before planets were ejected and after
	mean_masses = [np.mean(masses)*1e3, np.mean(masses_left)*1e3]	# In terms of Jupiter's mass
	# Get major axes before and after
	major_axes = [major_axes_init, major_axes_left]

	# Return values
	return num_planets_ej, mean_masses, major_axes, eccs

def plot_system(list_of_datas, save = False):
	'''
	list_of_datas:	list of the matrices containing orbital elements of planets  	(list)
 	save:		whether to save a png file of the plot				(boolean)
	
	Outputs the orbits and trajectories of the system
	'''
	# Counter for for-loop
	count = 1

	# Plot the orbit of each body (excluding the sun)
	for body in list_of_datas[1:]:
		x = body[:,1]		# x-position
		y = body[:,2]		# y-position
		a = body[:,4]		# semi-major axis
		e = body[:,5]		# eccentricity

		# If its a circular orbit
		if round(e[0], 1) == 0:
			xe = np.linspace(-a[0], a[0], 1001)
			ye = np.sqrt(a[0]**2 - xe**2)

		# Else its an elliptic orbit
		else:
			b = a[0]*np.sqrt(1 - e[0]**2)			# Minor-axis
			c = np.sqrt(a[0]**2 - b**2)			# Distance of focus
			xe = np.linspace(-(c + a[0]), a[0] - c, 1001)
			ye = b*np.sqrt(1 - (xe + c)**2/a[0]**2)

		# Plot starting ellipses
		# plt.plot(xe, ye, '--')
		# plt.plot(xe, -ye, '--')

		# Plot planet's positions
		plt.plot(x, y, '-',markersize = 2, label='_nolegend_',zorder=2)
		# Plot planet's initial position
		plt.plot(x[0], y[0], 'k*',markersize = 5,zorder=3)
		plt.annotate(str(count), (x[0], y[0]))			# Label
		# Plot planet's final position
		plt.plot(x[-1], y[-1], 'k^',markersize = 5,zorder=3)
		plt.annotate(str(count), (x[-1], y[-1]))		# Label
		# Specifiy starting and ending positions
		plt.legend(['Start Position', 'End Position'])
		
		count += 1

	# Get final time
	tf = list_of_datas[0][-1,0]

	# Plot star
	plt.plot(0, 0, 'ko')	
	plt.grid(zorder=1)
	plt.title('Orbits of bodies after {0:.0f} years'.format(tf))
	plt.xlabel('AU'); plt.ylabel('AU')
	plt.xlim([-35, 35]); plt.ylim([-35,35])
	
	if save:
		plt.savefig('Plot_Of_System.png', bbox_inches='tight')
	
	plt.show()
	plt.clf()
	plt.close()

def animate_system(list_of_datas, save = False):
	'''
	Arguments:
	list_of_datas:	list of the matrices containing orbital elements of planets  	(list)
 	save:		whether to save a png file of the plot				(boolean)

	Creates an animation of the planets in motion around the sun and saves it as a gif (optional)
	'''

	# Create figure
	plt.style.use('dark_background')

	fig, ax = plt.subplots()
	#plt.grid()

	# Size of plot in AU
	s = 31
	ax.set(xlim = [-s, s], ylim = [-s, s])
	ax.set_xlabel('Distance [AU]'); ax.set_ylabel('Distance [AU]'); 
	# Plot sun in the middle
	plt.plot(0,0,'wo',markersize = 8,label='_nolegend_')

	# List of x-positions for all planets (row 1: planet 1 and so on, while the columns are the x-values)
	xs = []
	# List of y-positions
	ys = []
	# List of planets
	planets = []
	# time
	t = list_of_datas[0][:,0]
	# num of planets
	N = len(list_of_datas[1:])
	# legend string for plot (labelling all planets and their masses)
	leg = []
	# counter for for-loop
	count = 1
	# Colors for planets
	colors = plt.cm.hsv(np.linspace(0.2,1,N + 0))
	# How many points to show on frame
	pts = 1

	# Go through each planet
	for planet in list_of_datas[1:]:
		
		# Get the semi-major axis
		a = planet[0,4]

		# Get starting orbit
		xe = np.linspace(-a, a, 1001)
		ye = np.sqrt(a**2 - xe**2)

		# Plot the starting orbit
		plt.plot(xe, ye, 'k--', lw=1,label='_nolegend_', color = colors[count - 1])
		plt.plot(xe, -ye, 'k--', lw=1,label='_nolegend_', color = colors[count - 1])

		# Get x-positions and y-positions
		xs.append(planet[:,1])
		ys.append(planet[:,2])

		# Create scatter for planet
		planets.append(ax.scatter(0,0, color = colors[count - 1]))

		# Get mass of planet in terms of jupiters mass
		m = planet[0,9]*1e3

		# Add string to leg
		leg.append('P{0}, M = {1:.2f} mj'.format(count, m))

		# Add count
		count += 1


	# convert from list to nparrays
	xs = np.array(xs)
	ys = np.array(ys)

	# Add legend to plot
	plt.legend(leg, loc= 'upper left')

	def update(frame):
		# for each frame, update the data stored on each artist.
		for i in range(len(planets)):
			data = np.stack([xs[i,frame:frame + pts], ys[i,frame:frame + pts]]).T
			planets[i].set_offsets(data)

		# Update title of plot
		plt.title('Solar System with {0:.0f} Planets \nTime = {1:.0f} years'.format(N, t[frame+2]))
		
		return planets

	# create animation and display it
	MAX_frames = int(t[-1]) - 3 	# Displays the entire evolution of the system
	num_frames = 400		# Custom number of frames to display (each frame is to (time output) years)
	ani = animation.FuncAnimation(fig=fig, func=update, frames= num_frames, interval=10)
	if save:
		ani.save('Solar_System.gif', writer = 'ffmpeg', fps = 20)
	plt.show()
	
def plot_planets_orbits(list_of_datas, save = False, initial_time = 0):
	'''
	list_of_datas:	list of the data matrices for each planet 	(list)
 	save:		whether to save a png file of the plot		(boolean)
	initial_time:	The initial time of plots (default is 0)	(float)[years]
	
	Creates a plot with 4 subplots and saves it for each planet.
	The subplots are the semi-major axis, eccentricity, inclination, 
	and orbit (x, y position) of the planet as a function of time
	'''
	
	# Counter for the for-loop
	count = 1

	# Get final time
	tf = list_of_datas[0][-1,0]
	
	# Plot the elements of each body (excluding star)
	for body in list_of_datas[1:]:
		
		# Get index where time is equal to start time
		ind, = np.where(body[:,0] == initial_time)
		ind = int(ind)

		# Get orbitals
		t = body[ind:, 0]	# time
		x = body[ind:, 1]	# x-positions
		y = body[ind:, 2]	# y-positions
		a = body[ind:, 4]	# Semi-major axes
		e = body[ind:, 5]	# Eccentricities
		i = body[ind:, 6]	# Inclinations
		mass = body[0, 9]	# Mass of body
		ai = body[0, 4] 	# Initial semi-major axis of body

		# Do a subplot
		plt.figure()
		plt.subplot(2, 2, 1); 
		plt.suptitle('Planet {0}, Mass = {1:.2f}mj, Initial SemiMajorAxis = {2:.2f} AU'.format(count, mass*1e3, ai) + \
			'\nTime = {0:.0f} yrs to {1:.0f} yrs'.format(initial_time, tf))

		# Plot semi-major axis
		plt.plot(t, a); plt.grid()
		plt.xlabel('Time [yr]'); plt.ylabel('Semi-major axis [AU]')
		
		# Plot eccentricity
		plt.subplot(2, 2, 2)
		plt.plot(t, e); plt.grid()
		plt.xlabel('Time [yr]'); plt.ylabel('Eccentricity')
		
		# Plot inclination
		plt.subplot(2, 2, 3)
		plt.plot(t, i); plt.grid()
		plt.xlabel('Time [yr]'); plt.ylabel('Inclination [deg]')
		
		# Get shape of orbit
		# If its a circle
		if round(e[0], 1) == 0:
			xe = np.linspace(-a[0], a[0], 1001)
			ye = np.sqrt(a[0]**2 - xe**2)

		# Else its an ellipse
		else:
			b = a[0]*np.sqrt(1 - e[0]**2)			# Minor-axis
			c = np.sqrt(a[0]**2 - b**2)			# Distance of focus
			xe = np.linspace(-(c + a[0]), a[0] - c, 1001)
			ye = b*np.sqrt(1 - (xe + c)**2/a[0]**2)

		# Plot orbit
		plt.subplot(2, 2, 4)
		# Plot initial orbit
		plt.plot(xe, ye, 'k--',label='_nolegend_')
		plt.plot(xe, -ye, 'k--',label='_nolegend_')
		# Plot positions
		plt.plot(x, y, '.',label='_nolegend_'); plt.grid()
		#Plot starting position
		plt.plot(x[0], y[0], 'm.')
		#Plot ending position
		plt.plot(x[-1], y[-1], 'r.')
		# Specifiy starting and ending positions
		plt.legend(['Start Position', 'End Position'])
		# Plot star
		plt.plot(0, 0, 'ko')	
		plt.xlabel('x [AU]'); plt.ylabel('y [AU]')
		plt.tight_layout()

		if save:
			plt.savefig('Planet{0}_Plot.png'.format(count))
		
		plt.show()
		plt.clf()
		plt.close()

		# Add to count
		count += 1	

def print_results(N, tf, num_planets_ej, mean_masses, major_axes, eccs):
	'''
	Arguments:
	N: 			number of planets						(int)
	tf:			Final time values						(float) [years]
	num_planets_ej:		Number of planets that have left the solar system 		(int) 	
	mean_masses:		Average mass of the system before and after ejection  		(list 1x2) [jupiter's mass]
	major_axes:		Semi-major axes of planets in system before and after 		(list 1x2) [AU]
	eccs:			Eccentricties of planets in system at final time 		(list)

	Creates two figures, displays and saves them as png files:
	The first figure contains two subplots. The first subplot displays a histogram
	of the semi-major axes of the planets in the initial setup. The second subplot
	displays the same histogram at the final time.

	The second figure contains one plot and a table. The plot displays a historam of 
	the eccentricities of the planets at the final time value. The table displays 
	the average initial and final mass of the system , the initial and final minimum 
	semi-major axis of the system, the average eccentricity of the final, and the number
	of planets ejected. 
	'''
	
	# Get number of simulations:
	num_sims = len(major_axes[0])/N

	# Plot figure 1 (two subplots)
	fig, axs = plt.subplots(2,1)

	# Plot initial distribution of semi-major axes
	axs[0].hist(major_axes[0],edgecolor = "black", bins=10)
	axs[0].set_title('Distribution of semi-major axes of planets at t = {0:.0f} yrs [AU]'.format(0))
	axs[0].set_xlabel('Semi-major axis [AU]')
	axs[0].set_ylabel('N (Number of planets)')
	
	# Plot final distribution of semi-major axes
	axs[1].set_title('Distribution of semi-major axes of planets at t = {0:.0f} yrs [AU]'.format(tf))
	axs[1].hist(major_axes[1],edgecolor = "black", bins=10)
	axs[1].set_xlabel('Semi-major axis [AU]')
	axs[1].set_ylabel('N (Number of planets)')
	
	plt.tight_layout()
	# plt.xticks(np.arange(1, 30, step=1))

	# Save figure and show it
	plt.savefig('Semi_major_axes_histograms.png')
	plt.show()

	# Plot figure 2 (1 plot and 1 table)
	fig, axs = plt.subplots(2,1)
	plt.title('Number of planets in system = {0:.0f} \nNumber of simulations = {1:.0f}'.format(N, num_sims))
	# Plot final distribution of eccentricities
	axs[0].hist(eccs,edgecolor = "black")
	axs[0].set_title('Distribution of e at t = {0:.0f} yrs'.format(tf))
	axs[0].set_xlabel('Eccentricity')
	axs[0].set_ylabel('N (Number of planets)')
	
	# Plot table
	collabel = ('t = {0:.0f} yrs'.format(0),'t = {0:.0f} yrs'.format(tf))
	rowlabel = ('Average mass of the system [mj]','Minimum semi-major axis [AU]',
		'Average eccentricity','Average num of planets ejected')
	table_data = [[round(mean_masses[0], 2), round(mean_masses[1], 2)],
			[round(min(major_axes[0]), 2), round(min(major_axes[1]), 2)],
			[0, round(np.mean(eccs), 2)],
			['N/A', num_planets_ej]]

	axs[1].axis('off')
	axs[1].axis('tight')
	axs[1].table(cellText=table_data,colLabels=collabel,rowLabels=rowlabel,loc='center')
	plt.tight_layout()

	# Save figure and show it
	plt.savefig('Eccentricity_and_table.png')
	plt.show()

def run_simulation(N, tf, ts, to, mlog, plot_orbits = False, animate = False):
	'''
	Arguments:
	N:  		 number of planets 					(int)
	tf: 		 Final time value 					(float) [years]
	ts:      	 Time step	 					(float) [years]
	to:     	 Output time step 					(float) [years]
	mlog:		 Whether masses are distrubted logarithmically 		(boolean)
	plot_orbits: 	 Save plots of planets if True  			(boolean)

	Run one simulation of the solar system for tf years with the following description:

	N planets with uninclined circular orbits are placed in a system with a star having the same size as the sun. 
	The semi-major axes of the planets are randomly picked from 5AU < r < 30AU (logarithmically distributed). 
	The masses are either randomly picked from a uniform distribution of 0 < mp < 4mj (mp is the planet's mass)
	(where mj is Jupiter's mass which in this case is taken to be equal to 0.001M where M is the sun's mass) 
	or from -1 < log(mp/mj) < 1 depending on whether the variable mlog is passed as True or False.

	Returns:
	num_planets_ej:	Number of planets that have left the solar system 	(int) 	
	mean_masses:	Average mass of the sysem before and after ejection  	(list 1x2)[mj]
	major_axes:	Semi-major axes of planets in system before and after 	(list 1x2)[AU]
	eccs:		Eccentricties of planets in system at final time 	(list)
	'''

	# Create initial conditions for N planets within 5AU < r < 30AU
	inits = create_inits(N, [5, 30], mlog)

	# Write file with initial conditions
	write_hnbfile(N, inits, tf, ts, to)

	# Run Nbody integrator
	run_hnbody()

	# Get data of planets
	list_of_datas = get_data_of_planets(N, tf, to)

	# Time duration to average the eccentricity (last 200 time steps):
	time_avg = 200*to

	# Get number of planets ejected and the mean masses,axes, and ecc with an averaging time of 10 yrs
	num_planets_ej, mean_mass, axes, eccs = get_data_of_remaining_planets(list_of_datas, time_avg, tf)

	# If plot_orbits is true then create plots of planets and system
	if plot_orbits:
		plot_system(list_of_datas)
		plot_planets_orbits(list_of_datas)
		
		# Print the results of the simulation
		print_results(N, tf, num_planets_ej, mean_mass, axes, eccs)

	if animate:
		animate_system(list_of_datas)
	
	# Returns results of simulation
	return num_planets_ej, mean_mass, axes, eccs

def run_Nsimulations(num_sims, N, tf, ts, to, mlog):
	'''
	Arguments:
	num_sims:	 number of simulations 	 			  	(int)
	N:  		 number of planets 					(int)
	tf: 		 Final time value 					(float) [years]
	ts:      	 Time step	 					(float) [years]
	to:     	 Output time step 					(float) [years]
	mlog:		 Whether masses are distrubted logarithmically 		(boolean)
	plot_orbits: 	 Save plots of planets if True  			(boolean)

	This function runs the simulation N times. The setup of the simulation is described
	in the function run_simulation.
	Prints the final results of the simulations combined.
	'''

	# Create lists containing results for each simulation:

	# List of number of planets ejected
	lst_nums_planets_ejected = []
	# List of mass averages mass (col1: initial, col2: final)
	lst_mean_masses = []
	# List of semi-major axes (col1: initial, col2: final)
	lst_major_axes = [[],[]]
	# List of eccentricities (final, initial is just zero by choice)
	lst_eccs = []

	# For each simulation
	for i in range(num_sims):
		
		# Print simulation number
		print('\n\n\n---------------------- Starting simulation {} --------------------\n\n\n'.format(i + 1))

		# Run simulation and get results
		num_planets_ej, mean_masses, axes, eccs = run_simulation(N, tf, ts, to, mlog)

		# Add to lists
		lst_nums_planets_ejected.append(num_planets_ej)
		lst_mean_masses.append(mean_masses)
		lst_major_axes[0].extend(axes[0])
		lst_major_axes[1].extend(axes[1])
		lst_eccs.extend(eccs)

	# Get the average of the mean masses:
	avg_mean_masses = [np.mean([row[0] for row in lst_mean_masses]),
						np.mean([row[1] for row in lst_mean_masses])]
	# Get the average of the number of planets ejected
	avg_num_planets_ej = np.mean(num_planets_ej)

	# Print results:
	print_results(N, tf, avg_num_planets_ej, avg_mean_masses, lst_major_axes, lst_eccs)


#################### Main ########################

# Variables
tf = 10000 		# Final time 				(float) [years]
to = 1  		# Output interval 			(float) [years]
ts = to		 	# Time step	 			(float) [years]
mlog = 1		# Log distribution of mass or not 	(boolean)
show_orbits = 0  	# Show orbits or not 			(boolean)
num_sims = 2 	 	# Number of simulations			(int)
N = 4 			# Number of planets			(int)
animate = 1 		# Animate or not 			(boolean)

# Run simulation(s)
run_simulation(N, tf, ts, to, mlog, show_orbits, animate)
# run_Nsimulations(num_sims, N, tf, ts, to, mlog)
