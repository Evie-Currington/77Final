Project outline:
Add several planets into a 2D plane with the initial positions, masses, and starting velocities.
Next, iterate through each planet and calculate its future position by using its velocity and the effect of gravity pulling it towards the other planets.
We will calculate the effect of gravity using the Newtonian formula, not gravitational relativity.
Plot a graph of the planets with their positions over time

Division of labor:
Evie: Calculating gravitational force and future locations
Federico: Interactivity with the audience so that the user can change the position, mass, and velocity of the planets
Alec: Visualization of the planets as a graph or a video

Tools/packages needed:
Numpy
Matplotlib
Github
Vpython

Potential functions needed:
def find_distance(a, b):
    	# return distance between two points
    	return np.sqrt((a[0]-b[0]) ** 2 + (a[1] - b[1]) ** 2)

Def find_gravity(a, b):
	#return gravitational pull between two points

def find_bearing(a, b):	
	#return angle between two points

Skeleton code:

Github repository:
https://github.com/thompalec/physics88final

Will post project description as well as versions as we finish coding.  



Applications:
Two-body problem
Three-body problem
Solar system simulation

pluto_d = 5.9*10**9
    
    
