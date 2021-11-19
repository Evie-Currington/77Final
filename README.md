# Physics 77 Final


#position of the Earth
position1 = [1,2]
#position of the sun
position2 = [2,1.5*10**8]
#mass of the Earth
Earth_mass = 5.9 * 10**24 
#mass of the sun
sun_mass = 1.989 * 10**30
#Earth's average orbital velocity, unit:Km/s
Earth_orbit_velocity = 30


def force_gravity(planet1,planet2):
    
    #gravitational constant, units: m^3 kg^-1 s^-2
    G = 6.67408 * 10**(-11)   
    #vector difference
    distance = [(position1[0]-position2[0]), (position1[1]-position2[1])]
    #distance between mass 1 and mass 2
    distance_magnitude = np.sqrt(((position1[0]-position1[1])**2)+(position2[0]-position2[1])**2) 
    #unit vector
    r_hat = distance/distance_magnitude
    #force magnitude
    force_magnitude = (G*planet1*planet2)/distance_magnitude**2
    #force vector
    force_vector = -force_magnitude*r_hat
    print('The distance between the planet and the sun is {:.1e} km'.format(distance_magnitude))
    #print('The direction of the force between the planet and the sun is {:.1e}'.format(r_hat))
    print('The magnitude of the gravitational force between the planet and the sun is {:.1e}'.format(force_magnitude))
    #print('The gravitational force vector between the planet and the sun is {:.1e} km'.format(force_vector))
    
    return force_vector
    
force_gravity(Earth_mass,sun_mass)


distance_magnitude = np.sqrt(((position1[0]-position1[1])**2)+(position2[0]-position2[1])**2) 

def angular_momentum(planet_mass, planet_velocity, distance):
        
    #planet angular momentum with respect to the Sun
    L = planet_mass * planet_velocity * distance
    print('The angular momentum of the planet is {:.1e} kg*m^2/s'.format(L))
    
angular_momentum(Earth_mass, Earth_orbit_velocity, distance_magnitude)

#planet masses, unit:kg
Mercury_m = 3.3*10**23
earth_m = 5.9*10**24
saturn_m = 5.7*10**26
mars_m = 6.4*10**23
venus_m = 4.8*10**24
neptune_m = 10**23
uranus_m = 8.6*10**25
jupiter_m = 1.8*10**25
pluto_m = 1.2*10**22

#distances from sun, unit:km
mercury_d = 5.7*10**7
earth_d = 1.5*10**8
venus_d = 10**8
mars_d = 2.2*10**8
jupiter_d = 7.7*10**8
saturn_d = 1.4*10**9
uranus_d = 2.8*10**9
neptune_d = 4.4*10**9
pluto_d = 5.9*10**9
    
    
