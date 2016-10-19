from AST1100SolarSystem import AST1100SolarSystem
seed=6392
mySolarSystem=AST1100SolarSystem(seed)
mass_star=mySolarSystem.starMass
radius_star=mySolarSystem.starRadius
no_of_planets=mySolarSystem.numberOfPlanets
temperature_planet=mySolarSystem.temperature

mass_launcher=1100
mass_lander=90
area_launcher=15
area_lander=6.2


print "-------------STAR-------------"
print "Mass of star: "+str(mass_star)
print "Radius star: "+ str(radius_star)
print "Temperature star: " + str(temperature_planet)
print "------------------------------"
print "Number of planets: " + str(no_of_planets)
print "------------------------------"


major_axis=mySolarSystem.a
radius_planets=mySolarSystem.radius
mass_planets=mySolarSystem.mass
period_of_planets=mySolarSystem.period
atmosphere_density=mySolarSystem.rho0


for i in range(no_of_planets):
	print "----------Planet " + str(i) + "--------"
	print "Major axis: " + str(major_axis[i])
	print "Radius of planet: "+ str(radius_planets[i])
	print "Mass of planet: " + str(mass_planets[i])
	print "Density of atmosphere: "+str(atmosphere_density[i])
print "--------------------------------"
