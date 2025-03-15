

#include <armadillo>
#include <vector>
#include <cmath> //for the cosinus in E-field time dependancy

#include "Particle.hpp"



	//Constructor
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity) {
		q = charge; //charge e
		m = mass; //mass u
		r = position; //position Cartesian micro meter
		v = velocity; //velocity Cartesian meter per second (micro meter per micro second)
	}
