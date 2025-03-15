#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>


class Particle
{
public:

	double q;
	double m;
	arma::vec r;
	arma::vec v; 


	Particle(double charge, double mass, arma::vec position, arma::vec velocity);

};



#endif