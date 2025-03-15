#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__


#include <armadillo>
#include <vector>
#include <Particle.hpp>



class PenningTrap
{
public:

	double B0;
	double V0;
	double d;
	std::vector<Particle> particles;
	bool part_inter;

	int N;
	double E_const;
	double ke;

	std::vector<double> t;
	std::vector<std::vector<arma::vec>> r;
	std::vector<std::vector<arma::vec>> v;

	double w_V;
	double f;


	PenningTrap(double B_field, double V_field, double charDim, std::vector<Particle> particleList, bool particle_interaction, double f_amplitude, double w_V_frequancy);


	std::vector<arma::vec> ElectricForce(std::vector<arma::vec> r, double t);


	std::vector<arma::vec> MagneticForce(std::vector<arma::vec> v);


	std::vector<arma::vec> ParticleForce(std::vector<arma::vec> r);


	std::vector<arma::vec> force(std::vector<arma::vec> r, std::vector<arma::vec> v, double t);


	void solveFE(double time_interval, int steps);


	void solveRK4(double time_interval, int steps);


	double particlesInsideTrap();

};


#endif