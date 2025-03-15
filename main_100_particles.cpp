

#include <iostream> //prefere not using namespace when I work with two types of vectors
#include <armadillo>
#include <vector>
#include <string>

#include "Particle.hpp"
#include "PenningTrap.hpp"











int main() {

	//
	//Create particles
	//

	double N = 100.0; //number of particles


	//Particle proparties
	double q = 1; //in e
	double m = 40; //in u

	//Penning trap proparties
	double d = 500; //micro m
	double V0 = 2.41 * pow(10, 6); //u micro m^2/(micro s ^2 e)
	double B0 = 96.5; //u/(micro s e)
	bool particle_interaction = false; //do we want to include electric particle interactions?


	//Declaring inital vectors
	arma::vec r;
	arma::vec v;
	
	std::vector<Particle> particleList; //particle list

	arma::arma_rng::set_seed(123); //setting seed

	for (int i = 0; i < N; i++) {

		//Assigning initial values
		r = arma:: vec(3).randn() * 0.1 * d;  //random initial position
		v = arma::vec(3).randn() * 0.1 * d;  //random initial velocity

		Particle particle_i = Particle(q, m, r, v); //create particles

		particleList.push_back(particle_i); //add particles to the list
	}
	
	//
	//Create penning traps
	//	

	double time_interval = 500.; //simulation duration in micro s
	int steps = 10; //number of steps in the simulation
	double inside_trap; //describes number of particles inside the trap at the end of the simulation
	double particle_in_factor; //describes the fraction of particles inside the trap at the end of simulation
	

	arma::vec f_list = arma::vec("0.2 0.4 0.7"); //list of f values
	double f; //Declare f


	double start = 0.2; //first value of w_V in MHz, to zoom in 2.1 was used
	double end = 2.5; //last value of w_Z in MHz to zoom in 2.4 was used
	double steps_w_V = 120; //Number of different w_V, 120 ensures difference of less than 0.02, to zoom in 15 was used
	arma::vec w_V_list = arma::linspace(start, end, steps_w_V);
	double w_V; //declare w_V
	

	//Loop for f
	for (int i = 0; i < f_list.size(); i++) {
		f = f_list.at(i);

		//Making text file with the results in order to plot them using Python
		std::string filename = "resonance_frequancy" + std::to_string(f) + ".txt"; //choose file name
		std::ofstream ofile;//create an "output file stream"
		ofile.open(filename); //open it/connect to our file name


		//Loop for w_V
		for (int j = 0; j < steps_w_V; j++) {
			w_V = w_V_list.at(j);


			//Creatign a penning trap (instance of PenningTrap)
			PenningTrap pen_trap(B0, V0, d, particleList, particle_interaction, f, w_V);

			//
			//Running the simulation
			//
			pen_trap.solveRK4(time_interval, steps); //solving the using either FE or RK4 (solveRK4 and solveFE) are valid methods.


			inside_trap = pen_trap.particlesInsideTrap(); //finding the number of particles inside the trap
			particle_in_factor = inside_trap / N;

			ofile << w_V << "   " << particle_in_factor << "   " << std::endl; //write result into the file
		}

		ofile.close(); //closing the output file
	}

	return 0;
}