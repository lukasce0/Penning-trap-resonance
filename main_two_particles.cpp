

#include <iostream> //prefere not using namespace when I work with two types of vectors
#include <armadillo>
#include <vector>

#include "Particle.hpp"
#include "PenningTrap.hpp"










int main() {
	
	//
	//Create particles
	//
	//Defining inital vectors
	arma::vec r1 = arma::vec(3);
	arma::vec v1 = arma::vec(3);
	arma::vec r2 = arma::vec(3);
	arma::vec v2 = arma::vec(3);

	//Assigning initial values
	//position particle 1:
	r1[0] = 20.; //in mico m
	r1[1] = 0.;
	r1[2] = 20.;
	//velocity particle 1:
	v1[0] = 0.; //in m/s
	v1[1] = 25.;
	v1[2] = 0.;

	//position particle 2:
	r2[0] = 25.;
	r2[1] = 25.;
	r2[2] = 0.;
	//velocity particle 2:
	v2[0] = 0.;
	v2[1] = 40.;
	v2[2] = 5.;

	//Creating two particles (instances of class Particle)
	double q = 1; //in e
	double m = 40; //in u
	Particle Ca1(q , m, r1, v1);
	Particle Ca2(q, m, r2, v2);
	
	//
	//Create penning trap
	// 
	//Specifying values of the penning trap
	double d = 500; //micro m
	double V0 = 2.41 * pow(10, 6); //u micro m^2/(micro s ^2 e)
	double B0 = 96.5; //u/(micro s e)
	bool particle_interaction = false; //do we want to include electric particle interactions?

	//Creating an array containing particles
	std::vector<Particle> particleList;
	particleList.push_back(Ca1);
	particleList.push_back(Ca2);
	
	//Creatign a penning trap (instance of PenningTrap)
	PenningTrap pen_trap(B0, V0, d, particleList, particle_interaction, 0., 0.2); //f = 0 means time independant electric potential

	//
	//Running the simulation
	//
	int N = particleList.size(); //number of particles
	double time_interval = 50.; //simulation lengdth in micro s
	int steps = 2000; //number of steps in the simulation
	pen_trap.solveFE(time_interval, steps); //solving the using either FE or RK4 (solveRK4 and solveFE) are valid methods.

	double inside_trap = pen_trap.particlesInsideTrap();
	std::cout << inside_trap;


	//
	//Making text file with the results in order to plot it using Python
	//
	std::string filename1 = "positions.txt"; //output file name
	std::string filename2 = "time.txt";
	std::string filename3 = "velocity.txt";
	std::ofstream ofile1; //create an "output file stream"
	std::ofstream ofile2;
	std::ofstream ofile3;
	ofile1.open(filename1); //open it/connect to our file name
	ofile2.open(filename2);
	ofile3.open(filename3);

	//Write to the file (columbs containing vectors of repsective particles) 
	for (int n = 0; n <= steps; n++) {
		//Time:
		ofile2 << pen_trap.t[n] << std::endl;

		//Positions:
		for (int i = 0; i < N; i++) {
			ofile1 << pen_trap.r[n][i][0] << " ";
			ofile3 << pen_trap.v[n][i][0] << " ";
		}
		ofile1 << std::endl;
		ofile3 << std::endl;
		for (int i = 0; i < N; i++) {
			ofile1 << pen_trap.r[n][i][1] << " ";
			ofile3 << pen_trap.v[n][i][1] << " ";
		}
		ofile1 << std::endl;
		ofile3 << std::endl;
		for (int i = 0; i < N; i++) {
			ofile1 << pen_trap.r[n][i][2] << " ";
			ofile3 << pen_trap.v[n][i][2] << " ";
		}
		ofile1 << std::endl << std::endl;
		ofile3 << std::endl;
	}

	ofile1.close(); //Closing the output file
	ofile2.close();
	ofile3.close();




	return 0;
}
