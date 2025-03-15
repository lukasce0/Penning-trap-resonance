

#include <armadillo>
#include <vector>
#include <cmath> //for the cosinus in E-field time dependancy

#include "PenningTrap.hpp"
#include "Particle.hpp"






	//Constructor:
PenningTrap::PenningTrap(double B_field, double V_field, double charDim, std::vector<Particle> particleList, bool particle_interaction, double f_amplitude, double w_V_frequancy)
	{
		//definiting attributes:
		B0 = B_field; //magnetic flux density
		V0 = V_field; //electroc potential
		d = charDim; //characteristic dimension micro meter
		particles = particleList; //array of Particle instances
		
		f = f_amplitude; //amplitude of the time dependant electric field
		w_V = w_V_frequancy; //frequancy of the time dependant electric field

		N = particles.size(); //number of particles
		E_const = -V0 / (2 * pow(d, 2)); //common constant for all E components
		ke = 1.38935333 * pow(10, 5); //coulomb constant in u * micro m^3/(micro s^2 * e^2)

		part_inter = particle_interaction; //do we want electrical interactions between functions? true-yes false-no

	}





	//Force from external Electric field:
	std::vector<arma::vec> PenningTrap::ElectricForce(std::vector<arma::vec> r, double t) {
		std::vector<arma::vec> F(N); //list of three-vector that we will return
		arma::vec F_i = arma::vec(3);
		double Eq_const; //E_const * q
		double time_dep_term = 1 + f * cos(w_V * t); //time dependancy term. Set f to 0 for time independancy

		
		for (int i = 0; i < N; i++) { //for particle number i
			Eq_const = E_const * particles.at(i).q * time_dep_term; //the Eq_const ratio do not change for time independant calculations, making it unnecessary to calculate it each time

			F_i[0] = Eq_const * (-2 * r[i][0]); //force at x
			F_i[1] = Eq_const * (-2 * r[i][1]); //force at y
			F_i[2] = Eq_const * (4 * r[i][2]);	//force at z
			F[i] = F_i;
		}
		return F; //in u micro m/(micro s^2 * e), generally unefficent to return array, insted save it as attribute. I do that multiple places
	}





	//Force from external magnetic field:
	std::vector<arma::vec> PenningTrap::MagneticForce(std::vector<arma::vec> v) {
		std::vector<arma::vec> F(N); //list of three-vector that we will return
		arma::vec F_i = arma::vec(3);
		double Bq_const;


		for (int i = 0; i < N; i++) {

			Bq_const = B0 * particles.at(i).q; //again do not change and unneccessary to calculate each time
			F_i[0] = Bq_const * v[i][1]; //force at x
			F_i[1] = -Bq_const * v[i][0]; //force at y
			F_i[2] = 0.; //force at z
			F[i] = F_i;
		}
		return F; //in u/(micro s * e)
	}





	//Force between particles:
	std::vector<arma::vec> PenningTrap::ParticleForce(std::vector<arma::vec> r) {
		std::vector<arma::vec> F(N); //list of three-vector that we will return
		arma::vec F_i = arma::vec(3);

		for (int i = 0; i < N; i++) {
			F[i] = F_i; //making F[i] a (0,0,0), may be done more efficently
		}

		double keq_i; //ke * q
		arma::vec r_i; //the vector i
		arma::vec r_j; //the vector j
		double distance; //distance between ri and rj
		arma::vec distance_vec; //unfortunatly I didn't find a way to save the distance directly as a double
		arma::vec r_relative; //vector from ri to rj

		for (int i = 0; i < N; i++) { //for particle number i

			keq_i = ke * particles.at(i).q; //again constants, may be only calculated onec
			r_i = r[i];

			for (int j = i + 1; j < N; j++) { //using Newtons 3. law, skipping same calculation twice and avoiding i=j
				r_j = r[j];
				r_relative = r_i - r_j;
				
				distance_vec = r_relative.t() * r_relative;
				distance = pow(distance_vec.at(0), 0.5);
				
				F_i = keq_i * particles.at(j).q * r_relative / pow(distance, 3);
				F[i] += F_i;
				F[j] -= F_i; //using Newtons 3. law
			}
		}
		return F; //in u micro m/(micro s^2 * e)
	}







	std::vector<arma::vec> PenningTrap::force(std::vector<arma::vec> r, std::vector<arma::vec> v, double t) {

		std::vector<arma::vec> F(N); //list of three-vector that we will return
		std::vector<arma::vec> F_E = ElectricForce(r, t); //finding all trap to particle electric forces
		std::vector<arma::vec> F_M = MagneticForce(v); //finding all trap to particle magnetic forces

		//May be generally done more efficently by just ignoring paticles outside the trap
		arma::vec distance_vec; //May be unefficent to declere variales each time, which I do a lot in whole class

		/*
		Here we check if we want particle interactions.
		This is highly inneficent as we check each time step, 
		but I didn't manage to assign a variable to a function 
		or make a vector of functions (within reasonable time).
		*/
		if (part_inter) {
			std::vector<arma::vec> F_P = ParticleForce(r);
			for (int i = 0; i < N; i++) {

				//Check if particle is the trap, note the inefficency, as we have already calculated the forces ont the particles, we just don't use them
				distance_vec = r[i].t() * r[i]; //Check if particle is within the trap
				if (pow(distance_vec.at(0), 0.5) > d) {
					F[i] = F_P[i]; //we still want particle interactions outside the trap (altought we don't use it)
					continue;
				}

				F[i] = F_E[i] + F_M[i] + F_P[i]; //total force
			}
		}

		else {
			for (int i = 0; i < N; i++) {

				//Check if the particle is in the trap
				distance_vec = r[i].t() * r[i]; //Check if particle is within the trap
				if (pow(distance_vec.at(0), 0.5) > d) {
					F[i] = arma::vec(3);
					continue;
				}

				F[i] = F_E[i] + F_M[i]; //electric and magnetic force only
			}
		}
		return F;
	}







	void PenningTrap::solveFE(double time_interval, int steps) { //time in micro s
		
		double dt = time_interval / steps;
		std::vector<arma::vec> F;

		//Note the user only may use one solve function per instance, without deleting old values:
		r = std::vector<std::vector<arma::vec>>(steps + 1);
		v = std::vector<std::vector<arma::vec>>(steps + 1);

		//Initial values:
		t.push_back(0.);
		for (int i = 0; i < N; i++) {
			r[0].push_back(particles.at(i).r); //At first step
			v[0].push_back(particles.at(i).v);
		}
		//Forward Euler:
		for (int n = 1; n <= steps; n++) {
			F = force(r[n - 1], v[n - 1], t[n - 1]);
			t.push_back(n * dt);
			for (int i = 0; i < N; i++) {
				r[n].push_back(r[n - 1][i] + v[n - 1][i] * dt);
				v[n].push_back(v[n - 1][i] + F[i] / particles.at(i).m * dt);
			}
		}
	}







	void PenningTrap::solveRK4(double time_interval, int steps) { //time in micro s

		double dt = time_interval / steps;
		std::vector<arma::vec> F;

		r = std::vector<std::vector<arma::vec>>(steps + 1);
		v = std::vector<std::vector<arma::vec>>(steps + 1);

		//Initial values:
		t.push_back(0.);
		for (int i = 0; i < N; i++) {
			r[0].push_back(particles.at(i).r);
			v[0].push_back(particles.at(i).v);
		}


		std::vector<arma::vec> kr1(N);
		std::vector<arma::vec> kr2(N);
		std::vector<arma::vec> kr3(N);
		std::vector<arma::vec> kv1(N);
		std::vector<arma::vec> kv2(N);
		std::vector<arma::vec> kv3(N);

		std::vector<arma::vec> r_halfStep(N); //these are the space steps at half times
		std::vector<arma::vec> v_halfStep(N);

		arma::vec kr_loc;
		arma::vec kv_loc;

		double dt2 = dt / 2.0; //here we may avoid calculating dt/2 multiple times
		double sixth = 1.0 / 6.0; //same for 1/6


		//RK4:
		for (int n = 1; n <= steps; n++) {
			//k1:
			F = force(r[n - 1], v[n - 1], t[n - 1]); //find force
			t.push_back(n * dt); //save the time


			for (int i = 0; i < N; i++) {
				kr_loc = dt * v[n - 1][i];
				kv_loc = dt * F[i] / particles.at(i).m;

				kr1[i] = kr_loc; //store k value
				kv1[i] = kv_loc;

				r_halfStep[i] = r[n - 1][i] + kr_loc / 2.0;
				v_halfStep[i] = v[n - 1][i] + kv_loc / 2.0;
			}

			//k2:
			F = force(r_halfStep, v_halfStep, t[n - 1] + dt2);

			for (int i = 0; i < N; i++) {
				kr_loc = dt * (v[n - 1][i] + F[i] / particles.at(i).m * dt2); //Prediction using FE for time t+1/2*dt
				kv_loc = dt * F[i] / particles.at(i).m;

				kr2[i] = kr_loc;
				kv2[i] = kv_loc;

				r_halfStep[i] = r[n - 1][i] + kr_loc / 2.0;
				v_halfStep[i] = v[n - 1][i] + kv_loc / 2.0;
			}

			//k3:
			F = force(r_halfStep, v_halfStep, t[n - 1] + dt2);

			for (int i = 0; i < N; i++) {
				kr_loc = dt * (v[n - 1][i] + F[i] / particles.at(i).m * dt2);
				kv_loc = dt * F[i] / particles.at(i).m;

				kr3[i] = kr_loc;
				kv3[i] = kv_loc;


				r_halfStep[i] = r[n - 1][i] + kr_loc;
				v_halfStep[i] = v[n - 1][i] + kv_loc;
			}

			//k4:
			F = force(r_halfStep, v_halfStep, t[n - 1] + dt);

			for (int i = 0; i < N; i++) {
				kr_loc = dt * (v[n - 1][i] + F[i] / particles.at(i).m * dt); //Prediction using FE for time t + dt
				kv_loc = dt * F[i] / particles.at(i).m;; //We dont need to save the k4 values

				r[n].push_back(r[n - 1][i] + sixth * (kr1[i] + 2.0 * kr2[i] + 2.0 * kr3[i] + kr_loc)); 
				v[n].push_back(v[n - 1][i] + sixth * (kv1[i] + 2.0 * kv2[i] + 2.0 * kv3[i] + kv_loc));

			}
		}
	}






	double PenningTrap::particlesInsideTrap() {

		arma::vec distance_vec;

		int steps = r.size(); //This is really steps+1 as we also have initial values in r
		double count = 0.0; //It's double to avoid integer division

		for (int i = 0; i < N; i++) {

			distance_vec = r[steps][i].t() * r[steps][i]; //Find the lengdth of position vector

			if (pow(distance_vec.at(0), 0.5) < d) {
				count += 1; //If the position vector is inside the penning trap add 1
			}
		}

		return count;
	}
