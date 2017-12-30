/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 80;//set the number of particles

	default_random_engine gen;

	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);

	for(int i=0;i<num_particles;i++)
	{
		weights.push_back(1);
		Particle par;
		par.id = i;
		par.x = dist_x(gen);
		par.y = dist_y(gen);
		par.theta = dist_theta(gen);
		par.weight=1;
		particles.push_back(par);
	}
	is_initialized = true;
	//debug
	// cout<<"finishi init"<<endl<<"the x:"<<endl;
	// for(int kk=0;kk<particles.size();kk++)
	// {cout<<particles[kk].x<<','<<weights[kk]<<endl;}
	
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//debug
	// cout<<'v'<<velocity<<",yawrate"<<yaw_rate<<endl;
	// for(int i=0;i<num_particles;i++)
	// {
	// 	cout<<"p:"<<particles[i].id<<",,"<<particles[i].x<<endl;
	// }

	default_random_engine gen;
	for(int i=0;i<num_particles;i++)
	{
		
		double xx = particles[i].x;
		double yy = particles[i].y;
		double tt = particles[i].theta;
		if(yaw_rate!=0)
		{
			xx = xx + velocity / yaw_rate *(sin(tt + yaw_rate*delta_t) - sin(tt));
			yy = yy + velocity / yaw_rate *(cos(tt) - cos(tt + yaw_rate*delta_t));
			tt = tt + yaw_rate*delta_t;
			//debug
			//cout<<xx<<','<<yy<<','<<tt<<','<<endl;
		}
		else
		{
			xx = xx + velocity*cos(tt)*delta_t;
			yy = yy + velocity*sin(tt)*delta_t;
			//debug
			//cout<<'0'<<xx<<','<<yy<<','<<tt<<','<<endl;
		}

		normal_distribution<double> dist_x(xx,std_pos[0]);
		normal_distribution<double> dist_y(yy,std_pos[1]);
		normal_distribution<double> dist_theta(tt,std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
	//debug
	// cout<<"finishi predict"<<endl;
	// for(int kk=0;kk<particles.size();kk++)
	// {cout<<particles[kk].weight<<"::"<<particles[kk].x<<endl;}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//get the sense data

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// get the predict sense
	for(int i=0; i<num_particles; i++){
		Particle par;
		par.id = particles[i].id;
		par.x= particles[i].x;
		par.y=particles[i].y;
		par.theta=particles[i].theta;
		for(int j=0; j<observations.size(); j++){
			LandmarkObs obm = observations[j];
			double sx = par.x + obm.x*cos(par.theta) - obm.y*sin(par.theta);
			double sy = par.y + obm.x*sin(par.theta) + obm.y*cos(par.theta);
			par.sense_x.push_back(sx);
			par.sense_y.push_back(sy);
			double min_dis=sensor_range+10.0;
			int asso_index=0;
			for(int k=0;k<map_landmarks.landmark_list.size();k++){
				if(dist(par.x,par.y,
					map_landmarks.landmark_list[k].x_f,
					map_landmarks.landmark_list[k].y_f)>sensor_range)
				{}
				else{
					double dis = dist(sx,sy,
						map_landmarks.landmark_list[k].x_f,
						map_landmarks.landmark_list[k].y_f);
					if(dis<min_dis){
						min_dis=dis;
						asso_index=k;
					}
				}
			}
			par.associations.push_back(asso_index+1);
		}
		double weight=1.0;
		for(int m=0;m<par.associations.size();m++){
			int mm = par.associations[m];
			weight*=MG_prob(par.sense_x[m],
				par.sense_y[m],
				map_landmarks.landmark_list[mm-1].x_f,
				map_landmarks.landmark_list[mm-1].y_f,std_landmark);
		}
		particles[i].weight=weight;
		particles[i].sense_x=par.sense_x;
		particles[i].sense_y=par.sense_y;
		particles[i].associations=par.associations;
		weights[i]=weight;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> ddd(weights.begin(), weights.end());
	std::vector<Particle> partemp;

	for(int i=0;i<num_particles;i++)
	{
		int k = ddd(gen);
		Particle par;
		par=particles[k];
		par.id=i;
		partemp.push_back(par);
	}
	particles=partemp;
	//debug
	// cout<<"finish resample"<<endl;
	// for(int kk=0;kk<particles.size();kk++)
	// {cout<<particles[kk].weight<<"::"<<particles[kk].x<<endl;}
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
