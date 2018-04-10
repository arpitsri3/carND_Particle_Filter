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
    if(!is_initialized){
	//"gen" is the random engine initialized here
    default_random_engine gen;
	
	//Initialize number of particles
	num_particles = 69;
    
	//Set the vector size for particles array
    //particles.resize(num_particles);
    
    //Set the vector size for weights array for all particles
    weights.resize(num_particles);

	// These lines create a normal (Gaussian) distribution for x , y and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

    for(int i=0; i<num_particles; i++){
       Particle par;
       par.id = i;
	   par.x = dist_x(gen);
	   par.y = dist_y(gen);
	   par.theta = dist_theta(gen);
	   par.weight = 1.0;
	   particles.push_back(par);
	}
	
	is_initialized = true;
	}
	cout<<"Initialized"<<"\n";
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//"gen" is the random engine initialized here
    default_random_engine gen;
    
	// These line create a normal (Gaussian) distribution for x , y and theta to help in adding noise
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
    
    //Using the equations discussed in classroom(and UKF) to update x, y and yaw_rate
	for(int i=0; i<num_particles; ++i){
	if(fabs(yaw_rate)>0.0001){
        particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + (yaw_rate*delta_t))- sin(particles[i].theta));
	    particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta)- cos(particles[i].theta + (yaw_rate*delta_t)));
	    particles[i].theta += yaw_rate*delta_t;
	}
	else{
		particles[i].x += velocity * delta_t * cos(particles[i].theta);
		particles[i].y += velocity * delta_t * sin(particles[i].theta);
	}

    //Add noise 
	particles[i].x += dist_x(gen);
	particles[i].y += dist_y(gen);
	particles[i].theta += dist_theta(gen);
   }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
	int sel_id;
	
	double x_1 , x_2 , y_1 , y_2 , place_holder , differ;

	for(int i=0; i<observations.size(); i++){
		place_holder = numeric_limits<double>::max();
		for(int j=0; j<predicted.size(); j++){
			x_1 = observations[i].x;
			x_2 = predicted[j].x;
			y_1 = observations[i].y;
			y_2 = predicted[j].y;
			differ = dist(x_1,y_1,x_2,y_2);//sqrt(((x_1 - x_2) * (x_1 - x_2)) + ((y_1 - y_2) * (y_1 - y_2)));
			if(differ < place_holder){
				sel_id = j;
				place_holder = differ;
			}
		}
		observations[i].id = sel_id;
	}
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
    
double std_x = std_landmark[0];
	double std_y = std_landmark[1];
    
	//cout<<"Point1"<<"\n";

    //Iterating over all the particles
	for(int i = 0; i<num_particles; i++){
     
       //Transformation from vehicle to map co-ordinates
	   vector<LandmarkObs> transformed;
	   //Using equation 3.33 for combining translation and rotation

	   for (int j=0; j<observations.size(); j++){
		   LandmarkObs obs;
		   obs.id = observations[j].id;
		   obs.x = particles[i].x +(cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y);
		   obs.y = particles[i].y +(sin(particles[i].theta)*observations[j].x) + (cos(particles[i].theta)*observations[j].y);
		   transformed.push_back(obs); 
	   }
	 
	   //Create a list of landmarks near for the particle
	   vector<LandmarkObs> near;
	   double x1, x2, y1, y2;//, difference;
       for(int j=0; j<map_landmarks.landmark_list.size(); j++){
           x1 = map_landmarks.landmark_list[j].x_f;
		   x2 = particles[i].x;
		   y1 = map_landmarks.landmark_list[j].y_f;
		   y2 = particles[i].y;
		   
		   if(dist(x1,y1,x2,y2) < sensor_range){
			   LandmarkObs mark;
			   mark.id = map_landmarks.landmark_list[j].id_i;
			   mark.x = map_landmarks.landmark_list[j].x_f;
			   mark.y = map_landmarks.landmark_list[j].y_f;
			   near.push_back(mark);
		   }
	   }
	 


	   //call dataAssociation function to Find the predicted measurement that is closest to each observed measurement and assign the 
	   //observed measurement to this particular landmark.
       dataAssociation(near, transformed);


	   //Calculate parameters for estabilishing associations for the particle
	   vector<int> associations;
	   vector<double> sense_x;
	   vector<double> sense_y;

	   for(int j=0; j<transformed.size(); j++){
		   sense_x.push_back(transformed[j].x);
		   sense_y.push_back(transformed[j].y);
		   associations.push_back(near[transformed[j].id].id);
	   }


	   SetAssociations(particles[i], associations, sense_x, sense_y);


	   //calculate normalization term
       double gauss_norm = 1/(2 * M_PI * std_x * std_y);
       //Taking product of the likelihoods over all measurements
       
	   //Initializing the particle's weight to 1.0
	   double wt = 1.0;
	   for(int j=0; j<observations.size(); j++){
		   double x_obs = transformed[j].x;
		   double mu_x = near[transformed[j].id].x;
		   double y_obs = transformed[j].y;
		   double mu_y = near[transformed[j].id].y;
			
		   //cout<<x_obs<<","<<mu_x<<" "<<y_obs<<" "<<mu_y<<endl;

           //calculate exponent
           double exponent= (((x_obs - mu_x)*(x_obs - mu_x))/(2 * std_x *std_x)) + (((y_obs - mu_y)*(y_obs - mu_y))/(2 * std_y *std_y));
           wt*= gauss_norm * exp(-(exponent));

		   //cout<<weight<<endl;

           transformed.clear();
		   near.clear();
	   }
	   
	   //cout<<particles[i].weight<<endl;

	   //Updating the weight in the weights array for the particle
	   particles[i].weight = wt;
	   weights[i] = wt;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    //Create vector for resampled particles
	vector<Particle> resampled_particles;

	//"gen" is the random engine initialized here
    default_random_engine gen;
    
	//Create a placeholder for storing the largest weight
    double placehold = numeric_limits<double>::min();

	for(int i=0; i<num_particles; i++){
		if(particles[i].weight > placehold){
			placehold = particles[i].weight;
		}
	}
    
	//Create an evenly distributed range for particles
	uniform_int_distribution<int> dis_1(0, num_particles-1);
    
	//Create an evenly distributed range for weights
    uniform_real_distribution<double> dis_2(0.0, placehold);

    int index = dis_1(gen);

    double beta = 0.0;

    for( int i=0; i<num_particles; i++){
         beta += dis_2(gen) * 2.0 * placehold;
         while(beta > weights[index]){
               beta -= weights[index];
               index = (index+1)%num_particles;
		 }
     resampled_particles.push_back(particles[index]);
	}
    
	particles = resampled_particles;
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

	return particle;
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
