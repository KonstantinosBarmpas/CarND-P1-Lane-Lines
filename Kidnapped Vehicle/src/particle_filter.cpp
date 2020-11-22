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
	
    //Random generator
    default_random_engine gen;
    
    //Number of particles  --- 500 is a good number for initial particles for this project
    num_particles = 500;
    
    //Create normal Gaussian distributions for x, y and theta using the GPS and standard deviations
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i=0; i<num_particles; i++){
        
        //call Particle constructor
        Particle particle;
        
        //set parameters to the instance of the class
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        particle.id = i;
        
        //Add them to the two vectors
        weights.push_back(1.0);
        particles.push_back(particle);
        
    }
    //Set to true. Initiallization done.
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	
    //Random generator
    default_random_engine gen;
    
    for (int i=0; i<num_particles; i++){
        
        //Use the formula to calculate the new x,y,theta
        if (yaw_rate != 0){
            particles[i].x = particles[i].x + ( velocity / yaw_rate )*(sin(particles[i].theta + delta_t*yaw_rate) - sin(particles[i].theta));
            particles[i].y = particles[i].y + ( velocity / yaw_rate )*(cos(particles[i].theta) - cos(particles[i].theta + delta_t*yaw_rate));
            particles[i].theta = particles[i].theta + yaw_rate * delta_t;
        }else{
            particles[i].x = particles[i].x + velocity * delta_t * cos( particles[i].theta);
            particles[i].y = particles[i].y + velocity * delta_t * sin( particles[i].theta);
        }
        
        //Create normal Gaussian distributions for x, y and theta using the GPS and standard deviations
        normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
        
        //Add the random Guassian noise
        particles[i].x = dist_x (gen);
        particles[i].y = dist_y (gen);
        particles[i].theta = dist_theta (gen);
    
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	
    //Iterate through all observations
    for (int i=0; i<observations.size(); i++){
        
        //Say that minimum is predicted[0] and associate the id's
        double minimum_distance = distance(observations[i].x,observations[i].y,predicted[0].x,predicted[0].y);
        observations[i].id = predicted[0].id;
        
        //Iterate through predicted to see if there is closer value
        for (int j=1; j<predicted.size(); j++){
            double dist = distance(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
            
            //If yes, change the id's
            if (dist < minimum_distance){
                minimum_distance = dist ;
                observations[i].id = predicted[j].id;
            }
            
        }
        
    }
}

//Helping function to calculate distance between two points -- Didnt have to define it here there is a similar one to helperFunctions.h
double ParticleFilter::distance (double x1,double y1,double x2,double y2){
    double distancex = (x2 - x1) * (x2 - x1);
    double distancey = (y2 - y1) * (y2 - y1);
    double distance = sqrt(distancex + distancey);
    return distance;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	
    //Iterate through all particles
    for (int number_of_particle=0; number_of_particle<num_particles; number_of_particle++){
    
        //Create the vector with the LandmarkObs predicted landmarks. No transformation needed.
        vector<LandmarkObs> predicted_landmarks;
        for (int i=0; i<map_landmarks.landmark_list.size(); i++){
        
            //Define the new LandmarkObs
            LandmarkObs newLandOb;
            newLandOb.id = map_landmarks.landmark_list[i].id_i;
            newLandOb.x = map_landmarks.landmark_list[i].x_f;
            newLandOb.y = map_landmarks.landmark_list[i].y_f;
        
            //If its inside the sensor range append it
            if (fabs(map_landmarks.landmark_list[i].x_f - particles[number_of_particle].x) <= sensor_range &&        fabs(map_landmarks.landmark_list[i].y_f - particles[number_of_particle].y)<=sensor_range){
                predicted_landmarks.push_back(newLandOb);
            }
        }
    
        //Create the observation Landmarks by transforming them
        std::vector<LandmarkObs> observed_landmarks;
        for (int i=0; i<observations.size(); i++){
            LandmarkObs newLandOb2;
            newLandOb2.id = observations[i].id;
            newLandOb2.x = particles[number_of_particle].x + cos(particles[number_of_particle].theta)*observations[i].x - sin(particles[number_of_particle].theta)*observations[i].y;
            newLandOb2.y = particles[number_of_particle].y + sin(particles[number_of_particle].theta)*observations[i].x + cos(particles[number_of_particle].theta)*observations[i].y;
            observed_landmarks.push_back(newLandOb2);
        }
        
        //Call the dataAssociation function
        dataAssociation(predicted_landmarks,observed_landmarks);
        
        //Declare the vectors we will fill in order to call setAssociations function.
        std::vector<int> associations;
        std::vector<double> sense_x;
        std::vector<double> sense_y;
        
        //Set the weigth initially to one in order to do all the mutiplications
        particles[number_of_particle].weight = 1.0;
        
        //Iterate through all observed and predicted landmarks. If matching ids do the multivariable Guassian distribution to update the weights.
        
        for (int i=0; i<observed_landmarks.size(); i++){
            
            for (int j=0; j<predicted_landmarks.size(); j++){
                
                if (observed_landmarks[i].id == predicted_landmarks[j].id){
                    
                    double gaussian_normalizer = (1 / (2*M_PI*std_landmark[0]*std_landmark[1]));
                    double guassian_multi_exp = pow(observed_landmarks[i].x-predicted_landmarks[j].x,2) / (2*pow(std_landmark[0],2)) + pow(observed_landmarks[i].y-predicted_landmarks[j].y,2) / (2*pow(std_landmark[1],2));
                    double gaussian_multiplier = gaussian_normalizer * exp (-1 * guassian_multi_exp);
                    if (gaussian_multiplier > 0.0){ // bigger we dont want to lose all info in case of an error
                        particles[number_of_particle].weight = particles[number_of_particle].weight * gaussian_multiplier;
                    }
                    associations.push_back(observed_landmarks[i].id);
                    sense_x.push_back(observed_landmarks[i].x);
                    sense_y.push_back(observed_landmarks[i].y);
                }
            }
        }
       //Update particles and weigths
       particles[number_of_particle] = SetAssociations (particles[number_of_particle],associations,sense_x,sense_y);
       weights[number_of_particle]= particles[number_of_particle].weight;
     }
        
}

void ParticleFilter::resample() {
	//Random generator
    default_random_engine gen;
    
    //probability proportional to their weight
    discrete_distribution<int> proportional_distribution(weights.begin(),weights.end());
    
    vector<Particle> newParticles;
    
    for (int i=0; i<num_particles; i++){
        newParticles.push_back(particles[proportional_distribution(gen)]);
    }
    
    particles=newParticles;

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
