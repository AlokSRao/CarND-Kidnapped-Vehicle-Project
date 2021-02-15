/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <random>

#include "helper_functions.h"


using namespace std;
using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  //safety check to ensure init is called only once
  if (is_initialized) {
    return;
  }

  // Initializing the number of particles
  num_particles = 100;

  // Extracting standard deviations
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Creating normal distributions
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  // Generate particles with normal distribution with mean on GPS values.
  for (int i = 0; i < num_particles; i++) {

    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
    
  }

  // The filter is now initialized.
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  // Extracting standard deviations
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // Creating normal distributions
  std::normal_distribution<double> dist_x(0, std_x);
  std::normal_distribution<double> dist_y(0, std_y);
  std::normal_distribution<double> dist_theta(0, std_theta);

  // Calculate new state.
  for (int i = 0; i < num_particles; i++) {

  	double theta = particles[i].theta;
	//if yaw doesnt change
    if ( fabs(yaw_rate) < 0.001) 
    {
      particles[i].x += velocity * delta_t * cos( theta );
      particles[i].y += velocity * delta_t * sin( theta );
     
    } 
    else {
      particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
      particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate * delta_t;
    }

    //Adding noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {

  int number_of_Observations = observations.size();
  int number_of_Predictions = predicted.size();

  // Loop through observations

  for (int i = 0; i < number_of_Observations; i++) 
  { 
    // Initialize min distance to maximum value
    double minDistance = numeric_limits<double>::max();

    //initialising id for observation
    int observation_id = -1;

    //looping through each predicition
    for (int j = 0; j < number_of_Predictions; j++ ) 
    {

      double xDistance = observations[i].x - predicted[j].x;
      double yDistance = observations[i].y - predicted[j].y;

      double distance = xDistance * xDistance + yDistance * yDistance;

      // If the "distance" is less than min, stored the id and update min.
      if ( distance < minDistance ) {
        minDistance = distance;
        observation_id = predicted[j].id;
      }
    }

    // Update the observation identifier.
    observations[i].id = observation_id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  double std_landmark_range = std_landmark[0];
  double std_landmark_bearing = std_landmark[1];

  for (int i = 0; i < num_particles; i++) 
  {

    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    // Find landmarks in particle's range.
    double sensor_range_2 = sensor_range * sensor_range;
    vector<LandmarkObs> inRangeLandmarks;
    for(int j = 0; j < map_landmarks.landmark_list.size(); j++) 
    {
      float landmarkX = map_landmarks.landmark_list[j].x_f;
      float landmarkY = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dX = x - landmarkX;
      double dY = y - landmarkY;
      if ( dX*dX + dY*dY <= sensor_range_2 ) {
        inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
      }
    }

    // Transform observation coordinates.
    vector<LandmarkObs> mappedObservations;
    for(int j = 0; j < observations.size(); j++) 
    {
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mappedObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
    }

    // Observation association to landmark.
    dataAssociation(inRangeLandmarks, mappedObservations);

    // Reseting weight.
    particles[i].weight = 1.0;
    // Calculate weights.
    for(int j = 0; j < mappedObservations.size(); j++) 
    {
      double observationX = mappedObservations[j].x;
      double observationY = mappedObservations[j].y;

      int landmarkId = mappedObservations[j].id;

      double landmarkX, landmarkY;
      int k = 0;
      int number_of_Landmarks = inRangeLandmarks.size();
      bool found = false;
      while( !found && k < number_of_Landmarks ) 
      {
        if ( inRangeLandmarks[k].id == landmarkId) {
          found = true;
          landmarkX = inRangeLandmarks[k].x;
          landmarkY = inRangeLandmarks[k].y;
        }
        k++;
      }

      // Calculating weight.
      double dX = observationX - landmarkX;
      double dY = observationY - landmarkY;

      double weight = ( 1/(2*M_PI*std_landmark_range*std_landmark_bearing)) * exp( -( dX*dX/(2*std_landmark_range*std_landmark_range) + (dY*dY/(2*std_landmark_bearing*std_landmark_bearing)) ) );
      if (weight == 0) 
      {
        particles[i].weight *= 0.0001;
      }
      else
      {
        particles[i].weight *= weight;
      }
    }
  }
}

void ParticleFilter::resample() {
  vector<double> weights;
  double maxWeight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > maxWeight ) {
      maxWeight = particles[i].weight;
    }
  }

  // Creating distributions.
  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generating index.
  int index = distInt(gen);

  double beta = 0.0;

  // wheel algorithm
  vector<Particle> resampledParticles;
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}