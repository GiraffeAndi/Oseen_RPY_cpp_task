#ifndef OSEEN_RPY_HPP
#define OSEEN_RPY_HPP

#include <vector>
#include "thermostat.hpp"

struct Data_task{
    
    std::vector<double> forces;
    std::vector<double> positions;
    std::vector<double> velocities;
    std::vector<std::vector<double>> diffusion_tensor;
    //CRND = Correlated random normal diviates
    std::vector<double> CRND;

    //constructor
    Data_task(std::size_t N) : 

        forces(6*N, 0.0),
        positions(6*N, 0.0),
        velocities(6*N, 0.0),           
        diffusion_tensor(6*N, std::vector<double>(6*N, 0.0)),
        CRND(6*N, 0.0)   
    {}

};

void calc_mobility_matrix(Data_task& data, const double sigma, const std::vector<double>& cons_tr, const std::vector<double>& cons_rot);
void calc_velocities(Data_task& data);
void dump_velocity_vector(const std::vector<double>& d);
void dump_diffusion_tensor(const std::vector<std::vector<double>>& d);

#endif