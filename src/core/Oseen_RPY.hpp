#include <vector>

struct Data_task{
    
    std::vector<std::vector<double>> forces;

    std::vector<std::vector<double>> positions;
    
    std::vector<std::vector<double>> diffusion_tensor;

    //CRND = Correlated random normal diviates
    std::vector<double> CRND;

    double pot_energ, virial, dt;

    //constructor
    Data_task(size_t N) : forces(N, std::vector<double>(3, 0.0)),
    
                          positions(N, std::vector<double>(3, 0.0)),
                                    
                          diffusion_tensor(3*N, std::vector<double>(3*N, 0.0)),

                          CRND(3*N, 0.0),

                          pot_energ(0.0), virial(0.0), dt(0.0)
    {}

};

//pot_energ and virial are part of data
void calc_forces(Data_task& data, const double sigma, const double rcut, const double consii, const double consij);

            
void covar(Data_task& data, double dt, BrownianThermostat const &brownian);

void move(Data_task& data, double dt, double temp);
