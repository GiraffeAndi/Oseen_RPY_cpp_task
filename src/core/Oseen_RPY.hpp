#include <vector>

struct Data_task{
    
    std::vector<std::vector<double>> forces;

    std::vector<std::vector<double>> positions;
    
    std::vector<std::vector<double>> diffusion_tensor;

    //CRND = Correlated random normal diviates
    std::vector<double> CRND;

    double pot_energ, virial, dt;

};

//pot_energ and virial are part of data
void calc_forces(Data_task& data, const double sigma, const double rcut, const double consii, const double consij);

            
void covar(Data_task& data, double dt);

void move(Data_task& data, double dt, double temp);
