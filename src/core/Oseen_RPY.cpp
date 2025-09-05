#include "Oseen_RPY.hpp"

#include "random.hpp"

#include "grid.hpp"

#include <vector>

#include <cmath>

void calc_forces(Data_task& data, const double sigma, const double rcut, const double consii, const double consij){

    const std::size_t N = data.positions.size();

    std::vector<double> RXI(N, 0.0);

    std::vector<double> RYI(N, 0.0);

    std::vector<double> RZI(N, 0.0);

    double FXI, FYI, FZI;

    double FIJ, FXIJ, FYIJ, FZIJ;

    double RXIJ, RYIJ, RZIJ;

    double RIJSQ;

    double RIJ;

    double RRIJSQ;

    double OIJ;

    
    //0 for Oseen
    double RPIJ;
    
    std::size_t IC, JC;
    
    const double rcutsq = rcut * rcut;

    double SR2, SR6, VIJ, WIJ;

    //has to be adjusted formula is totally different to fortran code
    const double sigmacub12 = (sigma * sigma *sigma)/12;

    const double SIGSQ = sigma * sigma;

    //increases readability later
    auto &d = data.diffusion_tensor;

    //setting initial forces to zero
    for(std::size_t k = 0; k < N; k++){

        data.forces[k][0] = 0.0;

        data.forces[k][1] = 0.0;

        data.forces[k][2] = 0.0; 

    };

    for(std::size_t i = 0; i < N - 1; i++){

        RXI[i] = data.positions[i][0];

        RYI[i] = data.positions[i][1];

        RZI[i] = data.positions[i][2];

        //all the same because we zeroed them before
        FXI = data.forces[i][0];

        FYI = data.forces[i][1];
        
        FZI = data.forces[i][2];

        //the counter has to be adjusted 
        IC = 3 * i;
        
        for(std::size_t j = i + 1; j < N; j++){

            Utils::Vector3d pos1 (RXI[i], RYI[i], RZI[i]);
            
            Utils::Vector3d pos2 (RXI[j], RYI[j], RZI[j]);
            
            //see grid.hpp file for box_geo
            Utils::Vector3d sys_dist = box_geo.get_mi_vector(pos1, pos2);

            RXIJ = sys_dist.x;

            RYIJ = sys_dist.y;

            RZIJ = sys_dist.z;
            
            RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ;

            JC = 3 * j;

            RIJ = std::sqrt(RIJSQ);

            //maybe check for 1/0 case?
            //not a priority rn
            RRIJSQ = 1.0/RIJSQ;

            OIJ = consij/RIJ;

            RPIJ = OIJ * sigmacub12 * RRIJSQ;

            //making the tensor symmetric in one loop
            //avoiding extra loop at the end

            d[IC][JC] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RXIJ * RXIJ * RRIJSQ;

            d[JC][IC] = d[IC][JC];

            d[IC+1][JC+1] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RYIJ * RYIJ * RRIJSQ;

            d[JC+1][IC+1] = d[IC+1][JC+1];

            d[IC+2][JC+2] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RZIJ * RZIJ * RRIJSQ;

            d[JC+2][IC+2] = d[IC+2][JC+2];

            d[IC][JC+1] = (OIJ - 3.0 * RPIJ) * RXIJ * RYIJ * RRIJSQ;

            d[JC+1][IC] = d[IC][JC+1];

            d[IC][JC+2] = (OIJ - 3.0 * RPIJ) * RXIJ * RZIJ * RRIJSQ;

            d[JC+2][IC] = d[IC][JC+2];

            d[IC+1][JC+2] = (OIJ - 3.0 * RPIJ) * RYIJ * RZIJ * RRIJSQ;

            d[JC+2][IC+1] = d[IC+1][JC+2];

            if(RIJSQ < rcutsq){

                SR2 = SIGSQ * RRIJSQ;

                SR6 = SR2 * SR2 * SR2;

                VIJ = SR6 * (SR6 - 1.0);

                WIJ = SR6 * (SR6 - 0.5);

                FIJ = WIJ * RRIJSQ;

                FXIJ = FIJ * RXIJ;

                FYIJ = FIJ * RYIJ;

                FZIJ = FIJ * RZIJ;

                data.pot_energ += VIJ;

                data.virial += WIJ;

                FXI += FXIJ;

                FYI += FYIJ;

                FZI += FZIJ;

                data.forces[j][0] -= FXIJ;

                data.forces[j][1] -= FYIJ;

                data.forces[j][2] -= FZIJ;
            }


        }

        data.forces[i][0] = FXI;

        data.forces[i][1] = FYI;

        data.forces[i][2] = FZI;
        
        
    };

    data.pot_energ *= 4.0;

    data.virial *= 16.0;

    for(std::size_t p = 0; p < N; p++){

        data.forces[p][0] *= 48.0;

        data.forces[p][1] *= 48.0;

        data.forces[p][2] *= 48.0;

        //counter again has to be adjusted
        IC = 3 * p;

        d[IC][IC] = consii;

        d[IC+1][IC+1] = consii;

        d[IC+2][IC+2] = consii;

        d[IC][IC+1] = 0.0; 

        d[IC+1][IC] = d[IC][IC+1];
        
        d[IC][IC+2] = 0.0; 

        d[IC+2][IC] = d[IC][IC+2];
        
        d[IC+1][IC+2] = 0.0;

        d[IC+2][IC+1] = d[IC+1][IC+2];
    };

}

void covar(Data_task& data, double dt, BrownianThermostat const &brownian){

    const std::size_t N = data.positions.size();

    std::vector<std::vector<double>> L;

    L.resize(3*N, std::vector<double>(3*N, 0.0));

    std::vector<double> XI(3 * N, 0.0);

    L[0][0] = std::sqrt(data.diffusion_tensor[0][0]);

    L[1][0] = data.diffusion_tensor[1][0]/L[0][0];

    L[1][1] = std::sqrt(data.diffusion_tensor[1][1] - L[1][0] * L[1][0]);

    for(std::size_t i = 2; i < 3 * N; i++){

        L[i][1] = data.diffusion_tensor[i][1]/L[1][1];

        for(std::size_t j = 1; j < i; j++){

            double sum = 0.0;

            for(std::size_t k = 0; k < j; k++){

                sum += L[i][k] * L[j][k];
            }

            L[i][j] = (data.diffusion_tensor[i][j] - sum)/L[j][j];

        }

        double sum = 0.0;

        for(std::size_t p = 0; p < i - 1; p++){

            sum += L[i][p] * L[i][p];
        }

        L[i][i] = std::sqrt(data.diffusion_tensor[i][i] - sum);
    }

    //logic taken from bd_random_walk in the ./thermostats/brownian_inline.hpp file
    auto random_number_seed = brownian.rng_seed();

    auto counter = brownian.rng_counter();

    for(std::size_t f = 0; f < 3 * N; f++){

        //ndrn = normally distributed random number
        //used template from random.hpp instead of seperate gauss function

        //changing the key every loop with loop index to provide different noise

        double ndrn = Random::noise_gaussian<RNGSalt::BROWNIAN_WALK>(counter, random_number_seed, f)[0];

        XI[f] = ndrn * std::sqrt(2.0 * dt);

        double sum = 0.0;

        for(std::size_t w = 0; w <= f; w++){

            sum += L[f][w] * XI[w];
        }

        data.CRND[f] = sum;

    }

}

void move(Data_task& data, double dt, double temp){

    const std::size_t N = data.positions.size();

    std::vector<double> F(3 * N, 0.0);

    for(std::size_t i = 0; i < N; i++){

        std::size_t IC = 3 * i;

        F[IC] = data.forces[i][0];

        F[IC + 1] = data.forces[i][1];

        F[IC + 2] = data.forces[i][2];

    }

    for(std::size_t j = 0; j < N; j++){

        std::size_t JC = 3 * j;

        double sumx = 0.0;

        double sumy = 0.0;

        double sumz = 0.0;

        for(std::size_t c = 0; c < 3 * N; c++){

            sumx += data.diffusion_tensor[JC][c] * F[c];

            sumy += data.diffusion_tensor[JC + 1][c] * F[c];

            sumz += data.diffusion_tensor[JC + 2][c] * F[c];
        }

        data.positions[j][0] += (sumx * dt)/temp + data.CRND[JC];

        data.positions[j][1] += (sumy * dt)/temp + data.CRND[JC + 1];

        data.positions[j][2] += (sumz * dt)/temp + data.CRND[JC + 2];
    }

}

