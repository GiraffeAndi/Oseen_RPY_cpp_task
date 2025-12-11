#include "Oseen_RPY.hpp"

#include "random.hpp"

#include "grid.hpp"

#include <vector>

#include <cmath>

void calc_forces(Data_task& data, const double sigma, const double rcut, const std::vector<double> consii, const std::vector<double> consij, const std::vector<double> consii_rot, const std::vector<double> consij_rot){

    const std::size_t N = data.positions.size();

    std::vector<double> RXI(N, 0.0);

    std::vector<double> RYI(N, 0.0);

    std::vector<double> RZI(N, 0.0);

    std::vector<std::vector<double>> rr_block(3, std::vector<double>(3, 0.0));

    std::vector<std::vector<double>> rr(3, std::vector<double>(3, 0.0));

    std::vector<std::vector<double>> levi_civita(3, std::vector<double>(3, 0.0));

    double FXI, FYI, FZI;

    double FIJ, FXIJ, FYIJ, FZIJ;

    double RXIJ, RYIJ, RZIJ;

    double RIJSQ;

    double RIJ;

    double RRIJSQ;

    double OIJ;

    //0 for Oseen
    double RPIJ;
    
    std::size_t IC, JC, IC6, JC6;
    
    const double rcutsq = rcut * rcut;

    double SR2, SR6, VIJ, WIJ;

    const double sigmacub12 = (sigma * sigma *sigma)/12;

    const double SIGSQ = sigma * sigma;

    //increases readability later
    auto &d = data.diffusion_tensor;

    //setting initial forces to zero setting the initial positions
    for(std::size_t k = 0; k < N; k++){

        data.forces[k][0] = 0.0;

        data.forces[k][1] = 0.0;

        data.forces[k][2] = 0.0;
        
        RXI[k] = data.positions[k][0];

        RYI[k] = data.positions[k][1];

        RZI[k] = data.positions[k][2];

    };

    for(std::size_t i = 0; i < N - 1; i++){

        //all the same because we zeroed them before
        FXI = data.forces[i][0];

        FYI = data.forces[i][1];
        
        FZI = data.forces[i][2];

        //the counter has to be adjusted 
        IC = 3 * i;
        
        for(std::size_t j = i + 1; j < N; j++){

            Utils::Vector3d pos1 {RXI[i], RYI[i], RZI[i]};
            
            Utils::Vector3d pos2 {RXI[j], RYI[j], RZI[j]};
            
            //see grid.hpp file for box_geo
            Utils::Vector3d sys_dist = box_geo.get_mi_vector(pos1, pos2);

            RXIJ = sys_dist[0];

            RYIJ = sys_dist[1];

            RZIJ = sys_dist[2];
            
            RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ;

            JC = 3 * j;

            RIJ = std::sqrt(RIJSQ);

            //potential error
            RRIJSQ = 1.0/RIJSQ;

            OIJ = consij[i]/RIJ;

            RPIJ = OIJ * sigmacub12 * RRIJSQ;

            //assessing Oseen Case
            
            #ifdef OSEEN

            RPIJ = 0.0;

            #endif

            //translational part

            //defining constants to maintain readability

            double const1 = OIJ + RPIJ;

            double const2 = (OIJ - 3.0 * RPIJ) * RRIJSQ;

            d[IC][JC] = const1 + const2 * RXIJ * RXIJ;

            d[JC][IC] = d[IC][JC];

            d[IC+1][JC+1] = const1 + const2 * RYIJ * RYIJ;

            d[JC+1][IC+1] = d[IC+1][JC+1];

            d[IC+2][JC+2] = const1 + const2 * RZIJ * RZIJ;

            d[JC+2][IC+2] = d[IC+2][JC+2];

            d[IC][JC+1] = const2 * RXIJ * RYIJ;

            d[JC+1][IC] = d[IC][JC+1];

            d[IC][JC+2] = const2 * RXIJ * RZIJ;

            d[JC+2][IC] = d[IC][JC+2];

            d[IC+1][JC+2] = const2 * RYIJ * RZIJ;

            d[JC+2][IC+1] = d[IC+1][JC+2];

            //mixed and rotational part

            IC6 = 6 * i;

            JC6 = 6 * j;

            double XIJ_hat = RXIJ/RIJ;

            double YIJ_hat = RYIJ/RIJ;

            double ZIJ_hat = RZIJ/RIJ;

            rr[0][0] = XIJ_hat * XIJ_hat;

            rr[0][1] = XIJ_hat * YIJ_hat;

            rr[0][2] = XIJ_hat * ZIJ_hat;

            rr[1][0] = YIJ_hat * XIJ_hat;

            rr[1][1] = YIJ_hat * YIJ_hat;

            rr[1][2] = YIJ_hat * ZIJ_hat;

            rr[2][0] = ZIJ_hat * XIJ_hat;

            rr[2][1] = ZIJ_hat * YIJ_hat;

            rr[2][2] = ZIJ_hat * ZIJ_hat;

            levi_civita[0][1] = ZIJ_hat;

            levi_civita[0][2] = - YIJ_hat;

            levi_civita[1][0] = - ZIJ_hat;

            levi_civita[1][2] = XIJ_hat;

            levi_civita[2][0] = YIJ_hat;

            levi_civita[2][1] = -XIJ_hat;

            for (std::size_t a = 0; a < 3; a++){

                for (std::size_t b = 0; b < 3; b++){

                    //rr
                    //a is particle radius => sigma/2
                    d[IC6 + 3 + a][JC6 + 3 + b] =  - consij_rot[i] * (sigma/2)/RIJ * (sigma/2)/RIJ * (sigma/2)/RIJ * (((a == b) ? 1.0 : 0.0)- 3 * rr[a][b]);

                    d[JC6 + 3 + b][IC6 + 3 + a] = d[IC6 + 3 + a][JC6 + 3 + b];

                    //rt
                    d[IC6 + 3 + a][JC6 + b] = levi_civita[a][b] * consij_rot[i] * (sigma/2 * sigma/2 * sigma/2)/(RIJ * RIJ);

                    d[JC6 + 3 + a][IC6 + b] = - levi_civita[b][a] * consij_rot[i] * (sigma/2 * sigma/2 * sigma/2)/(RIJ * RIJ);


                    //tr
                    d[IC6 + a][JC6 + 3 + b] = d[IC6 + 3 + a][JC6 + b];

                    d[JC6 + a][IC6 + 3 + b] = d[JC6 + 3 + a][IC6 + b];

                }
            }

            if(RIJSQ < rcutsq){

                SR2 = SIGSQ * RRIJSQ;

                SR6 = SR2 * SR2 * SR2;

                //WCA adaptation V(rcut) = 0
                VIJ = SR6 * (SR6 - 1.0) + 0.25;

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

        //Oseen has no rotations
        IC = 3 * p;

        IC6 = 6 * p;

        d[IC][IC] = consii[p];

        d[IC+1][IC+1] = consii[p];

        d[IC+2][IC+2] = consii[p];

        d[IC][IC+1] = 0.0; 

        d[IC+1][IC] = d[IC][IC+1];
        
        d[IC][IC+2] = 0.0; 

        d[IC+2][IC] = d[IC][IC+2];
        
        d[IC+1][IC+2] = 0.0;

        d[IC+2][IC+1] = d[IC+1][IC+2];

        d[IC6 + 3][IC6 + 3] = consii_rot[p];

        d[IC6 + 4][IC6 + 4] = consii_rot[p];

        d[IC6 + 5][IC6 + 5] = consii_rot[p];

        d[IC6 + 3][IC6 + 4] = 0.0;

        d[IC6 + 3][IC6 + 5] = 0.0;

        d[IC6 + 4][IC6 + 3] = 0.0;

        d[IC6 + 4][IC6 + 5] = 0.0;

        d[IC6 + 5][IC6 + 3] = 0.0;

        d[IC6 + 5][IC6 + 4] = 0.0;

    };

}

void covar(Data_task& data, double dt, BrownianThermostat const &brownian, double kT){

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

        XI[f] = ndrn * std::sqrt(2.0 * dt * kT);

        double sum = 0.0;

        for(std::size_t w = 0; w <= f; w++){

            sum += L[f][w] * XI[w];
        }

        data.CRND[f] = sum;

    }

}

void move(Data_task& data, double dt){

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

        data.positions[j][0] += sumx * dt + data.CRND[JC];

        data.positions[j][1] += sumy * dt + data.CRND[JC + 1];

        data.positions[j][2] += sumz * dt + data.CRND[JC + 2];
    }

}

