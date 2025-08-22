#include "Oseen_RPY.hpp"

#include "random.hpp"

#include "BoxGeometry.hpp"

#include <utils/Vector.hpp>

void calc_forces(Data_task& data, const double sigma, const double rcut, const double consii, const double consij){

    int N = data.positions.size();

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
    
    int IC, JC;
    
    double rcutsq = rcut * rcut;

    double SR2, SR6, VIJ, WIJ;

    //has to be adjusted formula is totally different to fortran code
    double sigmacub12 = (sigma * sigma *sigma)/12;

    double SIGSQ = sigma * sigma;

    //setting initial forces to zero
    for(int k = 0; k < N; k++){

        data.forces[k][0] = 0.0;

        data.forces[k][1] = 0.0;

        data.forces[k][2] = 0.0; 

        RXI[k] = data.positions[k][0];

        RYI[k] = data.positions[k][1];

        RZI[k] = data.positions[k][2];

    };

	//dummy values as placeholder until we know how to use the values from
    //system.box_length ...
    double So = 1.0
     
    double me = 2.0

    double thing = 3.0

    //not correct like this
    //we need the system information passed through the python interface
    //not sure where I can find it           
    BoxGeometry box (So, me, thing);



    for(int i = 0; i < N; i++){

        //all the same because we zeroed them before
        FXI = data.forces[i][0];

        FYI = data.forces[i][1];
        
        FZI = data.forces[i][2];

        //the counter has to be adjusted 
        IC = 3 * i;
        
        for(int j = i + 1; j < N; j++){

            Vector3d pos1 (RXI[i], RYI[i], RZI[i]);
            
            Vector3d pos2 (RXI[j], RYI[j], RZI[j]);
            
            Vector3d sys_dist = box.get_mi_vector(pos1, pos2);

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

            data.diffusion_tensor[IC][JC] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RXIJ * RXIJ * RRIJSQ;

            data.diffusion_tensor[IC+1][JC+1] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RYIJ * RYIJ * RRIJSQ;

            data.diffusion_tensor[IC+2][JC+2] = OIJ + RPIJ + (OIJ - 3.0 * RPIJ) * RZIJ * RZIJ * RRIJSQ;

            data.diffusion_tensor[IC][JC+1] = (OIJ - 3.0 * RPIJ) * RXIJ * RYIJ * RRIJSQ;

            data.diffusion_tensor[IC][JC+2] = (OIJ - 3.0 * RPIJ) * RXIJ * RZIJ * RRIJSQ;

            data.diffusion_tensor[IC+1][JC+2] = (OIJ - 3.0 * RPIJ) * RYIJ * RZIJ * RRIJSQ;

            data.diffusion_tensor[IC+1][JC] = data.diffusion_tensor[IC][JC+1];

            data.diffusion_tensor[IC+2][JC] = data.diffusion_tensor[IC][JC+2];

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

    data.virial *= 48.0/3.0;

    for(int p = 0; p < N; p++){

        data.forces[p][0] *= 48.0;

        data.forces[p][1] *= 48.0;

        data.forces[p][2] *= 48.0;

        IC = 3 * p;

        data.diffusion_tensor[IC][IC] = consii;

        data.diffusion_tensor[IC+1][IC+1] = consii;

        data.diffusion_tensor[IC+2][IC+2] = consii;

        data.diffusion_tensor[IC][IC+1] = 0.0;

        data.diffusion_tensor[IC][IC+2] = 0.0;

        data.diffusion_tensor[IC+1][IC+2] = 0.0;
    };

    for(int r = 0; r < 3 * N - 1; r++){

        for(int t = r + 1; t < 3 * N; t++){

            data.diffusion_tensor[t][r] = data.diffusion_tensor[r][t];
        }
    };



};

void covar(Data_task& data, double dt){

    int N = data.positions.size();

    std::vector<std::vector<double>> L;

    L.resize(3*N, std::vector<double>(3*N, 0.0));

    std::vector<double> XI(3 * N, 0.0);

    L[0][0] = std::sqrt(data.diffusion_tensor[0][0]);

    L[1][0] = data.diffusion_tensor[1][0]/L[0][0];

    L[1][1] = std::sqrt(data.diffusion_tensor[1][1] - L[1][0] * L[1][0]);

    for(int i = 2; i < 3 * N; i++){

        L[i][1] = data.diffusion_tensor[i][1]/L[1][1];

        for(int j = 1; j < i; j++){

            double sum = 0.0;

            for(int k = 0; k < j; k++){

                sum += L[i][k] * L[j][k];
            }

            L[i][j] = (data.diffusion_tensor[i][j] - sum)/L[j][j];

        }

        double sum = 0.0;

        for(int p = 0; p < i - 1; p++){

            sum += L[i][p] * L[i][p];
        }

        L[i][i] = std::sqrt(data.diffusion_tensor[i][i] - sum);
    }

    //randomly chosen by me
    int random_number_seed = 11;

    for(int f = 0; f < 3 * N; f++){

        //ndrn = normally distributed random number
        //used template from random.hpp instead of seperate gauss function

        //changing the key every loop with loop index to provide different noise

        double ndrn = Random::noise_gaussian<RNGSalt::BROWNIAN_WALK, 1>(dt, random_number_seed, f)[0];

        XI[f] = ndrn * std::sqrt(2.0 * dt);

        double sum = 0.0;

        for(int w = 0; w < f; w++){

            sum += L[f][w] * XI[w];
        }

        data.CRND[f] = sum;

    }

};

void move(Data_task& data, double dt, double temp){

    int N = data.positions.size();

    std::vector<double> F(3 * N, 0.0);

    for(int i = 0; i < N; i++){

        int IC = 3 * i;

        F[IC] = data.forces[i][0];

        F[IC + 1] = data.forces[i][1];

        F[IC + 2] = data.forces[i][2];

    }

    for(int j = 0; j < N; j++){

        int JC = 3 * j;

        double sumx = 0.0;

        double sumy = 0.0;

        double sumz = 0.0;

        for(int c = 0; c < 3 * N; c++){

            sumx += data.diffusion_tensor[JC][c] * F[c];

            sumy += data.diffusion_tensor[JC + 1][c] * F[c];

            sumz += data.diffusion_tensor[JC + 2][c] * F[c];
        }

        data.positions[j][0] += (sumx * dt)/temp + data.CRND[JC];

        data.positions[j][1] += (sumy * dt)/temp + data.CRND[JC + 1];

        data.positions[j][2] += (sumz * dt)/temp + data.CRND[JC + 2];
    }

};

