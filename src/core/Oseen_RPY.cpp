#include "Oseen_RPY.hpp"
#include "random.hpp"
#include "grid.hpp"
#include <vector>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <fstream>

void dump_diffusion_tensor(const std::vector<std::vector<double>>& d){
    
    std::ofstream out("tensor_analysis_file.csv", std::ios::app);
    out << std::setprecision(17);
    out << "# new call\n";

    for (std::size_t i = 0; i < 12; i++) {
        
        for (std::size_t j = 0; j < 12; j++) {
            
            out << d[i][j];
            
            if (j + 1 < 12){
                out << ",";
            } 
        }
        out << "\n";
    }

    out << "\n";
    out.close();
}

void dump_velocity_vector(const std::vector<double>& d){
    
    std::ofstream out("velocity_analysis_file.csv", std::ios::app);
    out << std::setprecision(17);
    out << "# new call\n";

    for (std::size_t i = 0; i < d.size(); i++) {

        out << d[i];

        if (i + 1 < d.size()){
            out << ",";
        }
    }

    out << "\n\n";
    out.close();
}

void calc_mobility_matrix(Data_task& data, const double sigma, const std::vector<double>& cons_tr, const std::vector<double>& cons_rot){

    const std::size_t N = data.positions.size()/6;

    std::vector<double> RXI(N, 0.0);
    std::vector<double> RYI(N, 0.0);
    std::vector<double> RZI(N, 0.0);

    std::vector<std::vector<double>> rr(3, std::vector<double>(3, 0.0));
    std::vector<std::vector<double>> levi_civita(3, std::vector<double>(3, 0.0));

    double RXIJ, RYIJ, RZIJ;
    double RIJ, RIJ3, RIJSQ;

    //unit matrix
    double I;

    std::size_t IC, JC, IC6, JC6;

    auto radi = sigma/2.0;
    const double coeff1 = 3.0 * radi/4.0;
    const double coeff3 = radi * radi * radi;
    const double coeff2 = 0.5 * coeff3;
    auto &d = data.diffusion_tensor;

    //setting initial forces to zero, setting the initial positions
    for(std::size_t k = 0; k < N; k++){
        
        auto k6 = 6 * k;

        RXI[k] = data.positions[k6];
        RYI[k] = data.positions[k6 + 1];
        RZI[k] = data.positions[k6 + 2];
    };

    for(std::size_t i = 0; i < N-1; i++){
        
        for(std::size_t j = i + 1; j < N; j++){

            Utils::Vector3d pos1 {RXI[i], RYI[i], RZI[i]};
            Utils::Vector3d pos2 {RXI[j], RYI[j], RZI[j]};
            
            //see grid.hpp file for box_geo
            Utils::Vector3d sys_dist = box_geo.get_mi_vector(pos2, pos1);

            RXIJ = sys_dist[0];
            RYIJ = sys_dist[1];
            RZIJ = sys_dist[2];
            
            RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ;
            RIJ = std::sqrt(RIJSQ);
            RIJ3 = RIJ * RIJ * RIJ;

            IC = 3 * i;
            JC = 3 * j;
            IC6 = 3 * N + IC;
            JC6 = 3 * N + JC;

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

            #ifdef OSEEN

            for (std::size_t a = 0; a < 3; a++){

                for (std::size_t b = 0; b < 3; b++){

                    //unit matrix values
                    I = (a == b) ? 1.0 : 0.0;

                    //tt
                    d[IC + a][JC + b] = cons_tr[i] * (coeff1/RIJ) * (I + rr[a][b]);
                    d[JC + b][IC + a] = d[IC + a][JC + b];                

                    //rr
                    d[IC6 + a][JC6 + b] = 0.0;
                    d[JC6 + b][IC6 + a] = 0.0;

                    //rt
                    d[IC6 + a][JC + b] = 0.0;
                    d[JC6 + b][IC + a] = 0.0;

                    //tr
                    d[IC + a][JC6 + b] = 0.0;
                    d[JC + b][IC6 + a] = 0.0;

                }
            }

            #else

            for(std::size_t a = 0; a < 3; a++){

                for(std::size_t b = 0; b < 3; b++){

                    levi_civita[a][b] = 0.0;
                }
            }

            levi_civita[0][1] = ZIJ_hat;
            levi_civita[0][2] = -YIJ_hat;
            levi_civita[1][0] = -ZIJ_hat;
            levi_civita[1][2] = XIJ_hat;
            levi_civita[2][0] = YIJ_hat;
            levi_civita[2][1] = -XIJ_hat;

            for (std::size_t a = 0; a < 3; a++){

                for (std::size_t b = 0; b < 3; b++){

                    //unit matrix values
                    I = (a == b) ? 1.0 : 0.0;

                    //tt
                    d[IC + a][JC + b] = cons_tr[i] * (coeff2/RIJ3 * (I - 3 * rr[a][b]) + (coeff1/RIJ) * (I + rr[a][b]));
                    d[JC + b][IC + a] = d[IC + a][JC + b];                

                    //rr
                    d[IC6 + a][JC6 + b] = -cons_rot[i] * coeff2/RIJ3 * (I - 3 * rr[a][b]);
                    d[JC6 + b][IC6 + a] = d[IC6 + a][JC6 + b];                 

                    //rt
                    d[IC6 + a][JC + b] = levi_civita[a][b] * cons_rot[i] * coeff3/(RIJ * RIJ);
                    d[JC6 + b][IC + a] = d[IC6 + a][JC + b];

                    //tr
                    d[IC + a][JC6 + b] = d[IC6 + a][JC + b];
                    d[JC + b][IC6 + a] = d[IC6 + a][JC + b];

                }
            }

            #endif
        }

    }

    #ifdef OSEEN
    for(std::size_t p = 0; p < N; p++){

        IC = 3 * p;
        IC6 = 3 * N + IC;

        d[IC][IC] = d[IC + 1][IC + 1] = d[IC + 2][IC + 2] = cons_tr[p];
       
    };

    #else
    for(std::size_t p = 0; p < N; p++){

        IC = 3 * p;
        IC6 = 3 * N + IC;

        d[IC][IC] = d[IC + 1][IC + 1] = d[IC + 2][IC + 2] = cons_tr[p];
        d[IC6][IC6] = d[IC6 + 1][IC6 + 1] = d[IC6 + 2][IC6 + 2] = cons_rot[p];
       
    };

    #endif

    dump_diffusion_tensor(d);

}

void calc_velocities(Data_task& data){
    
    const std::size_t N = data.positions.size()/6;
    
    for(std::size_t f = 0; f < N; f++){
        
        auto f6 = 6*f;

        for(std::size_t t = 0; t < 6; t++){

            data.velocities[f6 + t] = 0.0;

            for(std::size_t k = 0; k < 6*N; k++){

                data.velocities[f6 + t] += data.diffusion_tensor[f6 + t][k] * data.forces[k];
            }
        }
    }

    dump_velocity_vector(data.velocities);
}