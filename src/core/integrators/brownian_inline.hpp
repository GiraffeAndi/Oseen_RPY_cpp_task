/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INTEGRATORS_BROWNIAN_INLINE_HPP
#define INTEGRATORS_BROWNIAN_INLINE_HPP

#include "config.hpp"

#include "ParticleRange.hpp"
#include "integrate.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"
#include "thermostats/brownian_inline.hpp"

#ifdef Oseen_RPY

#include "Oseen_RPY.hpp"

#endif

#include <utils/math/sqr.hpp>

inline void brownian_dynamics_propagator(BrownianThermostat const &brownian,
                                         const ParticleRange &particles,
                                         double time_step, double kT, double sigma) {

#ifdef Oseen_RPY
    
    Data_task data;

    data.pot_energ = 0.0;

    data.virial = 0.0;

    //dummy value, has to be chosen correctly
    data.dt = 0.1;

    //value chosen by me maybe needed to be adjusted see manual or user-defined?
    double rcut = 3.2;

    //dummy values for consii and consij for now

    double consii = 1.0;

    double consij = 1.0;

    //consii and consij have to be defined still

    int N = particles.size();

    //randomly chosen by me
    double temp = 32.0;


    //resize to avoid errors

    data.positions.resize(N, std::vector<double>(3, 0.0));

    data.forces.resize(N, std::vector<double>(3, 0.0));

    data.diffusion_tensor.resize(3*N, std::vector<double>(3*N, 0.0));

    data.CRND.resize(3 * N, 0.0);


    for(int i=0; i < N; i++){

        //x-coordinates
        data.positions[i][0] = particles[i].pos()[0];
        
        //y-coordinates
        data.positions[i][1] = particles[i].pos()[1];

        //z-coordinates
        data.positions[i][2] = particles[i].pos()[2];
    }
    
    calc_forces(data, sigma, rcut, consii, consij);

    covar(data, dt);

    move(data, dt, temp);

    //updating pos for next step
    for (auto &p : particles){

            if (!p.is_virtual() or thermo_virtual){
                
                for(int i = 0; i < N; i++){
                    particles[i].pos()[0] = data.positions[i][0];

                    particles[i].pos()[1] = data.positions[i][1];

                    particles[i].pos()[2] = data.positions[i][2];

                };

#ifdef ROTATION

      if (!p.can_rotate())
        continue;
      convert_torque_to_body_frame_apply_fix(p);
      p.quat() = bd_drag_rot(brownian.gamma_rotation, p, time_step);
      p.omega() = bd_drag_vel_rot(brownian.gamma_rotation, p);
      p.quat() = bd_random_walk_rot(brownian, p, time_step, kT);
      p.omega() += bd_random_walk_vel_rot(brownian, p);
#endif // ROTATION
    }
  }
  increment_sim_time(time_step);


#endif // INTEGRATORS_BROWNIAN_INLINE_HPP
    
};

#else



  for (auto &p : particles) {
    // Don't propagate translational degrees of freedom of vs
    if (!p.is_virtual() or thermo_virtual) {
      p.pos() += bd_drag(brownian.gamma, p, time_step);
      p.v() = bd_drag_vel(brownian.gamma, p);
      p.pos() += bd_random_walk(brownian, p, time_step, kT);
      p.v() += bd_random_walk_vel(brownian, p);
#ifdef ROTATION
      if (!p.can_rotate())
        continue;
      convert_torque_to_body_frame_apply_fix(p);
      p.quat() = bd_drag_rot(brownian.gamma_rotation, p, time_step);
      p.omega() = bd_drag_vel_rot(brownian.gamma_rotation, p);
      p.quat() = bd_random_walk_rot(brownian, p, time_step, kT);
      p.omega() += bd_random_walk_vel_rot(brownian, p);
#endif // ROTATION
    }
  }

#endif

  increment_sim_time(time_step);
}

#endif // INTEGRATORS_BROWNIAN_INLINE_HPP
