#ifndef NUMERICAL_CPP_SIMULATION_H
#define NUMERICAL_CPP_SIMULATION_H

#include "Particle.h"
#include "ProblemParameters.h"
#include <vector>
#include <iostream>
#include <cassert>
#include <random>


#define PART_C_NUM_PARTICLES 1e6


static double random_double(double &min, double &max)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(min, max);
    return distr(eng);
}

class Simulation
{
public:
    std::vector<Particle> particles;
    ProblemParameters *params;
    unsigned long crash_counter;

    Simulation(int n_particles, ProblemParameters *params, bool random = false) : params(params), crash_counter(0)
    {
        //initialize particles
        for (int i = 0; i < n_particles; i++)
        {
            State initial_condition;
            if (random)
            {
                //select random initial condition
                double v_min = 0;  //todo change
                double v_max = 0.5;  //todo change
                double r_min = -1 * params->R;
                double r_max = params->R;
                initial_condition = {{random_double(r_min, r_max), 0},
                                     {0,                           random_double(v_min, v_max)}};
            }
            else
            {
                initial_condition = {{0, 0},
                                     {0, 3 * (params->E / params->B)}};
            }
            particles.emplace_back(Particle(initial_condition, params));
        }
    }

    void run()
    {
        for (auto &particle: particles)
        {
            for (Time t = params->dt; t < params->T; t += params->dt)
            {
                //advance particle
                particle.advance(t);

                //check if particle crashed into the filter
                if (particle.crashed)
                {
//                    std::cout << (--particle.history.end())->second;
                    crash_counter++;
                    break;
                }
            }
        }
    }

    void export_to_excel()
    {
        for (auto &particle: particles)
        {
            particle.export_history_to_excel();
        }
    }

    void analyze() const
    {
        std::cout << "\npercentage of survivors: " << (1 - ((long double) crash_counter / PART_C_NUM_PARTICLES)) << "\n";
    }

    void plot_final_velocity_histogram() const
    {

    }

    void get_passing_percentage() const
    {

    }
};


#endif //NUMERICAL_CPP_SIMULATION_H