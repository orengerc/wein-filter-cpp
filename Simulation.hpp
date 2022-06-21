#ifndef NUMERICAL_CPP_SIMULATION_HPP
#define NUMERICAL_CPP_SIMULATION_HPP

#include "Particle.hpp"
#include "ProblemParameters.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <random>
#include <map>

#define PART_C_NUM_PARTICLES 1e5
#define FIELDS_RATIO 3.09e7       //meters per second
#define MAX_VELOCITY 3.17e7       //meters per second
#define MIN_VELOCITY 3.01e7       //meters per second
#define PROTON_MASS 1.67e-27      //kilograms
#define PROTON_CHARGE 1.6e-19     //coulombs


static long double random_double(long double &min, long double &max)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<long double> distr(min, max);
    return distr(eng);
}

static long double get_max(const std::vector<long double> &arr)
{
    long double max = arr[0];
    for (const auto &elem: arr)
    {
        if (elem > max)
        {
            max = elem;
        }
    }
    return max;
}

static long double get_min(const std::vector<long double> &arr)
{
    long double max = arr[0];
    for (const auto &elem: arr)
    {
        if (elem < max)
        {
            max = elem;
        }
    }
    return max;
}

static std::vector<long double> normalize(std::vector<long double> arr){
    auto min = get_min(arr);
    std::vector<long double> ret;
    for(auto elem : arr){
        ret.emplace_back((elem - min));
    }
    return ret;
}

static std::map<double, int> get_histogram(std::vector<long double> data, const long double bin_width)
{
    std::sort(data.begin(), data.end());

    std::map<double, int> histogram;

    double bin = 0; //Choose your starting bin
    for (const auto &e: data)
    {
        e >= bin + bin_width ? bin += bin_width : false;
        ++histogram[bin];
    }

    return histogram;
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
                long double v_min = MIN_VELOCITY;
                long double v_max = MAX_VELOCITY;
                long double r_min = -1 * params->R;
                long double r_max = params->R;
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

    void run(char stage)
    {
        if (stage == 'b')
        {
            for (auto &particle: particles)
            {
                for (Time t = params->dt; t < params->T; t += params->dt)
                {
                    //advance particle
                    particle.advance(t);
                }
            }
        }
        else if (stage == 'c')
        {
            for (auto &particle: particles)
            {
                Time t = params->dt;
                while (true)
                {
                    //advance particle
                    particle.advance(t);

                    //check if particle crashed into the filter
                    if (particle.crashed)
                    {
                        crash_counter++;
                        break;
                    }

                    //check if particle passed the filter
                    if (particle.passed)
                    {
                        break;
                    }

                    t += params->dt;
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

    void print_initial_conditions(bool only_passed = false)
    {
        std::ostringstream oss;
        oss << "initial_conditions.csv";
        std::string filename = oss.str();

        std::ofstream output_csv;
        output_csv.open(filename);

        for (auto &particle: particles)
        {
            bool condition = true;
            if (only_passed) condition = particle.passed;
            if (condition)
            {
                auto y = ((particle.history.begin()->second.r.y()) / params->R);
                auto x = ((MAX_VELOCITY - FIELDS_RATIO) / (particle.history.begin()->second.v.z()));
                output_csv << x << "," << y << "\n";
            }
        }
        output_csv.close();
    }

    void print_one_passed_one_didnt()
    {
        bool done_passed = false;
        bool done_crashed = false;
        for (auto &particle: particles)
        {
            if (particle.passed)
            {
                particle.export_history_to_excel(" passed");
                done_passed = true;
            }
            else if (particle.crashed)
            {
                particle.export_history_to_excel(" crashed");
                done_crashed = true;
            }

            if (done_passed && done_crashed)
            {
                break;
            }
        }
    }

    long double print_passing_percentage() const
    {
        long double pass_percent = (1 - ((long double) crash_counter / PART_C_NUM_PARTICLES)) * 100;
        return pass_percent;
    }

    void print_final_velocity_histogram()
    {
        std::ostringstream oss;
        oss << "final_velocity_histogram.csv";
        std::string filename = oss.str();

        std::ofstream output_csv;
        output_csv.open(filename);
        output_csv << "from,to,num\n";

        std::vector<long double> arr;
        for (auto &particle: particles)
        {
            if (particle.passed)
            {
                auto final_velocity = (--particle.history.end())->second.v.z();
                arr.emplace_back((MAX_VELOCITY - FIELDS_RATIO) / final_velocity);
            }
        }
//        for(auto elem:arr){
//            std::cout << elem << "\n";
//        }
        arr = normalize(arr);

        const int n_bins = 20;
        long double bin_width = get_max(arr) / n_bins;
        auto histogram = get_histogram(arr, bin_width);
        for (const auto& x : histogram)
        {
            output_csv << x.first << "," << x.first + bin_width << "," << x.second << "\n";
        }

        output_csv.close();
    }

};


#endif //NUMERICAL_CPP_SIMULATION_HPP