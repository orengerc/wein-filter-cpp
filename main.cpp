#include "Simulation.hpp"
#include <string>
#include <set>

bool s_equals(const std::string &a, const std::string &b)
{
    return std::equal(a.begin(), a.end(),
                      b.begin(), b.end(),
                      [](char a, char b) {
                          return tolower(a) == tolower(b);
                      });
}

long double average(std::vector<long double> const &v)
{
    if (v.empty())
    {
        return 0;
    }

    auto const count = static_cast<long double>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}

long double distance(DoublePair x, DoublePair y)
{
    // Calculating distance
    auto x1 = x.pair.first;
    auto y1 = x.pair.second;
    auto x2 = y.pair.first;
    auto y2 = y.pair.second;
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

int main(int argc, char **argv)
{
    if (s_equals(argv[1], "b"))
    {
        for (const auto &method: std::set{TAYLOR, MIDPOINT, RUNGE_KUTTA})
        {
            std::map<long double, long double> errors;
            for (double dt = 0.01; dt > 1e-10; dt /= 2)
            {
                ProblemParameters params{};
                params.method = method;
                params.dt = dt;
                Simulation partB(1, &params);
                partB.run('b');
                auto final_pos = (--partB.particles[0].history.end())->second.r;
                auto analytical_solution = DoublePair(0, 2 * M_PI);
                std::cout << distance(final_pos, analytical_solution) << "\n";
                errors[dt] = distance(final_pos, analytical_solution);
                partB.export_to_excel();
            }
        }
    }

    else if (s_equals(argv[1], "c"))
    {
        ProblemParameters params{FIELDS_RATIO, 1, PROTON_MASS, PROTON_CHARGE, 1e-9, 0.003, 1,
                                 RUNGE_KUTTA};
        Simulation partC(PART_C_NUM_PARTICLES, &params, true);
        partC.run('c');
        partC.print_one_passed_one_didnt();
        partC.print_initial_conditions(true);
        partC.print_final_velocity_histogram();
        std::cout << "done\n";

        {
            std::vector<long double> percentages;
            for (int i = 0; i < 53; ++i)
            {
                ProblemParameters parameters{FIELDS_RATIO, 1, PROTON_MASS, PROTON_CHARGE, 1e-9, 0.003, 1,
                                             RUNGE_KUTTA};
                Simulation sim(PART_C_NUM_PARTICLES, &parameters, true);
                sim.run('c');
                percentages.emplace_back(sim.print_passing_percentage());
            }
            for (auto elem: percentages)
            {
                std::cout << elem << "\n";
            }
            std::cout << "average pass percentage for 20 runs of the simulation: " << average(percentages) << " %\n";
        }
    }

    return 0;
}
