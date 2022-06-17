#include "Simulation.h"
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

int main(int argc, char **argv)
{
    if (s_equals(argv[1], "b"))
    {
        ProblemParameters params{};
        for (const auto &method: std::set{TAYLOR, MIDPOINT, RUNGE_KUTTA})
        {
            params.method = method;
            Simulation partB(1, &params);
            partB.run();
            partB.export_to_excel();
        }
    }

    else if (s_equals(argv[1], "c"))
    {
        ProblemParameters params{1, 1, 1, 1, 1, 0.003, 1, MIDPOINT};  //todo change to RUNGE when it works
        Simulation partC(1e5, &params, true);
        partC.run();
        partC.analyze();
    }

    return 0;
}
