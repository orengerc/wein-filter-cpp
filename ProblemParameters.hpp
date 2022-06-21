#ifndef NUMERICAL_CPP_PROBLEMPARAMETERS_HPP
#define NUMERICAL_CPP_PROBLEMPARAMETERS_HPP

#include <string>
#include <iostream>
#include "cmath"
#include <map>

#define NUM_PERIODS 1
#define DEFAULT_DT 0.01

enum METHOD
{
    TAYLOR,
    MIDPOINT,
    RUNGE_KUTTA
};

class ProblemParameters
{
public:
    long double E;
    long double B;
    long double m;
    long double q;
    long double w;
    long double T;
    long double dt;
    long double R;
    long double L;
    METHOD method;

    ProblemParameters() : ProblemParameters(1, 1, 1, 1, DEFAULT_DT, 1, 1, TAYLOR)
    {}

    ProblemParameters(long double E, long double B, long double m, long double q, long double dt, long double R, long double L,
                      METHOD method) : E(E), B(B), m(m), q(q), dt(dt), R(R), L(L),
                                       method(method), w(q * B / m), T(NUM_PERIODS * 2 * M_PI / w)
    {}

    ProblemParameters(const ProblemParameters &other) = default;
};


#endif //NUMERICAL_CPP_PROBLEMPARAMETERS_HPP
