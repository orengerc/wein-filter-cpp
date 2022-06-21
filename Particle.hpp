#ifndef NUMERICAL_CPP_PARTICLE_HPP
#define NUMERICAL_CPP_PARTICLE_HPP

#include "ProblemParameters.hpp"
#include <unordered_map>
#include <map>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>

typedef long double Time;

enum CTYPE
{
    r,
    v
};

class DoublePair
{
public:
    std::pair<long double, long double> pair;

    DoublePair() : DoublePair(0, 0)
    {}

    DoublePair(long double a, long double b) : pair(std::make_pair(a, b))
    {}

    long double &y()
    {
        return pair.first;
    }

    long double &z()
    {
        return pair.second;
    }

    friend DoublePair &operator+(DoublePair &lhs, const DoublePair &rhs)
    {
        lhs.pair.first += rhs.pair.first;
        lhs.pair.second += rhs.pair.second;
        return lhs;
    }

    friend DoublePair &operator*(DoublePair &lhs, const DoublePair &rhs)
    {
        lhs.pair.first *= rhs.pair.first;
        lhs.pair.second *= rhs.pair.second;
        return lhs;
    }

    friend DoublePair &operator*(DoublePair &lhs, const double &n)
    {
        lhs.pair.first *= n;
        lhs.pair.second *= n;
        return lhs;
    }

    friend DoublePair &operator*(const double &n, DoublePair &rhs)
    {
        rhs.pair.first *= n;
        rhs.pair.second *= n;
        return rhs;
    }
};

class State
{
public:
    DoublePair r;
    DoublePair v;

    friend std::ostream &operator<<(std::ostream &os, State &s)
    {
        std::cout << "(y,z) = (" << s.r.y() << "," << s.r.z() << ")\n";
        std::cout << "(Vy,Vz) = (" << s.v.y() << "," << s.v.z() << ")\n";
    }
};

class Particle
{
public:
    std::map<Time, State> history;
    ProblemParameters *params;
    bool crashed = false;
    bool passed = false;

    Particle(State initial_condition, ProblemParameters *params) : params(params)
    {
        history[0] = initial_condition;
    }

    DoublePair acceleration(DoublePair &v) const
    {
        auto factor = params->q / params->m;
        auto ay = factor * (params->E - (params->B * v.z()));
        auto az = factor * params->B * v.y();
        return {ay, az};
    }

    DoublePair K1(const CTYPE &type, State last_state)
    {
        switch (type)
        {
            case r:
                return params->dt * last_state.v;
            case v:
                auto a = acceleration(last_state.v);
                return params->dt * a;
        }
    }

    DoublePair K2(const CTYPE &type, DoublePair k1, State last_state)
    {
        switch (type)
        {
            case r:
                return params->dt * (last_state.v + (k1 * 0.5));
            case v:
                auto a = k1 * 0.5;
                auto b = last_state.v + a;
                auto acc = acceleration(b);
                return acc * params->dt;
        }
    }

    DoublePair K3(const CTYPE &type, DoublePair k2, State last_state)
    {
        return K2(type, k2, last_state);
    }

    DoublePair K4(const CTYPE &type, DoublePair k3, State last_state)
    {
        switch (type)
        {
            case r:
                return params->dt * (last_state.v + k3);
            case v:
                auto acc = acceleration(k3 + last_state.v);
                return params->dt * acc;
        }
    }

    void advance_taylor(Time &t)
    {
        //take last state
        auto last_state = (--history.end())->second;

        //calc velocity
        auto a = acceleration(last_state.v);
        auto v = last_state.v + a * params->dt;

        //calc position
        auto r = last_state.r + last_state.v * params->dt;

        //insert new data
        history[t] = {r, v};
    }

    void advance_midpoint(Time &t)
    {
        //take last state
        auto last_state = (--history.end())->second;

        //calculate velocity
        auto k1 = K1(v, last_state);
        auto k2 = K2(v, k1, last_state);
        auto v = last_state.v + k2;

        //calculate position
        k1 = K1(r, last_state);
        k2 = K2(r, k1, last_state);
        auto r = last_state.r + k2;

        //insert new data
        history[t] = {r, v};
    }

    void advance_runge(Time &t)
    {
        //take last state
        auto last_state = (--history.end())->second;

        //calculate velocity
        auto k1 = K1(v, last_state);
        auto k2 = K2(v, k1, last_state);
        auto k3 = K3(v, k2, last_state);
        auto k4 = K4(v, k3, last_state);
        auto v = last_state.v + (((long double) 1 / 6) * (k1 + (k2 * 2) + (k3 * 2) + (k4 * 2)));

        //calculate position
        k1 = K1(r, last_state);
        k2 = K2(r, k1, last_state);
        k3 = K3(r, k2, last_state);
        k4 = K4(r, k3, last_state);
        auto r = last_state.r + (((long double) 1 / 6) * (k1 + (k2 * 2) + (k3 * 2) + (k4 * 2)));

        //insert new data
        history[t] = {r, v};
    }

    void advance(Time &t)
    {
        switch (params->method)
        {
            case TAYLOR:
                advance_taylor(t);
                break;
            case MIDPOINT:
                advance_midpoint(t);
                break;
            case RUNGE_KUTTA:
                advance_runge(t);
                break;
        }

        //check if crashed
        auto new_pos = (--history.end())->second;
        if ((std::abs(new_pos.r.y()) >= params->R) && (new_pos.r.z() <= params->L))
        {
            crashed = true;
            return;
        }

        //check if passed
        if ((new_pos.r.z() > params->L))
        {
            passed = true;
        }

    }

    void export_history_to_excel(const std::string &str = "")
    {
        std::vector<std::string> method_names{"TAYLOR", "MIDPOINT", "RUNGE_KUTTA"};
        std::ostringstream oss;
        oss << method_names[params->method] << str << ".csv";
        std::string filename = oss.str();

        std::ofstream output_csv;
        output_csv.open(filename);

        output_csv << "y" << "," << "z" << "," << "vy" << "," << "vz" << "\n";
        for (auto &moment: history)
        {
            output_csv << moment.second.r.y() << "," << moment.second.r.z() << ",";
            output_csv << moment.second.v.y() << "," << moment.second.v.z() << "\n";
        }
        output_csv.close();
    }
};


#endif //NUMERICAL_CPP_PARTICLE_HPP
