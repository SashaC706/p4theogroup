#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>

#include "./constants.hpp"

using namespace std;
using namespace project;

namespace project
{
    pair<vector<long double>, vector<long double>> solve_TOV(long double h, long double n, long double P_c)
    {
        bool solved = false;
        vector<long double> _r = {h}; // Avoid any possible singularity CHECK IF REQUIRED IN TOV
        vector<long double> _P = {P_c}; // Initial condition
        long double last_dP = -h / 3.0L;

        long double M = 0; // Initial condition

        while (!solved)
        {
            long double r = _r.back();
            long double P = _P.back();

            _r.push_back(r + h);

            long double k1, k2, k3, k4;
            long double dk1, dk2, dk3, dk4;
            
            if (r > 100) throw 1;

            long double u = 0;

            if (u <= 0 || isnan(u))
            {
                solved = true;
                _r.pop_back();
            }
            else
            {
                _P.push_back(u);
            }
        }

        return {_r, _P};
    }
    void MR_TOV_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs)
    {
        vector<long double> _xi;
        vector<long double> _theta;

        tie(_xi, _theta) = project::solve_TOV(h, n, 0);

        vector<long double> rho_c;

        for (int i = 0; i < density_runs; i++)
        {
            rho_c.push_back(rho_min * powl(rho_max / rho_min, static_cast<long double>(i) / (density_runs - 1)));
        }

        for (long double rho : rho_c)
        {
            long double alpha = sqrtl((n + 1.0L) * K * powl(rho, (1.0L / n - 1)) / (4.0L * __pi * __gravitation));
            long double R = alpha * _xi.back();
            long double M = 0.0L;
            for (int j = 1; j < _xi.size(); j++)
            {
                long double r = alpha * _xi[j];
                long double dr = alpha * (_xi[j] - _xi[j - 1]);
                long double rho_r = rho * powl(_theta[j], n);

                M += 4.0L * __pi * (r * r) * rho_r * dr;
            }

            cout << M << "," << R << endl;
        }
    }
}