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
    pair<vector<long double>, vector<long double>> solve_LE(long double h, long double n)
    {
        bool solved = false;
        vector<long double> xi = {h};
        vector<long double> theta = {1.0L - (h * h) / 6.0L}; // Avoid singularity by using xi0=h
        long double last_dtheta = -h / 3.0L;

        while (!solved)
        {
            xi.push_back(xi.back() + h);

            long double k1, k2, k3, k4;
            long double dk1, dk2, dk3, dk4;

            k1 = h * last_dtheta;
            dk1 = h * (-(2.0L * last_dtheta) / (xi.back())-powl(theta.back(), n));

            k2 = h * (last_dtheta + dk1 / 2.0L);
            dk2 = h * (-(2.0L * (last_dtheta + dk1 / 2.0L)) / (xi.back() + h / 2.0L) - powl((theta.back() + k1 / 2.0L), n));

            k3 = h * (last_dtheta + dk2 / 2.0L);
            dk3 = h * (-(2.0L * (last_dtheta + dk2 / 2.0L)) / (xi.back() + h / 2.0L) - powl((theta.back() + k2 / 2.0L), n));

            k4 = h * (last_dtheta + dk3);
            dk4 = h * (-(2.0L * (last_dtheta + dk3)) / (xi.back() + h) - powl((theta.back() + k3), n));

            long double u = theta.back() + (1.0L / 6.0L * k1) + (1.0L / 3.0L * k2) + (1.0L / 3.0L * k3) + (1.0L / 6.0L * k4);
            last_dtheta = last_dtheta + (1.0L / 6.0L * dk1) + (1.0L / 3.0L * dk2) + (1.0L / 3.0L * dk3) + (1.0L / 6.0L * dk4);

            if (xi.back() > 100) throw 1;

            if (u <= 0 || isnan(u))
            {
                solved = true;
                xi.pop_back();
            }
            else
            {
                theta.push_back(u);
            }
        }

        return {xi, theta};
    }
    void solve_LE_over_n(long double h, long double n_min, long double n_max, long double n_h, long double end_point)
    {
        vector<string> csv;
        vector<long double> xi;
        long double xi0 = h;

        for (long double last_xi = xi0; last_xi <= end_point; last_xi += h)
        {
            xi.push_back(last_xi);
            csv.push_back(to_string(last_xi));
        }
        for (long double i = n_min; i <= n_max; i += n_h)
        {
            long double last_theta = 1.0L - (xi0 * xi0) / 6.0L;
            long double last_dtheta = -xi0 / 3.0L;
            
            csv[0] += "," + to_string(last_theta);

            for (long double j = 1; j < xi.size(); j++)
            {
                long double k1, k2, k3, k4;
                long double dk1, dk2, dk3, dk4;

                k1 = h * last_dtheta;
                dk1 = h * (-(2.0L * last_dtheta) / (xi[j])-powl(last_theta, i));

                k2 = h * (last_dtheta + dk1 / 2.0L);
                dk2 = h * (-(2.0L * (last_dtheta + dk1 / 2.0L)) / (xi[j] + h / 2.0L) - powl((last_theta + k1 / 2.0L), i));

                k3 = h * (last_dtheta + dk2 / 2.0L);
                dk3 = h * (-(2.0L * (last_dtheta + dk2 / 2.0L)) / (xi[j] + h / 2.0L) - powl((last_theta + k2 / 2.0L), i));

                k4 = h * (last_dtheta + dk3);
                dk4 = h * (-(2.0L * (last_dtheta + dk3)) / (xi[j] + h) - powl((last_theta + k3), i));

                last_theta = last_theta + (1.0L / 6.0L * k1) + (1.0L / 3.0L * k2) + (1.0L / 3.0L * k3) + (1.0L / 6.0L * k4);
                last_dtheta = last_dtheta + (1.0L / 6.0L * dk1) + (1.0L / 3.0L * dk2) + (1.0L / 3.0L * dk3) + (1.0L / 6.0L * dk4);

                csv[j] += "," + to_string(last_theta);
            }
        }
        for (auto& entry : csv)
        {
            cout << entry << endl;
        }
    }
    void MR_LE_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs)
    {
        vector<long double> _xi;
        vector<long double> _theta;

        tie(_xi, _theta) = project::solve_LE(h, n);

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