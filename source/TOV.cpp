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
    // POSSIBLY BREAK THESE OUT INTO SEPARATE FILES AND HAVE THEM BE PARAMETERS IN FUNCTION CALL
    long double rho_PeS(long double P, long double n, long double K) { return powl(P/K,n / (n+1.0L)); }
    long double mass_continuity(long double r, long double rho)
    {
        return __4pi * r * r * rho;
    };
    long double tolman_oppenheimer_volkoff(long double r, long double P, long double M, long double rho)
    {
        long double term1 = (rho + P / __c2);
        long double term2 = (M + 4.0L * __pi * r * r * r * P / __c2);
        long double term3 = r * (r - 2.0L * __gravitation * M / __c2);
        return -(__gravitation * term1 * term2) / term3;
    };

    pair<vector<long double>, vector<long double>> solve_TOV(long double h, long double n, long double P_c, long double K)
    {
        // Boolean turns true when pressure drops below 0 and this determines the radial extent of star
        bool solved = false;

        // initial conditions
        // Radius vector starting from initial step size to avoid singularity CHECK IF REQUIRED IN TOV
        vector<long double> _r = {h};
        vector<long double> _P = {P_c};
        vector<long double> _M = {0.0L};

        while (!solved)
        {
            // Grab radius and pressure for THIS step
            long double r = _r.back();
            long double P = _P.back();

            _r.push_back(r + h); // Append radius of NEXT step THIS MAY BE INCORRECT, THINK r DECLARATION SHOULD FOLLOW THIS MATHEMATICALLY

            // Calculate density from EOS
            long double rho = powl((P / K), (n / (n + 1.0L)));

            // Initialise Runge-Kutta variables
            long double M_k1, M_k2, M_k3, M_k4;
            long double P_k1, P_k2, P_k3, P_k4;

            // Perform RK4 for mass
            M_k1 = mass_continuity(r         , rho) * h;
            M_k2 = mass_continuity(r + 0.5L*h, rho) * h;
            M_k3 = mass_continuity(r + 0.5L*h, rho) * h;
            M_k4 = mass_continuity(r + h     , rho) * h;

            // POSSIBLE IMPROVEMENT USING SUBSTEPS OF PRESSURE ASWELL INSTEAD OF USING DENSITY
            // M_k1 = __4pi * r * r                    * rho_PeS(P, n, K);
            // M_k2 = __4pi * powl((r + h/2.0L), 2.0L) * rho_PeS(P, n, K);
            // M_k3 = __4pi * powl((r + h/2.0L), 2.0L) * rho_PeS(P, n, K);
            // M_k4 = __4pi * powl((r + h), 2.0L)      * rho_PeS(P, n, K);

            // Update mass vector
            _M.push_back(_M.back() + (1.0L/6.0L) * (M_k1 + 2.0L*M_k2 + 2.0L*M_k3 + M_k4));
            
            // Perform RK4 for pressure
            P_k1 = tolman_oppenheimer_volkoff(r             , P               , _M.back(), rho) * h;
            P_k2 = tolman_oppenheimer_volkoff(r + 0.5L * h  , P + 0.5L * P_k1 , _M.back(), rho) * h;
            P_k3 = tolman_oppenheimer_volkoff(r + 0.5L * h  , P + 0.5L * P_k2 , _M.back(), rho) * h;
            P_k4 = tolman_oppenheimer_volkoff(r + h         , P + P_k3        , _M.back(), rho) * h;

            // Update pressure vector
            _P.push_back(P + (1.0L/6.0L) * (P_k1 + 2.0L * P_k2 + 2.0L * P_k3 + P_k4));

            // If pressure is negative then set to solved and remove next r value MAY ALSO REMOVE NEGATIVE PRESSURE VALUE
            if (_P.back() <= 0 || isnan(_P.back()))
            {
                solved = true;
                _r.pop_back();
            }
        }

        return {_r, _P};
    }
    // void MR_TOV_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs)
    // {
    //     vector<long double> _xi;
    //     vector<long double> _theta;

    //     tie(_xi, _theta) = project::solve_TOV(h, n, 0);

    //     vector<long double> rho_c;

    //     for (int i = 0; i < density_runs; i++)
    //     {
    //         rho_c.push_back(rho_min * powl(rho_max / rho_min, static_cast<long double>(i) / (density_runs - 1)));
    //     }

    //     for (long double rho : rho_c)
    //     {
    //         long double alpha = sqrtl((n + 1.0L) * K * powl(rho, (1.0L / n - 1)) / (4.0L * __pi * __gravitation));
    //         long double R = alpha * _xi.back();
    //         long double M = 0.0L;
    //         for (int j = 1; j < _xi.size(); j++)
    //         {
    //             long double r = alpha * _xi[j];
    //             long double dr = alpha * (_xi[j] - _xi[j - 1]);
    //             long double rho_r = rho * powl(_theta[j], n);

    //             M += 4.0L * __pi * (r * r) * rho_r * dr;
    //         }

    //         cout << M << "," << R << endl;
    //     }
    // }
}