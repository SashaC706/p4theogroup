#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>

using namespace std;

namespace project
{
    pair<vector<long double>, vector<long double>> solve_LE(long double h, long double n);
    void solve_LE_over_n(long double h, long double n_min, long double n_max, long double n_h, long double end_point);
    void MR_LE_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs);
    pair<vector<long double>, vector<long double>> solve_TOV(long double h, long double n);
    void MR_TOV_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs);
}