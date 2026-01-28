#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <mutex>
#include <atomic>
#include <thread>
#include <functional>
#include <filesystem>
#include <ctime>
#include <utility>

#include "./source/project.hpp"
#include "./source/constants.hpp"

using namespace std;
using namespace std::chrono;
using namespace project;

int main(int argc, char **argv)
{   
    // Constants for non-relativistic white dwarfs
    const long double mu_e_wd_nonR = 2.0L;
    const long double K_wd_nonR = 1.0L / 20.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * ((__planck * __planck) / __mass_electron) * powl((1.0L / (mu_e_wd_nonR * __mass_proton)), (5.0L / 3.0L));

    // Solve lane emden for n in (-2,10)
    freopen("./output/laneEmdenSolutions.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::solve_LE_over_n(0.01L, -2.0L, 10.0L, 0.5L, 20.0L);

    // Mass-Radius for White Dwarf using LE
    freopen("./output/massRadius1_5.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::MR_LE_over_density_logspace(1.5L, 0.01L, K_wd_nonR, 1000000, 1000000000000, 1000);

    // Mass-Radius for White Dwarf using TOV
    freopen("./output/massRadius1_5.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::MR_LE_over_density_logspace(1.5L, 0.01L, K_wd_nonR, 1000000, 1000000000000, 1000);

    // Get maximum thread count
    // unsigned int max_threads = thread::hardware_concurrency();
    // cout << max_threads << " threads available." << endl << setprecision(std::numeric_limits<long double>::digits10 + 1);

    return 0;
}