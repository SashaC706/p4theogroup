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

using namespace std;
using namespace std::chrono;

constexpr unsigned long long int million = 1000000ULL;
constexpr unsigned long long int billion = 1000000000ULL;
constexpr unsigned long long int trillion = 1000000000000ULL;

class Estimator
{
public:
    Estimator(long double initial_guess, unsigned long long int given_threads)
        : _2lt(initial_guess / 2.0L),
          _r((initial_guess / 4.0L) / 2.0L),
          _thread_count(given_threads),
          _start(steady_clock::now())
    {}
    // void relative() {cout << "Relative Err:  " << std::scientific << (pi - M_PI) / M_PI << endl;}
    void run()
    {
        // Print details on the log
        _print();

        // Initialise workers
        for (unsigned int t = 0; t < _thread_count; t++)
        {
            _threads.emplace_back(&Estimator::_worker, this);
        }

        // Join all threads, Should not run, Partially deprecated
        for (auto &t : _threads)
        {
            t.join();
        }
    }
    void rejoin(unsigned long long int N, unsigned long long int n)
    {
        _N += N;   // Recombine itteration totals
        _n += n;   // Recombine intersection totals
        _R++;      // Itterate thread run count

        _printr(); // Output state given by backup
    }
private:
    // Constants
    const unsigned long long int _batch_width = billion; // How many cycles each worker will run before resyncing
    const unsigned long long int _thread_count;          // Number of threads simulation runs on
    const long double _2lt;                              // Separation/Length ration multiplied by 2
    const long double _r;                                // Radius of match stick, avoids recomputation
    const steady_clock::time_point _start;               // Start time point
    // Random number generation
    random_device _rd; // Obtain a random number for seeding random number generator
    // Multithreadding
    atomic<bool> _stop_flag{false}; // To allow for watcher thread to halt workers
    vector<thread> _threads;        // Vector to hold all threads
    // Itterators
    atomic<unsigned long long int> _N{0ULL}; // Number of itterations
    atomic<unsigned long long int> _n{0ULL}; // Number of counted intersections
    atomic<unsigned long long int> _R{0ULL}; // Number of thread runs
    // Recombination function
    void _rejoin(unsigned long long int N, unsigned long long int n)
    {
        // lock_guard<mutex> lock(_mtx); // Ensure no conflict between threads DEPRECATED
        _N += N;                      // Recombine itteration totals
        _n += n;                      // Recombine intersection totals
        _R++;                         // Itterate thread run count

        // if (_R != 0 && _R % 10 == 0)
        // {
            _printr();               // Output update
        // }

        // if (_R.load() >= 24ULL) // Testing condition for timing improvements
        // {
        //     _stop_flag.store(true); // Halt all threads
        // }
    }
    void _worker()
    {
        mt19937_64 thread_gen(_rd() + hash<thread::id>{}(this_thread::get_id())); // Seed 64-bit standard mersenne twister generator with combination of seed and thread id
        uniform_real_distribution<long double> thread_x_distribution(0.0L, 0.5L); // Uniform displacement distribution, symmetric around the line so 0 -> t/2 is the same as -t/2 -> t/2, with use of length ratio t = 1
        uniform_real_distribution<long double> thread_v_distribution(0.0L, 1.0L); // Uniform random distribution for vectors
        unsigned long long int thread_N = 0ULL;                                            // Thread specific count of itterations
        unsigned long long int thread_n = 0ULL;                                            // Thread specific count of intersections

        while (!_stop_flag.load())
        {
            if (thread_N == _batch_width)
            {
                // Recombine then reset totals
                _rejoin(thread_N, thread_n);
                thread_N = 0ULL;
                thread_n = 0ULL;
            }

            long double x = thread_x_distribution(thread_gen);   // Random x displacement
            long double u_x = thread_v_distribution(thread_gen); // Random vector component
            long double u_y = thread_v_distribution(thread_gen); // Random vector component

            long double norm = sqrtl(u_x * u_x + u_y * u_y);

            while (norm > 1.0L || norm < 0.001L) // Reject vector if it lies outside unit circle and if it is problematically small
            {
                // Generate replacement vector
                u_x = thread_v_distribution(thread_gen);
                u_y = thread_v_distribution(thread_gen);
                norm = sqrtl(u_x * u_x + u_y * u_y);
            }

            thread_N++; // Itterate this threads itteration count
            if (x <= (long double)(_r * u_x / norm))
            {
                thread_n++; // Intersection thus itterate counter
            }
        }
    }
    void _print()
    {
        cout<< setprecision(1) << fixed
            << "time_taken"
            << ":" << "N_itterations"
            << ":" << "n_hits"
            << ":" << "2*l/t"
            << ":" << "pi_estimate"
            << endl;
    }
    void _printr()
    {
        long double P = (long double)(_n.load()) / (long double)(_N.load());
        long double pi = _2lt / (P);

        cout<< setprecision(1) << fixed
            << duration<double, std::milli>(steady_clock::now() - _start).count() / 1000
            << ":" << _N.load()
            << ":" << _n.load()
            << ":" << _2lt
            << ":" << setprecision(std::numeric_limits<long double>::digits10 + 1) << pi
            << endl;
    }
};

int main(int argc, char **argv)
{
    // Reroute console to file
    freopen("./.log", "w", stdout);

    // Initialise constants
    constexpr long double initial_pi = 3.0L; // Inital guess for pi, needed for recursive efficiency code

    // Get maximum thread count
    unsigned int max_threads = thread::hardware_concurrency();
    cout << "Running on " << max_threads << " threads." << endl;

    // Create simulation from class
    Estimator u(initial_pi, max_threads);

    // Check if partial run exists and restart from the runs log file
    // if (!fs::exists("./.log"))
    // {
    //     unsigned long long int N;
    //     unsigned long long int n;

    //     // unimplemented

    //     // Input backup counts to new simulation
    //     u.rejoin(N, n);
    // }

    // Run simulation
    u.run();

    return 0;
}