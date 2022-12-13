#include <iostream>
#include <cmath>
#include <random>

#include "mytimer.hpp"

#include "simulation.hpp"
#include "lj.hpp"

/*double U(double r, double e, double o, double r_cut)
{
    using namespace std;
    if(r < r_cut)
        return 4 * e * (pow(o/r, 12) - pow(o/r, 6));
    return 0;
}

enum class quantities
{
    SI,
    LJ
};

enum class potential
{
    LJ,
    Dz,
    SW,
    Uk
};*/


int main()
{
    MyTimer timer;
    IPotential &U = *(new LJ(1, 1));
    std::cout << " time: " << timer.fromStart() << "\n";
    Simulation MD("config.txt", U);
    std::cout << " time: " << timer.fromLast() << "\n";
    MD.Run();

    return 0;
}
