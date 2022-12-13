#include "lj.hpp"

LJ::LJ(double e, double q) : e(e), q(q) {}

double LJ::operator()(double r, double r_cut)
{
    if(r_cut == 0 || r < r_cut)
        return 4 * e * (qPow(q/r, 12) - qPow(q/r, 6));
    return 0;
}
