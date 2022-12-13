#ifndef LJ_H
#define LJ_H

#include <QtMath>
#include "ipotential.hpp"

class LJ : public IPotential
{
public:
    LJ(double e, double q);
    inline virtual double operator ()(double r, double r_cut);
private:
    double e;
    double q;
};

#endif // LJ_H
