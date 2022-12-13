#ifndef IPOTENTIAL_H
#define IPOTENTIAL_H


class IPotential
{
public:
    virtual double operator()(double r, double r_cut) = 0;
    virtual ~IPotential() {};
};

#endif // IPOTENTIAL_H
