#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "ipotential.hpp"
#include <random>
#include <thread>
#include <QVector>
#include <QFile>
#include <QString>
#include <QTextStream>

struct SimParameters
{
    u_int64_t particles_amount;
    double cell_sizes[3];
    bool cell_bounds[3];
    double T;
    double m;
    double k;
    double dt;
    double runtime;
    size_t steps;
    int ansamble;  // nve = 1, nvt = 2, npt = 3
    bool start_with_random;
    QString state_file;
};

struct SysParameters
{
    double E_k;
    double T;
    double V;
    double density;
    double P;
};

class Simulation
{
public:
    Simulation(QString config_file, IPotential &U);
    ~Simulation();
    void Run();
    QVector<double> RDF(int size = 100);
    QVector<double> MSD();

private:
    inline void InitStateSet(bool not_ready[9]);  //Метод для вызова в отдельном потоке, параллельно с чтением конфиг файла

    inline void RandomStateSet();
    inline void OrderStateSet();
    inline void LoadStateSet();

    inline void CoordCalc();
    inline void StateCalc();
    inline void ForcesCalc();
    inline void AccelCalc();
    inline void VelocCalc();

    inline void BoundConitApply();
    inline void PeriodicConditApply(const int &axis);
    inline void IsolatedConditApply(const int &axis);

    inline void StatAnsApply();
    inline void BarostatApply();
    inline void ThermostatApply();

    inline void SystemStateWrite();

    SimParameters params_;
    IPotential &U;
    double * cache_;
    size_t cache_space_;
    double * coord_[3];
    double * vel_[3];
    double * accel_[3];
    double ** pot_table_;
    double ** force_table_;
    SysParameters state_;
    QVector<SysParameters> dynamic_;

};

#endif // SIMULATION_HPP
