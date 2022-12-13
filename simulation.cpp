#include "simulation.hpp"
#include <iostream>

Simulation::Simulation(QString config_file, IPotential &U) : U(U)
{
    params_.k = 1.38e-23;
    QFile inputFile(config_file);
    if (inputFile.open(QIODevice::ReadOnly))
    {

       QTextStream in(&inputFile);
       bool searching[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
       std::thread init_state_set_trd(&Simulation::InitStateSet, this, searching);
       while (!in.atEnd())
       {
          QString line = in.readLine();
          QStringList content = line.split(" ");
          if(searching[0] && line.contains("particles_amount:"))
          {
              params_.particles_amount = content[1].toInt();
              searching[0] = false;
          }
          if(searching[1] && line.contains("cell_sizes:"))
          {
              searching[1] = false;
              params_.cell_sizes [0] = content[1].toDouble();
              params_.cell_sizes [1] = content[2].toDouble();
              params_.cell_sizes [2] = content[3].toDouble();
          }
          if(searching[2] && line.contains("cell_bounds:"))
          {
              searching[2] = false;
              params_.cell_bounds[0] = (bool)content[1].toInt();
              params_.cell_bounds[1] = (bool)content[2].toInt();
              params_.cell_bounds[2] = (bool)content[3].toInt();
          }
          if(searching[3] && line.contains("T:"))
          {
              searching[3] = false;
              params_.T = content[1].toDouble();
          }
          if(searching[4] && line.contains("m:"))
          {
              searching[4] = false;
              params_.m = content[1].toDouble();
          }
          if(searching[5] && line.contains("timestep:"))
          {
              searching[5] = false;
              params_.dt = content[1].toDouble();
          }
          if(searching[6] && line.contains("runtime:"))
          {
              searching[6] = false;
              params_.runtime = content[1].toDouble();
          }
          if(searching[7] && line.contains("ansamble:"))
          {
              searching[7] = false;
              if(content[1] == "NVE")
              {
                  params_.ansamble = 1;
              }
              else if(content[1] == "NVT")
              {
                  params_.ansamble = 2;
              }
              else if(content[1] == "NPT")
              {
                  params_.ansamble = 3;
              }
              else
              {
                  params_.ansamble = 0;
              }
          }
          if(searching[8] && line.contains("start_conf:"))
          {
              searching[8] = false;
              if(content[1] == "random")
              {
                  params_.start_with_random = true;
              }
              else
              {
                  params_.start_with_random = false;
                  params_.state_file = content[2];
              }
              //Нужна реализация для загрузки конфигурации из файла
          }
       }
       inputFile.close();
       params_.steps = static_cast<int>(params_.runtime / params_.dt);
       init_state_set_trd.join();
    }
    else
    {
        std::cout << "Error while file opening";
    }
}

Simulation::~Simulation()
{
    delete[] cache_;
}

void Simulation::InitStateSet(bool *not_ready)
{
    bool not_done[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    while(not_done[0] && not_done[8])
    {
        // Если количество частиц уе известно, выделяем память и раскидываем указатели
        if(not_done[0] && !not_ready[0])
        {
            cache_space_ = 9 * params_.particles_amount;
            cache_ = new double[cache_space_];
            for(int i = 0; i < 3; ++i)
            {
                coord_[i] = cache_ + i * params_.particles_amount;
                vel_[i] = cache_ + (3 + i) * params_.particles_amount;
                accel_[i] = cache_ + (6 + i) * params_.particles_amount;
            }
            not_done[0] = false;
        }
        // Если известно правило устновки начально структуры, устанавливаем её
        if(not_done[8] && !not_ready[8])
        {
            if(params_.start_with_random)
            {
                RandomStateSet();
            }
        }
    }
}

void Simulation::RandomStateSet()
{
    using namespace std;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> x_coord(0, params_.cell_sizes[0]);
    uniform_real_distribution<> y_coord(0, params_.cell_sizes[1]);
    uniform_real_distribution<> z_coord(0, params_.cell_sizes[2]);
    double vel_avearge = pow(2 * params_.k * params_.T / params_.m, 0.5) / 3;
    uniform_real_distribution<> vel(0.5 * vel_avearge, 1.5 * vel_avearge);
    for(size_t i = 0; i < params_.particles_amount; ++i)
    {
        coord_[0][i] = x_coord(gen);
        coord_[1][i] = y_coord(gen);
        coord_[2][i] = z_coord(gen);
        vel_[0][i] = vel(gen);
        vel_[1][i] = vel(gen);
        vel_[2][i] = vel(gen);
        accel_[0][i] = 0;
        accel_[1][i] = 0;
        accel_[2][i] = 0;
    }
}

void Simulation::Run()
{
    for(size_t step = 0; step < params_.steps; ++step)
    {
        CoordCalc();
        //StateCalc();
        //ForcesCalc();
        //AccelCalc();
        VelocCalc();

        BoundConitApply();
        //StatAnsApply();
        //SystemStateWrite();
    }
}

void Simulation::CoordCalc()
{
    for(size_t i = 0; i < params_.particles_amount; ++i)
    {
        for(int axis = 0; axis < 3; ++axis)
        {
            coord_[axis][i] += vel_[axis][i] * params_.dt + accel_[axis][i] * params_.dt * params_.dt / 2;
            vel_[axis][i] += accel_[axis][i] * params_.dt / 2;
        }
    }
}

void Simulation::VelocCalc()
{
    for(size_t i = 0; i < params_.particles_amount; ++i)
    {
        for(int axis = 0; axis < 3; ++axis)
        {
            vel_[axis][i] += accel_[axis][i] * params_.dt / 2;
        }
    }
}



void Simulation::BoundConitApply()
{
    for(int i = 0; i < 3; ++i)
    {
        if(params_.cell_bounds[i])
        {
            PeriodicConditApply(i);
        }
        else
        {
            IsolatedConditApply(i);
        }
    }
}


void Simulation::PeriodicConditApply(const int &axis)
{
    for(size_t i = 0; i < params_.particles_amount; ++i)
    {
        if(coord_[axis][i] > params_.cell_sizes[axis])
        {
            coord_[axis][i] -= static_cast<int>(coord_[axis][i] / params_.cell_sizes[axis]) * params_.cell_sizes[axis];
        }
        else if (coord_[axis][i] < 0)
        {
            coord_[axis][i] -= (static_cast<int>(coord_[axis][i] / params_.cell_sizes[axis]) - 1) * params_.cell_sizes[axis];
        }
    }

}


void Simulation::IsolatedConditApply(const int &axis)
{
    for(size_t i = 0; i < params_.particles_amount; ++i)
    {
        if(coord_[axis][i] > params_.cell_sizes[axis])
        {
            coord_[axis][i] -= (static_cast<int>(coord_[axis][i] / params_.cell_sizes[axis]) - 1) * params_.cell_sizes[axis];
        }
        else if (coord_[axis][i] < 0)
        {
            coord_[axis][i] -= static_cast<int>(coord_[axis][i] / params_.cell_sizes[axis]) * params_.cell_sizes[axis];
        }
        vel_[axis][i] = -vel_[axis][i];
    }
}

QVector<double> Simulation::RDF(int size)
{
    QVector<double> result;
    result.reserve(size);

    return result;
}
