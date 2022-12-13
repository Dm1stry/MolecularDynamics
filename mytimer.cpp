#include "mytimer.hpp"

MyTimer::MyTimer() : start(std::chrono::system_clock::now()), last(start) {}

double MyTimer::fromStart()
{
    last = std::chrono::system_clock::now();
    return std::chrono::duration<double>(last - start).count();
}

double MyTimer::fromLast()
{
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    double time = std::chrono::duration<double>(now - last).count();
    last = now;
    return time;
}
