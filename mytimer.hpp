#ifndef MYTIMER_HPP
#define MYTIMER_HPP

#include <chrono>

class MyTimer
{
public:
    MyTimer();
    double fromStart();  // Returns seconds since the object creation (like 0.00023s)
    double fromLast();  // Returns seconds since last call of "fromLast" or "fromStart" methodes
private:
    const std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> last;
};

#endif // MYTIMER_HPP
