// Timer.h
#ifndef TIMER_H
#define TIMER_H

#include "Core_Export.h"
#include <chrono>
#include <string>

namespace Poisson{

/// @brief 计时器
class POISSONCORE_API Timer {
public:
    Timer() = default;
    ~Timer() = default;

    void start();
    void pause();
    void reset();
    void print_time();
    std::string format() const;

    template <typename Unit = std::chrono::milliseconds>
    typename Unit::rep elapsed() const {
        if (m_running) {
            auto now = std::chrono::steady_clock::now();
            return std::chrono::duration_cast<Unit>(m_accumulated + (now - m_start)).count();
        }
        return std::chrono::duration_cast<Unit>(m_accumulated).count();
    }

private:
    std::chrono::steady_clock::time_point m_start;
    std::chrono::steady_clock::duration m_accumulated{0};
    bool m_running = false;
};

}// namespace Poisson

#endif // TIMER_H