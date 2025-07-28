// Timer.cpp
#include "Timer.h"
#include "Core_Export.h"
#include <iostream>

namespace Poisson{

/// @brief 计时器，记录事件开始的时间
void Timer::start() {
    if (m_running) {
        std::cerr << "[Timer] Restarting while running, previous measurement lost\n";
    }
    m_start = std::chrono::steady_clock::now();
    m_running = true;
}

/// @brief 暂停记录，计算事件经过的时间
void Timer::pause() {
    if (!m_running) return;
    
    auto now = std::chrono::steady_clock::now();
    m_accumulated += now - m_start;
    m_running = false;
}

/// @brief 重置时间
void Timer::reset() {
    m_accumulated = std::chrono::steady_clock::duration::zero();
    m_running = false;
}

/// @brief 输出事件运行时间
void Timer::print_time(){
    std::cout << "\n程序运行总时间为: " << Timer::format() << std::endl; 
}

/// @brief 获取格式化时间
/// @return 返回格式化时间字符串
std::string Timer::format() const {
    using namespace std::chrono;
    
    nanoseconds total = m_accumulated;
    if (m_running) {
        total += duration_cast<nanoseconds>(steady_clock::now() - m_start);
    }

    if (total < microseconds(1)) {
        return std::to_string(total.count()) + " ns";
    } else if (total < milliseconds(1)) {
        return std::to_string(duration_cast<microseconds>(total).count()) + " μs";
    } else if (total < seconds(1)) {
        return std::to_string(duration_cast<milliseconds>(total).count()) + " ms";
    } else {
        return std::to_string(duration_cast<seconds>(total).count()) + " sec";
    }
}

} //namespace Poisson