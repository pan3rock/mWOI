#ifndef TIMER_H_
#define TIMER_H_
#include <chrono>
#include <ctime>

class Timer {
public:
  void start() {
    m_StartTime = std::chrono::system_clock::now();
    m_bRunning = true;
  }

  void stop() {
    m_EndTime = std::chrono::system_clock::now();
    m_bRunning = false;
  }

  double elapsedNanoseconds() {
    std::chrono::time_point<std::chrono::system_clock> endTime;

    if (m_bRunning) {
      endTime = std::chrono::system_clock::now();
    } else {
      endTime = m_EndTime;
    }

    return std::chrono::duration_cast<std::chrono::nanoseconds>(endTime -
                                                                m_StartTime)
        .count();
  }

  double elapsedMilliseconds() { return elapsedNanoseconds() / 1000.0; }
  double elapsedSeconds() { return elapsedNanoseconds() / 1.0e6; }

private:
  std::chrono::time_point<std::chrono::system_clock> m_StartTime;
  std::chrono::time_point<std::chrono::system_clock> m_EndTime;
  bool m_bRunning = false;
};
#endif