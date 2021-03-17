#include <chrono>
#include <string>
#include <vector>
#include <iostream>

class Timer {
 public:
  Timer(const std::string& delimiter = "") : delimiter_(delimiter) {}

  std::chrono::time_point<std::chrono::high_resolution_clock> Mark(
      const std::string& marker) {
    time_points_.emplace_back(marker, std::chrono::high_resolution_clock::now());
  }

  std::vector<std::pair<std::string, double>> GetDurations() const {
    std::vector<std::pair<std::string, double>> durations;
    for (int i = 1; i < static_cast<int>(time_points_.size()); ++i) {
      auto [s1, t1] = time_points_[i - 1];
      auto [s2, t2] = time_points_[i];
      durations.emplace_back(
          s2, std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                  .count());
    }
    return durations;
  }

  void Report() const {
    std::vector<std::pair<std::string, double>> durations = GetDurations();

    double total_duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            time_points_.back().second - time_points_.front().second)
            .count();
    std::cout << delimiter_ << "Timings : {\n";
    for (auto [s, d] : durations) {
      std::cout << delimiter_ << "    " << s << ": " << d << " ms ( " <<
          (100 * d / total_duration) << " % )\n";
    }
    std::cout << delimiter_ << "}\n";
  }

 private:
  std::vector<std::pair<std::string,
      const std::chrono::time_point<std::chrono::high_resolution_clock>>>
      time_points_;
  std::string delimiter_;
};