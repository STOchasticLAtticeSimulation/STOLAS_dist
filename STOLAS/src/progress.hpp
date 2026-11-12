#ifndef INCLUDED_progress_hpp_
#define INCLUDED_progress_hpp_

#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

// Progress reporter for long-running loops (thread-safe, usable inside OpenMP loops).
//   interactive terminal : one line is overwritten in place at every 1% step
//   redirected / job log : a new line is printed every `logStep` percent
//
//   Progress prog("dN map", total);
//   prog.tick();       // from any thread, once per finished item
//   prog.set(done);    // or give the absolute count (serial loops)
//   prog.finish();     // prints the total elapsed time and ends the line
class Progress {
public:
  Progress(const std::string& label, long total, int logStep = 10)
    : label_(label), total_(total), logStep_(logStep > 0 ? logStep : 10),
      tty_(isatty(STDOUT_FILENO) != 0), t0_(std::chrono::steady_clock::now()) {}

  void tick(long n = 1) { report(count_.fetch_add(n) + n); }

  void set(long done) {
    count_.store(done);
    report(done);
  }

  void finish() {
    std::lock_guard<std::mutex> lock(mtx_);
    if (tty_) std::cout << "\r\033[K";
    std::cout << std::left << std::setw(kLabelWidth) << label_
              << ": done in " << format(elapsed()) << std::endl;
  }

private:
  static constexpr int kLabelWidth = 28;

  std::string label_;
  long total_;
  int logStep_;
  bool tty_;
  std::chrono::steady_clock::time_point t0_;
  std::atomic<long> count_{0};
  std::atomic<int> lastKey_{-1};
  std::mutex mtx_;

  double elapsed() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
  }

  static std::string format(double sec) {
    long s = (long)sec;
    std::ostringstream os;
    os << std::setfill('0') << std::setw(2) << s / 3600 << ':'
       << std::setw(2) << (s % 3600) / 60 << ':' << std::setw(2) << s % 60;
    return os.str();
  }

  void report(long done) {
    int pct = total_ > 0 ? (int)std::min<long>(100, 100 * done / total_) : 100;
    int key = tty_ ? pct : pct / logStep_ * logStep_;
    int prev = lastKey_.load();
    while (key > prev) {
      if (lastKey_.compare_exchange_weak(prev, key)) {
        if (tty_ || pct < 100) print(pct, done); // in log mode finish() reports completion
        return;
      }
    }
  }

  void print(int pct, long done) {
    double el = elapsed();
    std::ostringstream os;
    os << std::left << std::setw(kLabelWidth) << label_ << ": "
       << std::right << std::setw(3) << pct << "%"
       << " | elapsed " << format(el);
    if (done > 0 && done < total_) os << " | ETA " << format(el * (total_ - done) / done);

    std::lock_guard<std::mutex> lock(mtx_);
    if (tty_) std::cout << '\r' << os.str() << "\033[K" << std::flush;
    else std::cout << os.str() << std::endl;
  }
};

#endif
