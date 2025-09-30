//
// Created by Gregory Johnson on 7/15/25.
//

#ifndef CLIPROGRESSBAR_H
#define CLIPROGRESSBAR_H
#include <Rcpp.h>
#include <RcppThread.h>
#include "progress_bar.hpp"
class CliProgressBar : public ProgressBar {
public:
    CliProgressBar() { reset(); }

    ~CliProgressBar() = default;

public:
    void display() override { RcppThread::Rcout << "\033[37mComputing: "; }

    void update(const float progress) override {
        if (_firstTime) {
            _display_ticks(progress);
            std::time(&start);
            _firstTime = false;
        }

        std::time(&end);
        _update_ticks_display(progress);
        if (_ticks_displayed >= _max_ticks) _finalize_display();
    }

    void end_display() override {
        update(1);
        reset();
    }

    void reset() {
        // Rcpp::Environment cli = Rcpp::Environment::namespace_env("cli");
        // const Rcpp::Function console_width = cli["console_width"];

        // From fiddling around with it, it seems that dividing the console
        // width by 2 produces the best display.
        _max_ticks = 50;//std::floor(Rcpp::as<int>(console_width()) / 2);
        _ticks_displayed = 0;
        _finalized = false;
    }

protected:
    void _update_ticks_display(const float progress) {
        const int nb_ticks = _compute_nb_ticks(progress);
        const int delta = nb_ticks - _ticks_displayed;
        if (delta > 0) {
            RcppThread::Rcout << "\r" << std::flush;
            _ticks_displayed = nb_ticks;
            _display_ticks(progress);
        }
    }

    void _finalize_display() {
        if (_finalized) return;
        RcppThread::Rcout << std::endl << std::flush;
        _finalized = true;
    }

    int _compute_nb_ticks(const float progress) const { return static_cast<int>(progress * _max_ticks); }

    void _display_ticks(const double progress) const {
        RcppThread::Rcout << "\033[37mComputing ";
        // calculate passed time and remaining time (in seconds)
        const double pas_time = std::difftime(end, start);
        const double rem_time = (pas_time / progress) * (1 - progress);

        // convert seconds to time string

        std::string time_string = _time_to_string(rem_time);
        if (_firstTime)
            time_string = "-";
        for (int i = 0; i < _ticks_displayed; ++i) {
            RcppThread::Rcout << "\033[32m\u25A0" << std::flush;
        }
        for (int i = 0; i < (_max_ticks - _ticks_displayed); i++) {
            RcppThread::Rcout << " " << std::flush;
        }
        RcppThread::Rcout << "\033[37m | " << static_cast<int>(progress * 100) << "%" << "  ETA: "
        << time_string << "...";
    }
    static std::string _time_to_string(const double seconds) {

        int time = static_cast<int>(seconds);

        int hour = 0;
        int min = 0;
        int sec = 0;

        hour = time / 3600;
        time = time % 3600;
        min = time / 60;
        time = time % 60;
        sec = time;

        //std::stringstream time_strs;
        std::string timeString;
        if (hour != 0) timeString += std::to_string(hour) + "h ";
        if (min != 0) timeString += std::to_string(min) + "min ";
        if (sec != 0) timeString += std::to_string(sec) + "s ";

        // if (hour != 0) time_strs << hour << "h ";
        // if (min != 0) time_strs << min << "min ";
        // if (sec != 0) time_strs << sec << "s ";
        // std::string time_str = time_strs.str();

        return timeString;
    }


private:
    int _max_ticks;        // the total number of ticks to print
    int _ticks_displayed;  // the nb of ticks already displayed
    bool _finalized;
    bool _firstTime = true;
    time_t start, end;
};
#endif //CLIPROGRESSBAR_H
