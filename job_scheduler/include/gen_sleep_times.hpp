#pragma once
#include <vector>
#include <random>
#include <algorithm>

void generate_sleep_times(std::vector<int>& sleep_times, int num_tasks) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10);

    sleep_times.resize(num_tasks);
    std::generate(sleep_times.begin(), sleep_times.end(), [&]() { return dis(gen); });
}