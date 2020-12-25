#include "index.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#ifdef __linux__
#include <sched.h>
#define CHECK_S(ret, s)                                                                                                \
  if (ret == -1) {                                                                                                     \
    perror((s));                                                                                                       \
    exit(EXIT_FAILURE);                                                                                                \
  }
cpu_set_t orig_set;
pid_t pid = 0;
void set_affinity() {
  CHECK_S(sched_getaffinity(pid, sizeof(cpu_set_t), &orig_set), "Could not get CPU affinity");
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(0, &set);
  CHECK_S(sched_setaffinity(pid, sizeof(cpu_set_t), &set), "Could not set CPU affinity");
}
void reset_affinity() { CHECK_S(sched_setaffinity(pid, sizeof(cpu_set_t), &orig_set), "Could not set CPU affinity"); }
#else
void set_affinity() {}
void reset_affinity() {}
#endif

template <typename T> std::vector<T> read_data_binary(const std::string &path) {
  try {
    std::fstream in(path, std::ios::in | std::ios::binary);
    in.exceptions(std::ios::failbit | std::ios::badbit);

    size_t size = 0;
    in.read((char *)&size, sizeof(T));
    size/=20;
    std::vector<T> data(size);
    in.read((char *)data.data(), size * sizeof(T));
    return data;
  } catch (std::ios_base::failure &e) {
    std::cerr << "Error reading the input data: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
}

size_t query_time(const MyIndex &index, const std::vector<uint64_t> &queries) {
  auto start = std::chrono::high_resolution_clock::now();

  uint64_t cnt = 0;
  for (auto &q : queries)
    cnt += index.nextGEQ(q);
  asm volatile("" : : "r,m"(cnt) : "memory");

  auto stop = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / queries.size();
}

void test(const MyIndex &index, const std::vector<uint64_t> &data, const std::string &name) {
  // Lookup test
  for (size_t i = 0; i < data.size(); ++i) {
    if (i > 0 && data[i] == data[i - 1])
      continue;
    auto q = data[i];
    auto returned = index.nextGEQ(q);
    if (q != returned) {
      std::cerr << "Test failed on " << name << ". For the query key " << q << " the index returned " << returned
                << " but the correct answer is " << q << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // nextGEQ test
  std::uniform_int_distribution<uint64_t> distribution(data.front(), data.back() - 2);
  std::mt19937 g(42);
  for (size_t i = 0; i < 15000000; ++i) {
    auto q = distribution(g);
    auto correct = *std::lower_bound(data.begin(), data.end(), q);
    auto returned = index.nextGEQ(q);
    if (correct != returned) {
      std::cerr << "Test failed on " << name << ". For the query key " << q << " the index returned " << returned
                << " but the correct answer is " << correct << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

int main(int argc, char **argv) {
  std::cout << "dataset,nanoseconds,bytes" << std::endl;

  for (int argi = 1; argi < argc; ++argi) {
    auto path = std::string(argv[argi]);
    auto name = path.substr(path.find_last_of("/\\") + 1);
    auto data = read_data_binary<uint64_t>(path);

    MyIndex index(data);

    set_affinity();
    test(index, data, name);

    // Benchmark
    std::vector<uint64_t> queries(10000000);
    std::mt19937 g(42);
    std::uniform_int_distribution<uint64_t> distribution(data.front(), data.back() - 2);
    std::generate(queries.begin(), queries.end(), [&]() { return distribution(g); });
    std::cout << name << "," << query_time(index, queries) << "," << index.size_in_bytes() << std::endl;
    reset_affinity();
  }

  return 0;
}