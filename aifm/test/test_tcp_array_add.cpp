extern "C" {
#include <runtime/runtime.h>
}

#include "array.hpp"
#include "device.hpp"
#include "helpers.hpp"
#include "manager.hpp"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <string>

using namespace far_memory;
using namespace std;

constexpr static uint64_t kCacheSize = (128ULL << 20);
constexpr static uint64_t kFarMemSize = (4ULL << 30);
constexpr static uint32_t kNumGCThreads = 12;
constexpr static uint32_t kNumEntries =
    (8ULL << 20); // So the array size is larger than the local cache size.
constexpr static uint32_t kNumConnections = 300;

uint64_t raw_array_A[kNumEntries];

template <uint64_t N, typename T>
void copy_array(Array<T, N> *array, T *raw_array) {
  for (uint64_t i = 0; i < N; i++) {
    DerefScope scope;
    (*array).at_mut(scope, i) = raw_array[i];
  }
}

template <typename T, uint64_t N>
void add_array(Array<T, N> *array_C, Array<T, N> *array_A,
               Array<T, N> *array_B) {
  for (uint64_t i = 0; i < N; i++) {
    DerefScope scope;
    (*array_C).at_mut(scope, i) =
        (*array_A).at(scope, i) + (*array_B).at(scope, i);
  }
}

template <typename T, uint64_t N>
int search_array(Array<T, N> *array_A, uint64_t targetElement) {
  uint64_t low = 0;
  uint64_t high = N - 1;
  uint64_t middle;
  uint64_t middleElement;
  bool foundElement = false;
  while ((high >= low) && !foundElement) {
    DerefScope scope;
    middle = (high + low) / 2;
    middleElement = (*array_A).at(scope, middle);
    if (middleElement > targetElement) { // middleElement is too high
      high = middle - 1;
    } elseif (middleElement < targetElement) { // middleElement is too low
      low = middle + 1;
    } else {
      foundElement = true;
    }
  }
  if (foundElement) {
    return middle;
  }
  return -1;
}

void gen_random_array(uint64_t num_entries, uint64_t *raw_array) {
  std::random_device rd;
  std::mt19937_64 eng(rd());
  std::uniform_int_distribution<uint64_t> distr(2,257);
  uint64_t value = 0;
  for (uint64_t i = 0; i < num_entries; i++) {
    value += distr(eng);
    raw_array[i] = value;
  }
}

void do_work(FarMemManager *manager) {
  auto array_A = manager->allocate_array<uint64_t, kNumEntries>();

  gen_random_array(kNumEntries, raw_array_A);
  copy_array(&array_A, raw_array_A);
  add_array(&array_C, &array_A, &array_B);
  uint64_t result;

  for (uint64_t i = 0; i < 10000; i += 100) {
    DerefScope scope;
    result = search_array(&array_A, raw_array_A[i]);
    if (result != i) {
      goto fail;
    }
    result = search_array(&array_A, raw_array_A[i] + 1);
    if (result != -1) {
      goto fail;
    }
  }

  cout << "Passed" << endl;
  return;

fail:
  cout << "Failed" << endl;
}

int argc;
void _main(void *arg) {
  char **argv = static_cast<char **>(arg);
  std::string ip_addr_port(argv[1]);
  auto raddr = helpers::str_to_netaddr(ip_addr_port);
  std::unique_ptr<FarMemManager> manager =
      std::unique_ptr<FarMemManager>(FarMemManagerFactory::build(
          kCacheSize, kNumGCThreads,
          new TCPDevice(raddr, kNumConnections, kFarMemSize)));
  do_work(manager.get());
}

int main(int _argc, char *argv[]) {
  int ret;

  if (_argc < 3) {
    std::cerr << "usage: [cfg_file] [ip_addr:port]" << std::endl;
    return -EINVAL;
  }

  char conf_path[strlen(argv[1]) + 1];
  strcpy(conf_path, argv[1]);
  for (int i = 2; i < _argc; i++) {
    argv[i - 1] = argv[i];
  }
  argc = _argc - 1;

  ret = runtime_init(conf_path, _main, argv);
  if (ret) {
    std::cerr << "failed to start runtime" << std::endl;
    return ret;
  }

  return 0;
}
