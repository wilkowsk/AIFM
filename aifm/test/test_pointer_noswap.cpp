extern "C" {
#include <runtime/runtime.h>
}

#include "deref_scope.hpp"
#include "device.hpp"
#include "manager.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>

using namespace far_memory;
using namespace std;

struct Data512 {
  uint64_t data[513];
};

void do_work(FarMemManager *manager) {
  cout << "Running " << __FILE__ "..." << endl;

  auto far_mem_ptr_0 = manager->allocate_unique_ptr<Data512>();

  DerefScope scope;
  {
    auto raw_ptr_0 = far_mem_ptr_0.deref_mut(scope);
    for (uint64_t i = 0; i < 513; i++) {
      raw_ptr_0->data[i] = i;
    }
  }

  {
    const auto raw_ptr_0 = far_mem_ptr_0.deref(scope);
    for (uint64_t i = 0; i < 513; i++) {
      uint64_t low = 0;
      uint64_t high = 512;
      uint64_t middle;
      bool foundElement = false;
      while ((high >= low) && !foundElement) {
        middle = (low + high) / 2;
        if (raw_ptr_0->data[middle] > i) { // update high
          high = middle - 1;
        } else if (raw_ptr_0->data[middle] < i) { // update low
          low = middle + 1;
        } else {
          foundElement = true;
        }
      }
      if (raw_ptr_0->data[middle] != i) {
        goto fail;
      }
      if (!foundElement) {
        goto fail;
      }
      //if (raw_ptr_0->data[i] != static_cast<char>(i)) {
      //  goto fail;
      //}
    }
  }

  cout << "Passed" << endl;
  return;

fail:
  cout << "Failed" << endl;
  return;
}

void _main(void *arg) {
  uint64_t cache_size = 1ULL << 32;
  uint64_t far_mem_size = 1ULL << 33;
  uint8_t num_gc_threads = 12;

  auto manager = std::unique_ptr<FarMemManager>(FarMemManagerFactory::build(
      cache_size, num_gc_threads, new FakeDevice(far_mem_size)));
  do_work(manager.get());
}

int main(int argc, char *argv[]) {
  int ret;

  if (argc < 2) {
    std::cerr << "usage: [cfg_file]" << std::endl;
    return -EINVAL;
  }

  ret = runtime_init(argv[1], _main, NULL);
  if (ret) {
    std::cerr << "failed to start runtime" << std::endl;
    return ret;
  }

  return 0;
}
