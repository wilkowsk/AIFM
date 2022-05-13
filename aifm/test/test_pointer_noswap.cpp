extern "C" {
#include <runtime/runtime.h>
}

#include "array.hpp"
#include "helpers.hpp"

#include "deref_scope.hpp"
#include "device.hpp"
#include "manager.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>

#include <cstdint>
#include <cstring>
#include <random>
#include <string>

#include <typeinfo>

using namespace far_memory;
using namespace std;

struct Data512 {
  uint64_t data[513];
};

void do_work(FarMemManager *manager) {
  cout << "Running " << __FILE__ "..." << endl;

  auto far_mem_ptr_0 = manager->allocate_shared_ptr<Data512>();

  auto partArray = manager->allocate_array<SharedPtr<Data512>, 10000>();

  {
    DerefScope scope;
    auto raw_ptr_0 = far_mem_ptr_0.deref_mut(scope);
    for (uint64_t i = 0; i < 513; i++) {
      raw_ptr_0->data[i] = i;
    }
  }

  {
    DerefScope scope;
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


  for (uint64_t partId = 0; partId < 10000; partId++)
    {
      DerefScope scope;
      auto &pointer_loc = partArray.at_mut(scope,partId);
      auto pointer_loc_ptr = &pointer_loc;
      //auto raw_arrary_elem = pointer_loc_ptr->deref_mut(scope);
      auto part = manager->allocate_shared_ptr<Data512>();
      //cout<< "Type checking inside part: " <<typeid(part).name()<<"\n";
      Data512* raw_part = part.deref_mut(scope);
      raw_part->data[0] = partId;
      //cout<< "PartId <<<<<<<<<<<<<<<<<<< " << raw_part->partId <<"\n";
      *pointer_loc_ptr = part;
    }

  for (uint64_t partId = 0; partId < 10000; partId++) {
    DerefScope scope;
    uint64_t low = 0;
    uint64_t high = 9999;
    uint64_t middle;
    bool foundElement = false;
    while (high >= low && !foundElement) {
      middle = (high + low) / 2;
      auto &pointer_loc = partArray.at_mut(scope,middle);
      auto pointer_loc_ptr1 = &pointer_loc;
      const Data512* raw_part = pointer_loc.deref(scope);
      if (raw_part->data[0] < partId) {
        low = middle + 1;
      } else if (raw_part->data[0] > partId) {
        high = middle - 1;
      } else {
        foundElement = true;
      }
    }
    if (!foundElement) {
      goto fail;
    }
    if (middle != partId) {
      goto fail;
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
