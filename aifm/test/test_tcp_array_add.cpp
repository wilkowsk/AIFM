extern "C" {
#include <runtime/runtime.h>
}

#include "deref_scope.hpp"
#include "list.hpp"
#include "manager.hpp"

#include <iostream>

#include "dataframe_vector.hpp"
#include "device.hpp"
#include "helpers.hpp"

#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

using namespace far_memory;

constexpr uint64_t kCacheSize = (128ULL << 20);
constexpr uint64_t kFarMemSize = (4ULL << 30);
constexpr uint32_t kNumGCThreads = 12;
constexpr uint64_t kNumElements = 10000;

struct graph
{
    static int num_edges;
    static int num_nodes;

    // The node reached by vertex i's first outgoing edge is given by
    // outgoing_edges[outgoing_starts[i]].  To iterate over all
    // outgoing edges, please see the top-down bfs implementation.
    far_memory::DataFrameVector<int>* outgoing_starts;
    far_memory::DataFrameVector<int>* outgoing_edges; 

    far_memory::DataFrameVector<int>* incoming_starts;
    far_memory::DataFrameVector<int>* incoming_edges; 
};

namespace far_memory {

class FarMemTest {
public:
  void do_work(FarMemManager *manager) {
    {
      std::cout << "Running " << __FILE__ "..." << std::endl;
      DerefScope scope;
      List<int> list = FarMemManagerFactory::get()->allocate_list<int>(
          scope, /* enable_merge = */ true);
      List<int> list2 = FarMemManagerFactory::get()->allocate_list<int>(
          scope, /* enable_merge = */ true);
      List<int> outputList = FarMemManagerFactory::get()->allocate_list<int>(
          scope, /* enable_merge = */ true);

      for (uint64_t i = 0; i < kNumElements; i++) {
        list.push_back(scope, (1337*i));
      }
      for (uint64_t i = 0; i < kNumElements; i++) {
        list2.push_back(scope, (9001*i));
      }

      while (!list.empty() && !list2.empty()) {
        if (list.cfront(scope) < list2.cfront(scope)) {
          outputList.push_back(scope, list.cfront(scope));
          list.pop_front(scope);
        } else {
          outputList.push_back(scope, list2.cfront(scope));
          list2.pop_front(scope);
        }
      }
      while (!list.empty()) {
        outputList.push_back(scope, list.cfront(scope));
        list.pop_front(scope);
      }
      while (!list2.empty()) {
        outputList.push_back(scope, list2.cfront(scope));
        list2.pop_front(scope);
      }

      uint64_t element = 0;
      while (!outputList.empty()) {
        TEST_ASSERT(outputList.cfront(scope) >= element);
        element = outputList.cfront(scope);
        outputList.pop_front(scope);
      }
    }
    {
      DerefScope scope;
      graph besenj;
      far_memory::DataFrameVector<int> dataframe_vector = manager->allocate_dataframe_vector<int>();
      besenj.outgoing_starts;
      for (int i = 0; i < 5; i++) {
        int arraye[] = {0,4,6,7,8};
        besenj.outgoing_starts->push_back(scope, arraye[i]);
      }
      std::cout << std::endl;
      for (int i = 0; i < 5; i++) {
        int arraye[] = {0,2,3,5,7};
        besenj.incoming_starts->push_back(scope, arraye[i]);
      }
      std::cout << std::endl;
      for (int i = 0; i < 8; i++) {
        int arraye[] = {1,2,3,4,2,3,0,0};
        besenj.outgoing_edges->push_back(scope, arraye[i]);
      }
      std::cout << std::endl;
      for (int i = 0; i < 8; i++) {
        int arraye[] = {2,3,0,0,1,0,1,0};
        besenj.outgoing_starts->push_back(scope, arraye[i]);
      }
    }
    std::cout << "Passed" << std::endl;


    /*
    list.push_back(scope, 1);
    list.push_back(scope, 2);
    list.push_back(scope, 3);
    TEST_ASSERT(list.cfront(scope) == 1);
    TEST_ASSERT(list.cback(scope) == 3);
    TEST_ASSERT(list.size() == 3);
    TEST_ASSERT(list.empty() == false);

    int idx = 0;
    for (auto iter = list.begin(scope); iter != list.end(scope);
         iter.inc(scope)) {
      TEST_ASSERT(++idx == iter.deref(scope));
    }
    for (auto iter = list.rbegin(scope); iter != list.rend(scope);
         iter.inc(scope)) {
      TEST_ASSERT(idx-- == iter.deref(scope));
    }

    list.pop_front(scope);
    list.pop_back(scope);
    TEST_ASSERT(list.cfront(scope) == 2);
    TEST_ASSERT(list.cback(scope) == 2);
    TEST_ASSERT(list.size() == 1);
    TEST_ASSERT(list.empty() == false);

    list.pop_front(scope);
    TEST_ASSERT(list.size() == 0);
    TEST_ASSERT(list.empty() == true);

    list.push_back(scope, 1);
    list.push_back(scope, 2);
    list.push_back(scope, 3);

    auto iter = list.begin(scope);
    idx = 0;
    while (iter != list.end(scope)) {
      TEST_ASSERT(++idx == iter.deref(scope));
      iter = list.erase(scope, iter);
    }

    for (uint32_t i = 1; i <= 2 * List<int>::kMaxNumNodesPerChunk; i++) {
      list.push_back(scope, i);
    }

    iter = list.begin(scope);
    idx = 0;
    while (iter != list.end(scope)) {
      TEST_ASSERT(++idx == iter.deref(scope));
      iter = list.erase(scope, iter);
    }
    TEST_ASSERT(idx == 2 * List<int>::kMaxNumNodesPerChunk);

    for (uint32_t i = 1; i <= 2 * List<int>::kMaxNumNodesPerChunk; i++) {
      list.push_front(scope, i);
    }

    auto rev_iter = list.rbegin(scope);
    idx = 0;
    while (rev_iter != list.rend(scope)) {
      TEST_ASSERT(++idx == rev_iter.deref(scope));
      rev_iter = list.erase(scope, rev_iter);
    }
    TEST_ASSERT(idx == 2 * List<int>::kMaxNumNodesPerChunk);

    for (uint32_t i = 1; i <= List<int>::kMaxNumNodesPerChunk - 1; i++) {
      list.push_back(scope, i);
    }
    iter = list.begin(scope);
    iter.inc(scope);
    iter.inc(scope);
    list.insert(scope, &iter, 0);

    idx = 0;
    iter = list.begin(scope);
    while (iter != list.end(scope)) {
      auto list_data = iter.deref(scope);
      TEST_ASSERT(++idx == list_data);
      iter = list.erase(scope, iter);
      if (unlikely(list_data == 2)) {
        TEST_ASSERT(0 == iter.deref(scope));
        iter = list.erase(scope, iter);
      }
    }
    TEST_ASSERT(idx == List<int>::kMaxNumNodesPerChunk - 1);
    TEST_ASSERT(list.empty());

    for (uint32_t i = 1; i <= List<int>::kMaxNumNodesPerChunk - 1; i++) {
      list.push_back(scope, i);
    }
    iter = list.end(scope);
    iter.dec(scope);
    iter.dec(scope);
    list.insert(scope, &iter, 0);

    idx = 0;
    iter = list.begin(scope);
    while (iter != list.end(scope)) {
      auto list_data = iter.deref(scope);
      TEST_ASSERT(++idx == list_data);
      iter = list.erase(scope, iter);
      if (unlikely(list_data == List<int>::kMaxNumNodesPerChunk - 3)) {
        TEST_ASSERT(0 == iter.deref(scope));
        iter = list.erase(scope, iter);
      }
    }
    TEST_ASSERT(idx == List<int>::kMaxNumNodesPerChunk - 1);

    for (uint32_t i = 1; i <= List<int>::kMaxNumNodesPerChunk + 1; i++) {
      list.push_back(scope, i);
    }

    iter = list.end(scope);
    for (uint32_t i = 0; i <= List<int>::kMaxNumNodesPerChunk / 2; i++) {
      iter.dec(scope);
    }
    TEST_ASSERT(33 == iter.deref(scope));
    for (uint32_t i = 0; i <= List<int>::kMaxNumNodesPerChunk / 4; i++) {
      iter = list.erase(scope, iter);
    }
    iter = list.begin(scope);
    for (uint32_t i = 0; i < List<int>::kMaxNumNodesPerChunk / 4; i++) {
      iter = list.erase(scope, iter);
    }
    TEST_ASSERT(list.local_list_.size() == 3);
    for (uint32_t i = 0; i < List<int>::kMaxNumNodesPerChunk / 4; i++) {
      TEST_ASSERT(i + 17 == static_cast<uint32_t>(iter.deref(scope)));
      iter = list.erase(scope, iter);
    }
    for (uint32_t i = 0; i < List<int>::kMaxNumNodesPerChunk / 4; i++) {
      TEST_ASSERT(i + 50 == static_cast<uint32_t>(iter.deref(scope)));
      iter = list.erase(scope, iter);
    }
    TEST_ASSERT(list.empty());

    std::cout << "Passed" << std::endl;
    */
  }
};
} // namespace far_memory

void _main(void *arg) {
  std::unique_ptr<FarMemManager> manager =
      std::unique_ptr<FarMemManager>(FarMemManagerFactory::build(
          kCacheSize, kNumGCThreads, new FakeDevice(kFarMemSize)));
  FarMemTest test;
  test.do_work(manager.get());
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
