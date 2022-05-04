extern "C" {
#include <runtime/runtime.h>
}

#include "deref_scope.hpp"
#include "list.hpp"
#include "manager.hpp"

#include <iostream>
#include <memory>

using namespace far_memory;

struct Data {
  uint32_t data;
  uint32_t dummy[1023];

  Data(uint32_t _data) : data(_data) {}
};

constexpr uint64_t kCacheSize = (256ULL << 20);
constexpr uint64_t kFarMemSize = (8ULL << 30);
constexpr uint32_t kNumGCThreads = 12;
constexpr uint32_t kNumDataEntries = 8 * kCacheSize / sizeof(Data);
constexpr uint32_t kScopeResetInterval = 256;

namespace far_memory {

class FarMemTest {
public:
  void do_work(FarMemManager *manager) {
    std::cout << "Running " << __FILE__ "..." << std::endl;
    DerefScope scope;
    Queue<Data> queue =
        FarMemManagerFactory::get()->allocate_queue<Data>(scope);
    uint32_t length = 0;
    Queue<Data> queue2 =
        FarMemManagerFactory::get()->allocate_queue<Data>(scope);
    uint32_t length2 = 0;
    Queue<Data> outputQueue =
        FarMemManagerFactory::get()->allocate_queue<Data>(scope);

    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
      	scope.renew();
      }
      queue.push(scope, Data(i*1337));
      length++;
    }
    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
      	scope.renew();
      }
      queue2.push(scope, Data(i*9001));
      length2++;
    }
	  
    uint32_t outputLength = 0;
    while ((length != 0) && (length2 != 0)) {
      if (unlikely(outputLength % kScopeResetInterval == 0)) {
      	scope.renew();
      }
      if (queue.cfront(scope).data < queue2.cfront(scope).data) {
        outputQueue.push(scope, Data(queue.cfront(scope).data));
        length--;
      } else {
        outputQueue.push(scope, Data(queue2.cfront(scope).data));
        length2--;
      }
      outputLength++;
    }
    while (length != 0) {
      outputQueue.push(scope, Data(queue.cfront(scope).data));
      length--;
    }
    while (length2 != 0) {
      outputQueue.push(scope, Data(queue2.cfront(scope).data));
      length2--;
    }

    uint32_t element = 0;
    for (uint32_t i = 0; i < outputLength; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
      	scope.renew();
      }
      TEST_ASSERT(outputQueue.cfront(scope).data >= element);
      element = outputQueue.cfront(scope).data;
      outputQueue.pop(scope);
    }

    std::cout << "Passed" << std::endl;
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
