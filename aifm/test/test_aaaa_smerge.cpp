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
    Stack<Data> stack =
        FarMemManagerFactory::get()->allocate_stack<Data>(scope);
    uint32_t length = 0;
    Stack<Data> stack2 =
        FarMemManagerFactory::get()->allocate_stack<Data>(scope);
    uint32_t length2 = 0;
    Stack<Data> outputStack =
        FarMemManagerFactory::get()->allocate_stack<Data>(scope);
    uint32_t outputLength = 0;
    
    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      stack.push(scope, Data(1337*i));
      length++;
    }
    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      stack2.push(scope, Data(9001*i));
      length2++;
    }
    while (length != 0 && length2 != 0) {
      if (unlikely(outputLength % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      if (stack.ctop(scope).data > stack2.ctop(scope).data) {
        outputStack.push(scope, Data(stack.ctop(scope).data));
        length--;
      } else {
        outputStack.push(scope, Data(stack2.ctop(scope).data));
        length2--;
      }
      outputLength++;
    }
    
    while (length != 0) {
      if (unlikely(outputLength % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      outputStack.push(scope, Data(stack.ctop(scope).data));
      length--;
      outputLength++;
    }
    while (length2 != 0) {
      if (unlikely(outputLength % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      outputStack.push(scope, Data(stack2.ctop(scope).data));
      length2--;
      outputLength++;
    }
    
    uint32_t element = outputStack.ctop(scope).data;
    while (outputLength != 0) {
      if (unlikely(outputLength % kScopeResetInterval == 0)) {
	      scope.renew();
      }
      TEST_ASSERT(outputStack.ctop(scope).data >= element);
      element = outputStack.ctop(scope).data;
      outputStack.pop(scope);
      outputLength--;
    }
    
    std::cout << "Passed" << std::endl;
	  /*
    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
	scope.renew();
      }
      stack.push(scope, Data(i));
    }

    for (uint32_t i = 0; i < kNumDataEntries; i++) {
      if (unlikely(i % kScopeResetInterval == 0)) {
	scope.renew();
      }
      TEST_ASSERT(stack.ctop(scope).data == kNumDataEntries - i - 1);
      stack.pop(scope);
    }

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
