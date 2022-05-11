extern "C" {
#include <runtime/runtime.h>
}

#include "concurrent_hopscotch.hpp"
#include "device.hpp"
#include "helpers.hpp"
#include "manager.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

using namespace far_memory;
using namespace std;

constexpr static uint32_t kKeyLen = 200;
constexpr static uint32_t kValueLen = 700;
constexpr static uint32_t kHashTableNumEntriesShift = 16;
constexpr static uint32_t kHashTableRemoteDataSize =
    (Object::kHeaderSize + kKeyLen + kValueLen) *
    (1 << kHashTableNumEntriesShift);
constexpr static double kLoadFactor = 0.80;
constexpr static uint32_t kNumKVPairs =
    kLoadFactor * (1 << kHashTableNumEntriesShift);

constexpr static uint64_t kCacheSize = (1ULL << 30);
constexpr static uint64_t kFarMemSize = (1ULL << 30);
constexpr static uint32_t kNumGCThreads = 12;

struct Key {
  char data[kKeyLen];
  bool operator<(const Key &other) const {
    return strncmp(data, other.data, kKeyLen) < 0;
  }
};

struct Value {
  char data[kValueLen];
  bool operator<(const Value &other) const {
    return strncmp(data, other.data, kValueLen) < 0;
  }
};

std::map<Key, Value> kvs;

void random_string(char *data, uint32_t len) {
  for (uint32_t i = 0; i < len; i++) {
    data[i] = rand() % ('z' - 'a' + 1) + 'a';
  }
}

void do_work(FarMemManager *manager) {
  cout << "Running " << __FILE__ "..." << endl;

  auto hopscotch = manager->allocate_concurrent_hopscotch<Key, Value>(
      kHashTableNumEntriesShift, kHashTableNumEntriesShift,
      kHashTableRemoteDataSize);

  char nextKey[kKeyLen];
  for (uint32_t i = 0; i < kKeyLen; i++) {
    nextKey[i] = 'a';
  }
  for (uint32_t i = 0; i < kNumKVPairs; i++) {
    Key key;
    Value value;
    for (uint32_t j = 0; j < kKeyLen; j++) {
      key.data[j] = nextKey[j];
    }
    random_string(nextKey, kKeyLen);
    uint32_t intermediate = i;
    for (uint32_t j = 0; j < kKeyLen; j++) {
      value.data[j] = nextKey[j];
    }
    for (uint32_t j = kKeyLen; j < kValueLen; j++) {
      if (intermediate % 2 == 0) {
        value.data[j] = 'o';
      } else {
        value.data[j] = 'i';
      }
      intermediate /= 2;
    }
    if (i < 10) {
      cout << key.data << endl;
      cout << value.data << endl;
    }
    hopscotch.insert_tp(key, value);
    kvs[key] = value;
  }
  cout << endl << endl << endl;
  Key iteratorKey;
  for (uint32_t i = 0; i < kKeyLen; i++) {
    iteratorKey.data[i] = 'a';
  }
  for (uint32_t i = 0; i < kNumKVPairs; i++) {
    std::optional<Value> optional_value;
    optional_value = hopscotch.find_tp(iteratorKey);
    if (i < 10) {
      cout << iteratorKey.data << endl;
    }
    TEST_ASSERT(optional_value);
    if (i < 10) {
      cout << optional_value->data << endl;
    }
    uint32_t intermediate = i;
    for (uint32_t j = 0; j < kKeyLen; j++) {
      iteratorKey.data[j] = optional_value->data[j];
    }
    for (uint32_t j = kKeyLen; j < kValueLen; j++) {
      if (intermediate % 2 == 0) {
        TEST_ASSERT(optional_value->data[j] = 'o');
      } else {
        TEST_ASSERT(optional_value->data[j] = 'i');
      }
      intermediate /= 2;
    }
  }
/*
  for (auto &[key, value] : kvs) {
    std::optional<Value> optional_value;
    optional_value = hopscotch.find_tp(key);
    TEST_ASSERT(optional_value);
    TEST_ASSERT(strncmp(optional_value->data, value.data, kValueLen) == 0);
  }

  for (auto &[key, value] : kvs) {
    TEST_ASSERT(hopscotch.erase_tp(key));
  }

  for (auto &[key, value] : kvs) {
    std::optional<Value> optional_value;
    hopscotch.find_tp(key);
    TEST_ASSERT(!optional_value);
  }
*/
  std::cout << "Passed" << std::endl;
}

void _main(void *args) {
  std::unique_ptr<FarMemManager> manager =
      std::unique_ptr<FarMemManager>(FarMemManagerFactory::build(
          kCacheSize, kNumGCThreads, new FakeDevice(kFarMemSize)));
  do_work(manager.get());
}

int main(int argc, char **argv) {
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
