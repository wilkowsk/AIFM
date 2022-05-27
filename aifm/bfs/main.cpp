/*
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string>
#include <getopt.h>

#include "CycleTimer.h"
#include "graph.h"
#include "bfs.h"


#define USE_BINARY_GRAPH 1


void usage(const char* progname) {
    printf("Usage: %s [options] graph_file_name\n", progname);
    printf("Program Options:\n");
    printf("  -t  --threads <N>  Use T threads\n");
    printf("  -?  --help         This message\n");
}

int main(int argc, char** argv) {

    int  num_threads = -1;
    std::string graph_filename;

    // parse commandline options ////////////////////////////////////////////
    int opt;
    static struct option long_options[] = {
        {"threads", 1, 0, 't'},
        {"help", 0, 0, '?'},
        {0 ,0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "t:?", long_options, NULL)) != EOF) {

        switch (opt) {
        case 't':
        {
            num_threads = atoi(optarg);
            break;
        }
        case '?':
        default:
            usage(argv[0]);
            return 1;
        }
    }
    // end parsing of commandline options


    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }

    graph_filename = argv[optind];

    graph g;

    //printf("----------------------------------------------------------\n");
    //printf("OMP Max system threads = %d\n", omp_get_max_threads());
    //if (num_threads > 0)
    //    omp_set_num_threads(num_threads);
    //printf("OMP will use at most %d threads.\n", omp_get_max_threads());
    //printf("----------------------------------------------------------\n");

    printf("Loading graph (this can take some time for the bigger graphs)...\n");
    load_graph_binary(graph_filename.c_str(), &g);

    printf("Graph stats:\n");
    printf("  Edges: %d\n", g.num_edges);
    printf("  Nodes: %d\n", g.num_nodes);

    solution sol1;
    sol1.distances = (int*)malloc(sizeof(int) * g.num_nodes); // sol1.distances = manager->allocate_array<int, g.num_nodes>();
    solution sol2;
    sol2.distances = (int*)malloc(sizeof(int) * g.num_nodes); // sol2.distances = manager->allocate_array<int, g.num_nodes>();
    solution sol3;
    sol3.distances = (int*)malloc(sizeof(int) * g.num_nodes); // sol3.distances = manager->allocate_array<int, g.num_nodes>();

    // execute bottom-up version. this can be unco,mmented.

    //start_time = CycleTimer::currentSeconds();
    //bfs_bottom_up(&g, &sol1);
    //end_time = CycleTimer::currentSeconds();
    //printf("Bottom up BFS time: %.3f sec\n", end_time - start_time);

    // execute top-down version

    double start_time = CycleTimer::currentSeconds();
    bfs_top_down(&g, &sol2);
    double end_time = CycleTimer::currentSeconds();
    printf("Top down BFS time: %.3f sec\n", end_time - start_time);

    // execute hybrid version

    //start_time = CycleTimer::currentSeconds();
    //bfs_hybrid(&g, &sol3);
    //end_time = CycleTimer::currentSeconds();
    //printf("Hybrid BFS time: %.3f sec\n", end_time - start_time);


    //for (int i=0; i<g.num_nodes; i++) {
    //  if (sol1.distances[i] != sol2.distances[i]) {
    //        fprintf(stderr, "*** Distance results disagree at node %d: %d, %d\n", i, sol1.distances[i], sol2.distances[i]);
    //        exit(1);
    //  }
   // }

    printf("Bottom-up and top-down distance results agree.\n");

    return 0;
}
*/
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
    } else if (middleElement < targetElement) { // middleElement is too low
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
  // add_array(&array_C, &array_A, &array_B);
  uint64_t result;

  for (uint64_t i = 0; i < 10000; i += 100) {
    DerefScope scope;
    result = search_array(&array_A, raw_array_A[i]);
    if (result != i) {
      goto fail;
    }
    result = search_array(&array_A, raw_array_A[i] + 1);
    if (result != -1ULL) {
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

