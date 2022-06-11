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

#include <random>

using namespace far_memory;

constexpr uint64_t kCacheSize = (128ULL << 20);
constexpr uint64_t kFarMemSize = (4ULL << 30);
constexpr uint32_t kNumGCThreads = 12;
constexpr uint64_t kNumElements = 10000;

#define ROOT_NODE_ID 0

namespace far_memory {

class FarMemTest {
public:
  template <uint64_t N, typename T>
  int graphAt(
  DataFrameVector<int>* graph,
  DerefScope scope,
  int x,
  int y) {
    return graph->at_mut(scope, x).at_mut(scope, y);
  }

  void graphSet(
  DataFrameVector<int>* graph,
  DerefScope scope,
  int x,
  int y,
  int target) {
    graph->at_mut(scope, x).at_mut(scope, y) = target;
  }

  template <uint64_t N, typename T>
  void top_down_sort(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  Array<T, N>* graph) {
    DerefScope scope;

    int num_edges = graph->at_mut(scope, 0).at_mut(scope, 0); // int num_edges = graphAt(graph, scope, 0, 0);
    int num_nodes = graph->at_mut(scope, 1).at_mut(scope, 0); // int num_nodes = graphAt(graph, scope, 1, 0);
    
    for (int i = 0; i < num_nodes; i++) {
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);
      } else {
        distances->push_back(scope, -1);
      }
    }

    int count = 1;
    int newCount = 0;

    auto vertex_set = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertex_set;
      frontier->push_back(scope, ROOT_NODE_ID);
      
    auto vertex_set2 = manager->allocate_dataframe_vector<int>();
    auto new_frontier = &vertex_set2;

    //std::cout << std::endl;
    while (count != 0) { // while we haven't run out of nodes:
      newCount = 0;
      new_frontier->clear(); // clear the new frontier
      {
        for (int i = 0; i < count; i++) {           // for each node in the frontier:
          int node = frontier->at(scope, i);          // establish the node
          //std::cout << "Working on node " << node << std::endl;
          int start_edge;
          {
            auto temp = &graph->at_mut(scope, 2); // start_edge = graphAt(graph, scope, 2, node);
            start_edge = temp->at(scope, node);       // establish the edge range: start
          }
          int end_edge;
          if (node == num_nodes - 1) {
            end_edge = num_edges;
          } else {
            auto temp = &graph->at_mut(scope, 2); // end_edge = graphAt(graph, scope, 2, node + 1);
            end_edge = temp->at(scope, node + 1);     // establish the edge range: end
          }
          for (int neighbor = start_edge; neighbor < end_edge; neighbor++) { // for each edge:
            int outgoing;
            {
              auto temp = &graph->at_mut(scope, 3); // outgoing = graphAt(graph, scope, 3, neighbor);
              outgoing = temp->at(scope, neighbor);   // establish the neighboring nodes
            }
            if (distances->at(scope, outgoing) == -1) {  // when the node hasn't been visited yet,
              distances->at_mut(scope, outgoing) = distances->at(scope, node) + 1;  // write its distance
              newCount++;                             // increment the new count
              new_frontier->push_back(scope, outgoing);  // add the node to the new frontier
              //std::cout << "Linked node " << outgoing << std::endl;
            }
          }
        }
      }
      //std::cout << "------" << std::endl;
      {
        auto temp = frontier;
        frontier = new_frontier;
        new_frontier = temp;
      }
      count = newCount;
    }
  }

  template <uint64_t N, typename T>
  void bottom_up_sort(
    FarMemManager *manager,
    DataFrameVector<int>* distances,
    Array<T, N>* graph) {
      DerefScope scope;

      int num_edges = graph->at_mut(scope, 0).at_mut(scope, 0); // int num_edges = graphAt(graph, scope, 0, 0);
      int num_nodes = graph->at_mut(scope, 1).at_mut(scope, 0); // int num_nodes = graphAt(graph, scope, 1, 0);

      for (int i = 0; i < num_nodes; i++) {
        if (i == ROOT_NODE_ID) {
          distances->push_back(scope, 0);
        } else {
          distances->push_back(scope, -1);
        }
      }
      auto vertex_set = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertex_set;

      int count = num_nodes - 1;
      int newCount = count;
      int iterator = 0;

      for (int i = 0; i < num_nodes; i++) {
        if (i != ROOT_NODE_ID) {
          frontier->push_back(scope, i);
        }
      }

      while (count != 0) {
        newCount = count;
        for (int i = count - 1; i >= 0; i--) {
          int node = frontier->at(scope, i);          // establish the node
          //std::cout << "Working on node " << node << std::endl;
          int start_edge;
          {
            auto temp = &graph->at_mut(scope, 4); // start_edge = graphAt(graph, scope, 4, node);
            start_edge = temp->at(scope, node);       // establish the edge range: start
          }
          int end_edge;
          if (node == num_nodes - 1) {
            end_edge = num_edges;
          } else {
            auto temp = &graph->at_mut(scope, 4); // end_edge = graphAt(graph, scope, 4, node + 1);
            end_edge = temp->at(scope, node + 1);     // establish the edge range: end
          }
          for (int neighbor = start_edge; neighbor < end_edge; neighbor++) { // for each edge:
            int outgoing;
            {
              auto temp = &graph->at_mut(scope, 5); // outgoing = graphAt(graph, scope, 5, neighbor);
              outgoing = temp->at(scope, neighbor);   // establish the neighboring nodes
            }
            if (distances->at(scope, outgoing) == iterator) {  // when the node is on the frontier,
              distances->at_mut(scope, node) = iterator + 1;  // write its distance
              newCount--;                             // decrement the new count
              {
                auto temp = frontier->at(scope, i);
                frontier->at_mut(scope, i) = frontier->at(scope, newCount);
                frontier->at_mut(scope, newCount) = temp;
              }
              frontier->pop_back(scope);
              //std::cout << "Linked to node " << outgoing << std::endl;
              break;
            }
          }
        }
        //std::cout << "------" << std::endl;
        if (newCount == count) {
          //std::cout << "This graph is not connected. ERROR." << std::endl;
          break;
        }
        count = newCount;
        iterator++;
      }
      
      auto vertex_set2 = manager->allocate_dataframe_vector<int>();
      auto new_frontier = &vertex_set2;
      return;
    }

  void do_work(FarMemManager *manager) {
    std::cout << "Running " << __FILE__ "..." << std::endl;


    //auto element = (manager->allocate_dataframe_vector<int>());
    //for (int i = 0; i < 2000; i++) {
    //  DerefScope scope;
    //  element.push_back(scope, i);
    //}
    //for (int i = 0; i < 2000; i++) {
    //  DerefScope scope;
    //  std::cout << element.at_mut(scope, 6) << std::endl;
    //}
    
    auto graph = manager->allocate_array<DataFrameVector<int>, 6>();
    auto distances = manager->allocate_dataframe_vector<int>();
    auto distances2 = manager->allocate_dataframe_vector<int>();
    int num_edges;
    int num_nodes;
    for (int repeats = 0; repeats < 1; repeats++) {
      {
        DerefScope scope;
        for (int i = 0; i < 6; i++) {
          graph.at_mut(scope, i).clear();
        }
        distances.clear();
        distances2.clear();

        int arraySize = 100;

        std::random_device rd;
        std::mt19937_64 eng(rd());
        std::uniform_int_distribution<uint64_t> distr(0, 75); // an edge has only a small chance to be present.

        //int edges[arraySize][arraySize];    // each entry indicates whether an edge exists between two nodes.
        //for (int i = 0; i < arraySize; i++) {
        //  edges[i][i] = 0;
        //  for (int j = 0; j < arraySize; j++) { // for each distinct pair of nodes:
        //    int random_num = distr(eng);
        //    if (random_num != 1 || j == i) {
        //      random_num = 0;           // if the random num wasn't 1, set it to 0.
        //    }
        //    if (j == i + 1) {
        //      random_num = 1;           // if the two nodes are consecutive, add the edge anyway.
        //    }
        //    edges[i][j] = random_num;
        //    // edges[j][i] = random_num;
        //  }
        //}
//
        //int countEdges = 0;
        //for (int i = 0; i < arraySize; i++) {
        //  for (int j = 0; j < arraySize; j++) {
        //    std::cout << edges[i][j];
        //    if (edges[i][j] == 1) {
        //    countEdges++;
        //    }
        //  }
//
        //  std::cout << std::endl;
        //}

        
        //graph.at_mut(scope, 0).push_back(scope, countEdges); // graph[0].push_back(scope, countEdges);
        //graph.at_mut(scope, 1).push_back(scope, arraySize); // graph[1].push_back(scope, arraySize);
        //int runningOffset = 0;
        //for (int i = 0; i < arraySize; i++) {
        //  graph.at_mut(scope, 2).push_back(scope, runningOffset); // graph[2].push_back(scope, runningOffset);
        //  for (int j = 0; j < arraySize; j++) {
        //    if (edges[i][j] == 1) {
        //      graph.at_mut(scope, 3).push_back(scope, j); // graph[3].push_back(scope, j);
        //      runningOffset++;
        //    }
        //  }
        //}
        
        graph.at_mut(scope, 1).push_back(scope, arraySize); // graph[1].push_back(scope, arraySize);
        int runningOffset = 0;
        for (int i = 0; i < arraySize; i++) {
          graph.at_mut(scope, 2).push_back(scope, runningOffset); // graph[2].push_back(scope, runningOffset);
          std::cout << "Node " << i << " running offset: " << runningOffset << "     ";
          std::cout << "Linked nodes:";
          for (int j = 0; j < arraySize; j++) {
            int random_num = distr(eng);
            if (random_num != 1 || j == i) {
              random_num = 0;           // if the random num wasn't 1, set it to 0.
            }
            if (j == i + 1) {
              random_num = 1;           // if the two nodes are consecutive, add the edge anyway.
            }
            if (random_num == 1) {
              std::cout << " " << j;
              graph.at_mut(scope, 3).push_back(scope, j); // graph[3].push_back(scope, j);
              runningOffset++;
            }
          }
          std::cout << std::endl;
        }
        graph.at_mut(scope, 0).push_back(scope, runningOffset); // graph[0].push_back(scope, runningOffset);
        

        //runningOffset = 0;
        //for (int i = 0; i < arraySize; i++) {
        //  graph.at_mut(scope, 4).push_back(scope, runningOffset);
        //  for (int j = 0; j < arraySize; j++) {
        //    if (edges[j][i] == 1) {
        //      graph.at_mut(scope, 5).push_back(scope, j);
        //      runningOffset++;
        //    }
        //  }
        //}

        //// 0: number of edges
        //// 1: number of nodes
        //// 2: outgoing indices/starts
        //// 3: outgoing destins/edges
        //// 4: incoming indices/starts
        //// 5: incoming destins/edges
        //graph.at_mut(scope, 1).push_back(scope, arraySize);
        //int runningOffset = 0;
        //for (int i = 0; i < arraySize; i++) {
        //  graph.at_mut(scope, 2).push_back(scope, runningOffset);
        //  graph.at_mut(scope, 3).push_back(scope, i + 1);
        //  runningOffset++;
        //  int arraye[2];
        //  for (int j = 0; j < 2; j++) {
        //    int potential = i;
        //    bool usePotential = false;
        //    while (!usePotential) {
        //      std::cout << "besenj";
        //      potential = distr(eng);
        //      usePotential = (potential != i) && (potential != i + 1);
        //      for (int k = 0; (k < j) && !usePotential; k++) {
        //        if (potential == arraye[k]) {
        //          usePotential = false;
        //        }
        //      }
        //    }
        //    graph.at_mut(scope, 3).push_back(scope, potential);
        //    arraye[j] = potential;
        //    runningOffset++;
        //  }
        //}
        //graph.at_mut(scope, 0).push_back(scope, runningOffset);

        
        num_edges = graph.at_mut(scope, 0).at_mut(scope, 0); // num_edges = graphAt(graph, scope, 0, 0);
        num_nodes = graph.at_mut(scope, 1).at_mut(scope, 0); // num_nodes = graphAt(graph, scope, 1, 0);

        { // for incoming edges
          for (int i = 0; i < num_edges; i++) {
            graph.at_mut(scope, 5).push_back(scope, -1); // graph[5].push_back(scope, -1);
          }
          auto node_counts = manager->allocate_dataframe_vector<int>();
          auto node_scatter = manager->allocate_dataframe_vector<int>();
          for (int i = 0; i < num_nodes; i++) {
            node_scatter.push_back(scope, 0);
            node_counts.push_back(scope, 0);
          }
          for (int i = 0; i < num_nodes; i++) {
            int start_edge;
            {
              auto temp = &graph.at_mut(scope, 2); // start_edge = graphAt(graph, scope, 2, i);
              start_edge = temp->at(scope, i);       // establish the edge range: start
            }
            int end_edge;
            if (i == num_nodes - 1) {
              end_edge = num_edges;
            } else {
              auto temp = &graph.at_mut(scope, 2); // end_edge = graphAt(graph, scope, 2, i + 1);
              end_edge = temp->at(scope, i + 1);     // establish the edge range: end
            }
            for (int j = start_edge; j < end_edge; j++) {
              int target_node;
              {
                auto temp = &graph.at_mut(scope, 3); // target_node = graphAt(graph, scope, 3, j);
                target_node = temp->at(scope, j);   // establish the neighboring nodes
              }
              node_counts.at_mut(scope, target_node)++;
            }
          }
          graph.at_mut(scope, 4).push_back(scope, 0); // graph[4].push_back(scope, 0);
          for (int i = 1; i < num_nodes; i++) {
            graph.at_mut(scope, 4).push_back(scope, graph.at_mut(scope, 4).at(scope, i - 1) + node_counts.at_mut(scope, i - 1));
            // graph[4].push_back(scope, graph[4].at(scope, i - 1) + node_counts.at(scope, i - 1));
          }
          for (int i = 0; i < num_nodes; i++) {
            int start_edge;
            {
              auto temp = &graph.at_mut(scope, 2); // start_edge = graphAt(graph, scope, 2, i);
              start_edge = temp->at(scope, i);       // establish the edge range: start
            }
            int end_edge;
            if (i == num_nodes - 1) {
              end_edge = num_edges;
            } else {
              auto temp = &graph.at_mut(scope, 2); // end_edge = graphAt(graph, scope, 2, i + 1);
              end_edge = temp->at(scope, i + 1);     // establish the edge range: end
            }
            for (int j = start_edge; j < end_edge; j++) {
              int target_node;
              {
                auto temp = &graph.at_mut(scope, 3); // target_node = graphAt(graph, scope, 3, j);
                target_node = temp->at(scope, j);   // establish the neighboring nodes
              }
              int desiredElement = graph.at_mut(scope, 4).at_mut(scope, target_node); // int desiredElement = graphAt(graph, scope, 4, target_node);
              desiredElement += node_scatter.at_mut(scope, target_node);
              graph.at_mut(scope, 5).at_mut(scope, desiredElement) = i; // graphSet(graph, scope, 5, desiredElement, i);
              node_scatter.at_mut(scope, target_node)++;
            }
          }
        }

        //std::cout << graph.at_mut(scope, 3).at(scope, 6) << std::endl;

        for (int i = 0; i < num_nodes; i++) {
          int starting = graph.at_mut(scope, 2).at(scope, i); // int starting = graphAt(graph, scope, 2, i);
          int ending = (i == num_nodes - 1) ? (num_edges) : (graph.at_mut(scope, 2).at(scope, i+1));
          // int ending = (i == num_nodes - 1) ? (num_edges) : (graphAt(graph, scope, 2, i + 1));
          std::cout << "Info ~ Node " << i << " [[ " << starting << " ]] is linked to nodes:";
          for (int j = starting; j < ending; j++) {
            std::cout << " " << graph.at_mut(scope, 3).at(scope, j); // std::cout << " " << graphAt(graph, scope, 3, j);
          }
          std::cout << std::endl;
        }

        //for (int i = 0; i < 5; i++) {
        //  int arraye[] = {0,4,6,7,8};
        //  graph.at_mut(scope, 2).push_back(scope, arraye[i]);
        //}
        //for (int i = 0; i < 5; i++) {
        //  int arraye[] = {0,2,3,5,7};
        //  graph.at_mut(scope, 4).push_back(scope, arraye[i]);
        //}
        //for (int i = 0; i < 8; i++) {
        //  int arraye[] = {1,2,3,4,2,3,0,0};
        //  graph.at_mut(scope, 3).push_back(scope, arraye[i]);
        //}
        //for (int i = 0; i < 8; i++) {
        //  int arraye[] = {2,3,0,0,1,0,1,0};
        //  graph.at_mut(scope, 5).push_back(scope, arraye[i]);
        //}
        //for (int i = 0; i < 6; i++) {
        //  int arraye[] = {1,1,5,8,5,8};
        //  for (int j = 0; j < arraye[i]; j++) {
        //    auto temp = &graph.at_mut(scope, i);
        //    std::cout << temp->at(scope, j);
        //  }
        //  std::cout << std::endl;
        //}
        

        //graph.at_mut(scope, 1).at_mut(scope, 0) = 9;



        //for (int i = 0; i < num_nodes; i++) {
        //  std::cout << distances.at_mut(scope, i);
        //}
        //std::cout << std::endl;
      }
        top_down_sort(manager, &distances, &graph);
        bottom_up_sort(manager, &distances2, &graph);
      {
        DerefScope scope;
        for (int i = 0; i < num_nodes; i++) {
          std::cout << distances.at_mut(scope, i);
        }
        std::cout << std::endl;
        for (int i = 0; i < num_nodes; i++) {
          std::cout << distances2.at_mut(scope, i);
        }
        std::cout << std::endl;
        bool isOkay = true;
        for (int i = 0; i < num_nodes; i++) {
          if (distances.at_mut(scope, i) != distances2.at_mut(scope, i)) {
            std::cout << "X";
            isOkay = false;
            break;
          } else {
            std::cout << "-";
          }
        }
        std::cout << std::endl;
        if (isOkay != true) {
          break;
        }
        std::cout << "Iteration " << repeats << " okay." << std::endl;
      }
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
