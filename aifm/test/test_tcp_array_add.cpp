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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#define GRAPH_HEADER_TOKEN 0xDEADBEEF

using namespace far_memory;

constexpr uint64_t kCacheSize = (128ULL << 20);
constexpr uint64_t kFarMemSize = (4ULL << 30);
constexpr uint32_t kNumGCThreads = 12;
constexpr uint64_t kNumElements = 10000;

#define ROOT_NODE_ID 0
#define DO_PRINT true

namespace far_memory {

class FarMemTest {
public:
  int graphLoc(
  DataFrameVector<int>* graph,
  DerefScope* scope,
  int x,
  int y) {
    // return (x+(6*y));
    int loc = 0;                              // initialize loc to 0
    int num_nodes = graph->at_mut(*scope, 1); // num_nodes is at (1,0)
    int num_edges = graph->at_mut(*scope, 0); // num_edges is at (0,0)
    switch(x) {
      case 5:
        loc += num_nodes;   // total: loc = (2*num_nodes) + num_edges + 2
      case 4:
        loc += num_edges;   // total: loc = num_nodes + num_edges + 2
      case 3:
        loc += num_nodes;   // total: loc = num_nodes + 2
      case 2:
        loc++;              // total: loc = 2;
      case 1:
        loc++;              // total: loc = 1;
      case 0:
                            // total: loc = 0;
      break;
    }
    loc += y;               // add y to the offset determined above
    //std::cout << " { " << x << " " << y << " " << loc << " } "; // test statement, can remove
    return loc;             // return loc
  }

  int graphAt(
  DataFrameVector<int>* graph,
  DerefScope* scope,
  int x,
  int y) {
    int funcOutput = graph->at_mut(*scope, graphLoc(graph, scope, x, y)); // return the element at (x,y)
    return funcOutput;
  }

  void graphSet(
  DataFrameVector<int>* graph,
  DerefScope* scope,
  int x,
  int y,
  int target) {
    graph->at_mut(*scope, graphLoc(graph, scope, x, y)) = target; // set the element at (x,y) to target
  }

  void top_down_sort(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    DerefScope scope;

    int num_edges = graphAt(graph, &scope, 0, 0);   // num_edges is at (0,0)
    int num_nodes = graphAt(graph, &scope, 1, 0);   // num_nodes is at (1,0)
    
    for (int i = 0; i < num_nodes; i++) { // initialise the distances of each node...
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);   // to 0 on the root node...
      } else {
        distances->push_back(scope, -1);  // and sentinel -1 on all else.
      }
    }

    int count = 1;    // the number of nodes on the previous layer.
    int newCount = 0; // the number of nodes on this layer.

    auto vertex_set = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertex_set;                  // make a frontier of the last layer of nodes...
      frontier->push_back(scope, ROOT_NODE_ID);     // currently containing only the root
      
    auto vertex_set2 = manager->allocate_dataframe_vector<int>();
    auto new_frontier = &vertex_set2;               // make a frontier of the current layer of nodes
                                                    // this will remain empty for now.

    if (DO_PRINT) {
      std::cout << std::endl;
    }
    while (count != 0) {  // while there are still nodes to sort through...
      newCount = 0;           // we're working on a new layer, so it's empty for now
      new_frontier->clear();  // therefore, clear the new frontier
      {
        for (int i = 0; i < count; i++) {           // for each node in the frontier:
          int node = frontier->at(scope, i);          // establish the node
          if (DO_PRINT) {
            std::cout << "Working on node " << node << std::endl;
          }
          int start_edge = graphAt(graph, &scope, 2, node); // establish the edge range: start
          int end_edge;
          if (node == num_nodes - 1) {
            end_edge = num_edges;
          } else {
            end_edge = graphAt(graph, &scope, 2, node + 1); // establish the edge range: end
          }
          for (int neighbor = start_edge; neighbor < end_edge; neighbor++) { // for each edge:
            int outgoing = graphAt(graph, &scope, 3, neighbor);               // establish the neighboring nodes
            if (distances->at(scope, outgoing) == -1) {                       // when the node hasn't been visited yet,
              distances->at_mut(scope, outgoing) = distances->at(scope, node) + 1;  // write its distance
              newCount++;                                                           // increment the new count
              new_frontier->push_back(scope, outgoing);                             // add the node to the new frontier
              if (DO_PRINT) {
                std::cout << "Linked node " << outgoing << std::endl;
              }
            }
          }
        }
      }
      if (DO_PRINT) {
        std::cout << "------" << std::endl;
      }
      {
        auto temp = frontier;       // swap the frontiers -- we're moving to the next layer now
        frontier = new_frontier;
        new_frontier = temp;
      }
      count = newCount;  // the counts are swapped as well, for the same reason.
    }
  }

  void bottom_up_sort(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    DerefScope scope;

    int num_edges = graphAt(graph, &scope, 0, 0);   // num_edges is at (0,0)
    int num_nodes = graphAt(graph, &scope, 1, 0);   // num_nodes is at (1,0)

    for (int i = 0; i < num_nodes; i++) {           // initialise the distances of each node...
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);             // to 0 for the root node...
      } else {
        distances->push_back(scope, -1);            // and sentinel -1 for all else
      }
    }
    auto vertex_set = manager->allocate_dataframe_vector<int>();
    auto frontier = &vertex_set;                    // make a frontier of the unvisited nodes
                                                    // we'll be initialising this later

    int count = num_nodes - 1;    // the number of unvisited nodes (updated after each iteration)
    int newCount = count;         // the number of unvisited nodes (updated more frequently)
    int iterator = 0;             // the iteration count.

    for (int i = 0; i < num_nodes; i++) {
      if (i != ROOT_NODE_ID) {                      // except for the root...
        frontier->push_back(scope, i);              // put every node into the "unvisited" frontier
      }
    }

    while (count != 0) {  // while there are still nodes to sort through...
      newCount = count;                           // TODO: test whether this line is unnecessary.
      for (int i = count - 1; i >= 0; i--) {      // for each node...
        int node = frontier->at(scope, i);          // establish the node
        if (DO_PRINT) {
          std::cout << "Working on node " << node << std::endl;
        }
        int start_edge = graphAt(graph, &scope, 4, node); // establish the edge range: start
        int end_edge;
        if (node == num_nodes - 1) {
          end_edge = num_edges;
        } else {
          end_edge = graphAt(graph, &scope, 4, node + 1); // establish the edge range: end
        }
        for (int neighbor = start_edge; neighbor < end_edge; neighbor++) { // for each edge:
          int outgoing = graphAt(graph, &scope, 5, neighbor); // establish the neighboring nodes
          if (distances->at(scope, outgoing) == iterator) {  // when the node is on the frontier,
            distances->at_mut(scope, node) = iterator + 1;      // write its distance
            newCount--;                                         // decrement the new count
            {
              auto temp = frontier->at(scope, i);   // move the node to the end of the list
              frontier->at_mut(scope, i) = frontier->at(scope, newCount);
              frontier->at_mut(scope, newCount) = temp;
            }
            frontier->pop_back(scope);              // and delete it
            if (DO_PRINT) {
              std::cout << "Linked to node " << outgoing << std::endl;
            }
            break;
          }
        }
      }
      if (DO_PRINT) {
        std::cout << "------" << std::endl;
      }
      if (newCount == count) {
        //std::cout << "This graph is not connected. ERROR." << std::endl;
        break;
      }
      count = newCount; // update the count
      iterator++;       // and the iterator
    }
    
    auto vertex_set2 = manager->allocate_dataframe_vector<int>(); // TODO: are these unnecessary?
    auto new_frontier = &vertex_set2;
    return;
  }

  void do_work(FarMemManager *manager) {
    std::cout << "Running " << __FILE__ "..." << std::endl;

    auto graph = manager->allocate_dataframe_vector<int>();       // allocate our vectors
    auto distances = manager->allocate_dataframe_vector<int>();
    auto distances2 = manager->allocate_dataframe_vector<int>();
    auto usedNumbers = manager->allocate_dataframe_vector<int>();
    int num_edges;                                                // don't forget our ints, allocate them too!
    int num_nodes;

    { // input subsection.
      unsigned int inputVar;
      int phase_i = 0;  // main phases
      int phase_j = 0;  // subphases
      int inputNodes = 0; // # nodes
      int inputEdges = 0; // # edges
      DerefScope scope;
      while (inputVar != -2) {
        switch (phase_i) {
          case 0: // constants
            switch (phase_j) {
              case 0: // graph header token
                std::cout << "Recite the magic value" << std::endl; // to mimick the graph header token
                std::cout << 0xDEADBEEF << std::endl; // this converts to a decimal number in cout.
                break;
              case 1: // number of nodes
                std::cout << "Enter the number of nodes" << std::endl;
                break;
              case 2: // number of edges
                std::cout << "Enter the number of edges" << std::endl;
                break;
              default: // it should never get here.
                std::cout << "ERROR" << std::endl;
            }
            break;
          case 1: // edge offsets, for each node
            std::cout << "Enter node #" << phase_j << std::endl;
            break;
          case 2: // edge destinations, for each edge
            std::cout << "Enter edge #" << phase_j << std::endl;
            break;
          default: // it should never get here.
            std::cout << "ERROR" << std::endl;
        }
        std::cin >> inputVar; // input. replace with file reading on port.
        switch (phase_i) {
          case 0: // constants
            switch (phase_j) {
              case 0: // graph header token
                if (inputVar != 0xDEADBEEF) { // match the magic value...
                  std::cout << "ERROR" << std::endl; // or else you'll get an error
                }
                phase_j++; // there's another subphase...
                break;
              case 1: // number of nodes
                inputNodes = inputVar;          // we've inputted the number of nodes.
                for (int k = 0; k < 2; k++) {
                  graph.push_back(scope, -1);   // expand the vector to accomodate the new constants
                }
                graphSet(&graph, &scope, 1, 0, inputNodes); // num_nodes is at (1,0), set it to inputNodes
                phase_j++; // there's another subphase...
                break;
              case 2: // number of edges
                inputEdges = inputVar;          // we've inputted the number of edges.
                graphSet(&graph, &scope, 0, 0, inputEdges); // num_edges is at (0,0), set it to inputEdges
                phase_i++;   // move to phase 1, 
                phase_j = 0;         // subphase 0.
                break;
              default: // it should never get here.
                std::cout << "ERROR" << std::endl;
                  inputVar = -2;  // to exit infinite loops.
                  break;
            }
            break;
          case 1: // node offsets, for each node
            for (int k = 0; k < 2; k++) {
              graph.push_back(scope, -1);   // expand the vector to accomodate the new node
            }
            graphSet(&graph, &scope, 2, phase_j, inputVar); // put the node offset into the vector
            phase_j++; // there's another subphase...
            if (phase_j >= inputNodes) { // unless there isn't...
              phase_i = 2;  // then, move to phase 2...
              phase_j = 0;                // subphase 0.
            }
            break;
          case 2: // edge offsets, for each edge
            for (int k = 0; k < 2; k++) {
              graph.push_back(scope, -1);   // expand the vector to accomodate the new edge
            }
            graphSet(&graph, &scope, 3, phase_j, inputVar); // put the edge destination into the vector.
            phase_j++; // there's another subphase...
            if (phase_j >= inputEdges) { // unless there isn't...
              inputVar = -2; // then exit the loop.
            }
            break;
            break;
          default: // it should never get here.
            std::cout << "ERROR" << std::endl;
            inputVar = -2; // to exit infinite loops.
            break;
        }
      }
    }
    
    for (int repeats = 0; repeats < 1; repeats++) { // for testing purposes. set to repeats<1 for normal use.
      {
        DerefScope scope;

        //graph.clear();          // clear the vectors from the previous iteration
        distances.clear();
        distances2.clear();
        usedNumbers.clear();
        int arraySize = 10000000; // establish the size of the array
        std::random_device rd;  // set up the RNG
        std::mt19937_64 eng(rd());
        std::uniform_int_distribution<uint64_t> distr(0, arraySize - 1); // an edge has only a small chance to be present.

        { // - 
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

          //int countEdges = 0;
          //for (int i = 0; i < arraySize; i++) {
          //  for (int j = 0; j < arraySize; j++) {
          //    std::cout << edges[i][j];
          //    if (edges[i][j] == 1) {
          //    countEdges++;
          //    }
          //  }

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
        }
/*
        { // graphgen 
          int runningOffset = 0;  // keep track of how many edges we've gone through so far
          int outDegree = 4;
          for (int k = 0; k < (2*arraySize*(1+outDegree))+2; k++) {
            graph.push_back(scope, -1);   // expand the vector to accomodate the new edge
          }
          graphSet(&graph, &scope, 1, 0, arraySize);      // add the graph's size (in nodes and in edges) to the graph
          graphSet(&graph, &scope, 0, 0, arraySize*outDegree);
          for (int i = 0; i < arraySize; i++) { // for each node...
            if (DO_PRINT) {
              std::cout << "Node " << i << " running offset: " << runningOffset << "     ";
              std::cout << "Linked nodes:";
            }
            int oldOffset = runningOffset;  // keep an old value, for later
            if (i == arraySize - 1) {
              int temp = distr(eng);                // add an edge to the current node
              usedNumbers.push_back(scope, temp);   // to a random node if we're at the last node
            } else {
              usedNumbers.push_back(scope, i + 1);  // to the next node otherwise
            }
            int currentOutDegree = 1;
            while (currentOutDegree < outDegree) {  // while we still haven't gotten all edges...
              int random_num = distr(eng);                    // select a random node
              bool isUnique = true;
              for (int j = 0; j < currentOutDegree; j++) {
                if (usedNumbers.at(scope, j) == random_num) {
                  isUnique = false;                           // try again if the edge there already exists
                }
              }
              if (isUnique && random_num != i) {              // if the edge isn't redundant...
                usedNumbers.push_back(scope, random_num);       // add it
                currentOutDegree++;                             // count it
              }
            }
            for (int j = 0; j < outDegree; j++) {   // for each edge from this node...
              //int random_num = distr(eng);
              //if (random_num != 1 || j == i) {
              //  random_num = 0;           // if the random num wasn't 1, set it to 0.
              //}
              //if (j == i + 1) {
              //  random_num = 1;           // if the two nodes are consecutive, add the edge anyway.
              //}
              //if (random_num == 1) {
              //  if (DO_PRINT) {
              //    std::cout << " " << j;
              //  }
              //  for (int k = 0; k < 6; k++) {
              //    graph.push_back(scope, -1);
              //  }
              //  graphSet(&graph, &scope, 3, runningOffset, j); // graph[3].push_back(scope, j);
              //  runningOffset++;
              //}
              if (DO_PRINT) {
                std::cout << " " << usedNumbers.at(scope, j);
              }
              graphSet(&graph, &scope, 3, runningOffset, usedNumbers.at(scope, j)); // add the edge to the graph
              runningOffset++;                                                      // update the running offset
            }
            graphSet(&graph, &scope, 2, i, oldOffset);  // add the offset to the graph
            if (DO_PRINT) {
              std::cout << std::endl;
            }
            usedNumbers.clear();
          }
        }
*/
        { // - 
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
        }

        num_edges = graphAt(&graph, &scope, 0, 0);  // num_edges at 0,0
        num_nodes = graphAt(&graph, &scope, 1, 0);  // num_nodes at 1,0

        { // for incoming edges
          auto node_counts = manager->allocate_dataframe_vector<int>();   // for the indegree of each node
          auto node_scatter = manager->allocate_dataframe_vector<int>();  // to adjust the offset of each node
          for (int i = 0; i < num_nodes; i++) {   // at each node...
            node_scatter.push_back(scope, 0);       // prefill the scatter and counts to 0
            node_counts.push_back(scope, 0);
          }
          for (int i = 0; i < num_nodes; i++) {   // for each node...
            int start_edge;
            {
              start_edge = graphAt(&graph, &scope, 2, i);
              // establish the edge range: start
            }
            int end_edge;
            if (i == num_nodes - 1) {
              end_edge = num_edges;
            } else {
              end_edge = graphAt(&graph, &scope, 2, i + 1);
              // establish the edge range: end
            }
            for (int j = start_edge; j < end_edge; j++) {
              int target_node;
              {
                target_node = graphAt(&graph, &scope, 3, j);
                // establish the neighboring nodes
              }
              node_counts.at_mut(scope, target_node)++; // update the indegree of the target node
            }
          }
          graphSet(&graph, &scope, 4, 0, 0);    // the zeroth node has offset zero
          for (int i = 1; i < num_nodes; i++) { // for each other node...
            graphSet(&graph, &scope, 4, i, graphAt(&graph, &scope, 4, i - 1) + node_counts.at_mut(scope, i - 1));
            // graph[4].push_back(scope, graph[4].at(scope, i - 1) + node_counts.at(scope, i - 1));
            // the offset is the sum of the previous offset and the indegree of the previous node
          }
          for (int i = 0; i < num_nodes; i++) { // for each node...
            int start_edge;
            {
              start_edge = graphAt(&graph, &scope, 2, i);
              // establish the edge range: start
            }
            int end_edge;
            if (i == num_nodes - 1) {
              end_edge = num_edges;
            } else {
              end_edge = graphAt(&graph, &scope, 2, i + 1);
              // establish the edge range: end
            }
            for (int j = start_edge; j < end_edge; j++) {
              int target_node;
              {
                target_node = graphAt(&graph, &scope, 3, j);
                // establish the neighboring nodes
              }
              int desiredElement = graphAt(&graph, &scope, 4, target_node);   // start with the offset for the target node
              desiredElement += node_scatter.at_mut(scope, target_node);      // and add the scatter for the target node
              graphSet(&graph, &scope, 5, desiredElement, i);                 // use that as an index to add the edge
              node_scatter.at_mut(scope, target_node)++;                      // and increment the scatter, to not overwrite it
            }
          }
        }

        { // graph print 
          if (DO_PRINT) {
            for (int i = 0; i < num_nodes; i++) {
              int starting = graphAt(&graph, &scope, 2, i);
              int ending = (i == num_nodes - 1) ? (num_edges) : (graphAt(&graph, &scope, 2, i + 1));  // TODO: use this syntax for other ending nodes
              // int ending = (i == num_nodes - 1) ? (num_edges) : (graphAt(graph, scope, 2, i + 1));
              std::cout << "Info ~ Node " << i << " [[ " << starting << " ]] is linked to nodes:";    // print debug info: node # and offset
              for (int j = starting; j < ending; j++) {
                std::cout << " " << graphAt(&graph, &scope, 3, j);                                    // print debug info: list of neighbouring nodes
              }
              std::cout << std::endl;
            }
          }
        }

        { // - 
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
      }
        top_down_sort(manager, &distances, &graph);   // call each BFS algorithm
        bottom_up_sort(manager, &distances2, &graph);
      { // compare 
        DerefScope scope;
        if (DO_PRINT) {
          for (int i = 0; i < num_nodes; i++) {
            std::cout << distances.at_mut(scope, i);  // print debug info: topdown distances
          }
          std::cout << std::endl;
          for (int i = 0; i < num_nodes; i++) {
            std::cout << distances2.at_mut(scope, i); // print debug info: bottomup distances
          }
          std::cout << std::endl;
        }
        bool isOkay = true;
        for (int i = 0; i < num_nodes; i++) {
          if (distances.at_mut(scope, i) != distances2.at_mut(scope, i)) {
            std::cout << "X"; // throw error if the distances don't agree
            isOkay = false;
            break;
          } else {
            if (DO_PRINT) {
              std::cout << "-";
            }
          }
        }
        if (DO_PRINT) {
          std::cout << std::endl;
        }
        if (isOkay != true) {
          break;
        }
        std::cout << "Iteration " << repeats + 1 << " okay." << std::endl;
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
