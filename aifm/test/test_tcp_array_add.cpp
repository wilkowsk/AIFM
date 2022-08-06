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
constexpr uint64_t kFarMemSize = (4ULL << 32);
constexpr uint32_t kNumGCThreads = 12;
constexpr uint64_t kNumElements = 10000;

#define ROOT_NODE_ID 0
#define DO_PRINT false
#define USER_INPUT true

#define TOP_DOWN_TRAV true
#define BOTTOM_UP_TRAV true
#define HYBRID_TRAV true

#define INPUT_NUM_NODES 8388608LL*2
#define INPUT_RATIO 64
#define INPUT_NUM_EDGES INPUT_NUM_NODES*INPUT_RATIO

#define SCOPE_COUNTER 32768LL

#define INPUT_ARRAY_SIZE 8388608/8

#include <chrono>
using namespace std::chrono;

struct data1k {
  int data[1024];
};

namespace far_memory {

class FarMemTest {
public:
  long int graphLoc(
  DataFrameVector<int>* graph,
  long int x,
  long int y,
  long int n,
  long int e) {
    // return (x+(6*y));
    long int loc = 0;                             // initialize loc to 0
    //int numNodes = graph->at_mut(*scope, 1); // numNodes is at (1,0)
    //int numEdges = graph->at_mut(*scope, 0); // numEdges is at (0,0)
    long int numNodes = n;
    long int numEdges = e;
    switch(x) {
      case 5:
        loc += numNodes;    // total: loc = (2*numNodes) + numEdges + 2
      case 4:
        loc += numEdges;    // total: loc = numNodes + numEdges + 2
      case 3:
        loc += numNodes;    // total: loc = numNodes + 2
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

  long int graphAt(
  DataFrameVector<int>* graph,
  long int x,
  long int y,
  long int n,
  long int e,
  DerefScope* scope = NULL) {
    long int funcOutput;
    if (scope == NULL) {
      DerefScope backupScope;
      funcOutput = graph->at_mut(backupScope, graphLoc(graph, x, y, n, e)); // return the element at (x,y)
    } else {
      funcOutput = graph->at_mut(*scope, graphLoc(graph, x, y, n, e)); // return the element at (x,y)
    }
    return funcOutput;
  }

  void graphSet(
  DataFrameVector<int>* graph,
  long int x,
  long int y,
  long int n,
  long int e,
  long int target,
  DerefScope* scope = NULL) {
    if (scope == NULL) {
      DerefScope backupScope;
      graph->at_mut(backupScope, graphLoc(graph, x, y, n, e)) = target; // set the element at (x,y) to target
    } else {
      graph->at_mut(*scope, graphLoc(graph, x, y, n, e)) = target; // set the element at (x,y) to target
    }
  }

  void gtepsDisplay(
  std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1, 1000000000> > > start,
  long int traversedEdges,
  const char* name) {
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "{{{{ " << name << " traversal }}}}" << std::endl;
    std::cout << "{{{{ " << duration.count() << " ns }}}}" << std::endl;
    std::cout << "{{{{ " << traversedEdges << " edges }}}}" << std::endl;
    std::cout << "{{{{ " << (traversedEdges * 1.0) / duration.count() << " gteps }}}}" << std::endl;
  }

  int topDownStep(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  DataFrameVector<int>* frontier,
  int* count,
  DataFrameVector<int>* newFrontier) {
    int numEdges = graphAt(graph, 0, 0, -1, -1);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int newCount = 0;      // we're working on a new layer, so it's empty for now
    newFrontier->clear();  // therefore, clear the new frontier
    {
      int savedJ = 0;
      int edgeOffset = 0;
      for (int i = 0; i < *count; i = savedJ) {           // for each node in the frontier:
        DerefScope scope;
        int scopeUsage = 0;
        for (int j = i; (j < *count) && (scopeUsage < SCOPE_COUNTER); j++) {
          int node;
          {
            node = frontier->at(scope, j);          // establish the node
          }
          if (DO_PRINT) {
            std::cout << "Working on node " << node << std::endl;
          }
          int startEdge = graphAt(graph, 2, node, numNodes, numEdges, &scope);     // establish the edge range: start
          startEdge += edgeOffset;
          edgeOffset = 0;
          int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, 2, node + 1, numNodes, numEdges, &scope); // establish the edge range: end
          for (int neighbor = startEdge; neighbor < endEdge; neighbor++) {        // for each edge:
            scopeUsage++;
            if (scopeUsage >= SCOPE_COUNTER) {
              edgeOffset += neighbor - graphAt(graph, 2, node, numNodes, numEdges, &scope);
              savedJ--;
              //std::cout << "{" << savedJ << "/" << edgeOffset << "} "; // 
              break;
            }
            int outgoing = graphAt(graph, 3, neighbor, numNodes, numEdges, &scope);  // establish the neighboring nodes
            traversedEdges++;
            {
              if (distances->at(scope, outgoing) == -1) {                            // when the node hasn't been visited yet,
                distances->at_mut(scope, outgoing) = distances->at(scope, node) + 1;  // write its distance
                newCount++;                                                           // increment the new count
                newFrontier->push_back(scope, outgoing);                             // add the node to the new frontier
                if (DO_PRINT) {
                  std::cout << "Linked node " << outgoing << std::endl;
                }
              }
            }
          }
          savedJ++;
        }
      }
    }
    if (DO_PRINT) {
      std::cout << "------" << std::endl;
    }
    *count = newCount;      // the counts are swapped, because we're going from one layer to the next.
    return traversedEdges;  // return # traversed edges.
  }

  void topDownTraverse(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    int numEdges = graphAt(graph, 0, 0, -1, -1);;   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);;   // numNodes is at (1,0)

    int traversedEdges = 0;
    
    for (int i = 0; i < numNodes; i += SCOPE_COUNTER) {  // initialise the distances of each node...
      DerefScope scope;
      for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
        if (j == ROOT_NODE_ID) {
          distances->push_back(scope, 0);     // to 0 on the root node...
        } else {
          distances->push_back(scope, -1);    // and sentinel -1 on all else.
        }
      }
    }

    int count = 1;    // the number of nodes on the previous layer.
    //int newCount = 0; // the number of nodes on this layer.

    auto vertexSet = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertexSet;                  // make a frontier of the last layer of nodes...
    {
      DerefScope scope;
      frontier->push_back(scope, ROOT_NODE_ID);    // currently containing only the root
    }
      
    auto vertexSet2 = manager->allocate_dataframe_vector<int>();
    auto newFrontier = &vertexSet2;               // make a frontier of the current layer of nodes
                                                  // this will remain empty for now.

    if (DO_PRINT) {
      std::cout << std::endl;
    }
    //std::cout << "ZONECK TD 1 ";
    while (count != 0) {  // while there are still nodes to sort through...
      traversedEdges += topDownStep(manager, distances, graph, frontier, &count, newFrontier);
      //std::cout << "ZONECK TD 2 {" << traversedEdges << "} ";
      {
        auto temp = frontier;   // swap the frontiers -- we're moving to the next layer now
        frontier = newFrontier;
        newFrontier = temp;
      }
    }
    gtepsDisplay(start, traversedEdges, "Top-down");
  }

  int bottomUpStep(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  int offset,
  DataFrameVector<int>* frontier,
  int* count,
  int* iterator) {
    int numEdges = graphAt(graph, 0, 0, -1, -1);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int downCounter = 0;

    int newCount = *count;
    {
      int savedJ = *count - 1;
      int edgeOffset = 0;
      for (int i = *count - 1; i >= 0; i = savedJ) {       // for each node...
        DerefScope scope;
        int scopeUsage = 0;
        for (int j = i; j >= 0 && (scopeUsage < SCOPE_COUNTER); j--) {
          int node;
          {
            node = frontier->at(scope, j + offset);  // establish the node
          }
          if (DO_PRINT) {
            std::cout << "Working on node " << node << std::endl;
          }
          int startEdge = graphAt(graph, 4, node, numNodes, numEdges, &scope); // establish the edge range: start
          startEdge += edgeOffset;
          edgeOffset = 0;
          int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, 4, node + 1, numNodes, numEdges, &scope);
          for (int neighbor = startEdge; neighbor < endEdge; neighbor++) {        // for each edge:
            scopeUsage += 4;
            if (scopeUsage >= SCOPE_COUNTER) {
              edgeOffset = neighbor - graphAt(graph, 4, node, numNodes, numEdges, &scope);
              savedJ++;
              //std::cout << "{" << savedJ << "/" << edgeOffset << "} "; // 
              break;
            }
            int outgoing = graphAt(graph, 5, neighbor, numNodes, numEdges, &scope);  // establish the neighboring nodes
            traversedEdges++;
            {
              if (distances->at(scope, outgoing) == *iterator) { // when the node is on the frontier,
                distances->at_mut(scope, node) = *iterator + 1;    // write its distance
                newCount--;                                         // decrement the new count
                {
                  //auto temp = frontier->at(scope, i);              // move the node to the end of the list
                  frontier->at_mut(scope, j + offset) = frontier->at(scope, newCount + offset);
                  //frontier->at_mut(scope, newCount) = temp;
                }
                frontier->pop_back(scope);                         // and delete it
                if (DO_PRINT) {
                  std::cout << "Linked to node " << outgoing << std::endl;
                }
                break;
              }
            }
          }
          savedJ--;
          if (savedJ != downCounter - 1) {
            //std::cout << "{ " << savedJ << "|" << downCounter << "} ";
          }
          downCounter = savedJ;
          //std::cout << j << " ";
        }
      }
    }
    if (DO_PRINT) {
      std::cout << "------" << std::endl;
    }
    if (newCount == *count) {
      std::cout << "This graph is not connected. ERROR." << std::endl;
      //break;
      *count = 0;             // break from infinite loop. todo: does this work? 
    } else {
      *count = newCount;      // update the count
    }
    return traversedEdges;  // return # traversed edges
  }

  void bottomUpTraverse(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    
    int numEdges = graphAt(graph, 0, 0, -1, -1);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);   // numNodes is at (1,0)

    int traversedEdges = 0;

    auto vertexSet = manager->allocate_dataframe_vector<int>();
    auto frontier = &vertexSet;                   // make a frontier of the unvisited nodes
                                                  // we'll be initialising this later

    int count = numNodes - 1;    // the number of unvisited nodes (updated after each iteration)
    //int newCount = count;        // the number of unvisited nodes (updated more frequently)
    int iterator = 0;            // the iteration count.

    for (int i = 0; i < numNodes; i += SCOPE_COUNTER) {      // initialise the distances of each node...
      DerefScope scope;
      for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
        if (j == ROOT_NODE_ID) {
          distances->push_back(scope, 0);         // to 0 for the root node...
        } else {
          distances->push_back(scope, -1);        // and sentinel -1 for all else
        }
        if (j != ROOT_NODE_ID) {        // except for the root...
          frontier->push_back(scope, j);  // put every node into the "unvisited" frontier
        }
      }
    }

    while (count != 0) {  // while there are still nodes to sort through...
      traversedEdges += bottomUpStep(manager, distances, graph, 0, frontier, &count, &iterator);
      std::cout << count;
      iterator++;           // update the iterator
    }
    
    //auto vertexSet2 = manager->allocate_dataframe_vector<int>(); // TODO: are these unnecessary?
    //auto newFrontier = &vertexSet2;
    gtepsDisplay(start, traversedEdges, "Bottom-up");
    return;
  }

  int hybridStep(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  int offset,
  DataFrameVector<int>* frontier,
  int* count) {
    int numEdges = graphAt(graph, 0, 0, -1, -1);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int newCount = 0;                               // we're working on a new layer, so it's empty for now
    {
      for (int i = 0; i < *count; i += SCOPE_COUNTER) {            // for each node in the frontier:
        DerefScope scope;
        for (int j = i; (j < *count) && (j < i+SCOPE_COUNTER); j++) {
          int node;
          {
            node = frontier->at(scope, j + offset);  // establish the node
          }
          if (DO_PRINT) {
            std::cout << "Working on node " << node << std::endl;
          }
          int startEdge = graphAt(graph, 2, node, numNodes, numEdges, &scope);     // establish the edge range: start
          int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, 2, node + 1, numNodes, numEdges, &scope); // establish the edge range: end
          for (int neighbor = startEdge; neighbor < endEdge; neighbor++) {        // for each edge:
            int outgoing = graphAt(graph, 3, neighbor, numNodes, numEdges, &scope);  // establish the neighboring nodes
            traversedEdges++;
            {
              if (distances->at(scope, outgoing) == -1) {                            // when the node hasn't been visited yet,
                distances->at_mut(scope, outgoing) = distances->at(scope, node) + 1;  // write its distance
                {
                                                                            // now i'm gonna warn you, what's about to happen is a bit odd
                  int backTraverser = outgoing;                             // make a traverser into the list of nodes...
                  int frontTraverser = frontier->at(scope, backTraverser); // starting at the index corresponding to the number of the target node.
                  while (frontTraverser != outgoing) {                      // while we haven't arrived at the target node yet...
                    //std::cout << backTraverser << " < ";
                    backTraverser = frontTraverser;
                    frontTraverser = frontier->at(scope, backTraverser);     // advance the traverser into the index corresponding to the number of the current node.
                  }                                                         // this while loop is guaranteed to halt, and runs in amortized constant time.
                  //std::cout << backTraverser << std::endl;
                  frontier->at_mut(scope, backTraverser) = frontier->at(scope, offset + *count + newCount); // copy from the first untraversed node in the vector to the location of the target node
                  frontier->at_mut(scope, offset + *count + newCount) = outgoing;  // copy the target node to replace the first untraversed node in the vector
                  if (DO_PRINT) {
                    for (int k = 0; k < numNodes; k++) {
                      if (k < 10) {
                        std::cout << 0;
                      }
                      std::cout << k;
                      if (k == offset + *count - 1 || k == offset + *count + newCount) {
                        std::cout << "|";
                      } else {
                        std::cout << " ";
                      }
                    }
                    std::cout << std::endl;

                    for (int k = 0; k < numNodes; i++) {
                      if (frontier->at(scope, k) < 10) {
                        std::cout << 0;
                      }
                      std::cout << frontier->at(scope, k);
                      if (k == offset + *count - 1 || k == offset + *count + newCount) {
                        std::cout << "|";
                      } else {
                        std::cout << " ";
                      }
                    }
                    std::cout << std::endl;
                  }
                }
                newCount++; // increment the new count
                if (DO_PRINT) {
                  std::cout << "Linked node " << outgoing << std::endl;
                }
              }
            }
          }
        }
      }
    }
    if (DO_PRINT) {
      std::cout << "------" << std::endl;
    }
    *count = newCount;      // the counts are swapped as well, because we're moving from one layer to the next.
    return traversedEdges;  // return # traversed edges
  }

  void hybridTraverse (
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    
    int numEdges = graphAt(graph, 0, 0, -1, -1);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, 1, 0, -1, -1);   // numNodes is at (1,0)

    int traversedEdges = 0;

    for (int i = 0; i < numNodes; i += SCOPE_COUNTER) {      // initialise the distances of each node...
      DerefScope scope;
      for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
        if (j == ROOT_NODE_ID) {
          distances->push_back(scope, 0);         // to 0 for the root node...
        } else {
          distances->push_back(scope, -1);        // and sentinel -1 for all else
        }
      }
    }
    auto bottomVertexSet = manager->allocate_dataframe_vector<int>(); // for petsdirbyh
    auto bottomFrontier = &bottomVertexSet;   // make a frontier of the unvisited nodes
                                              // we'll be initialising this later

    int bottomCount = 0;              // the number of unvisited nodes (updated after each iteration)
    //int bottomNewCount = bottomCount; // the number of unvisited nodes (updated more frequently)
    int iterator = 0;                 // the iteration count.

    int count = 1;    // the number of nodes on the previous layer.
    //int newCount = 0; // the number of nodes on this layer.

    auto vertexSet = manager->allocate_dataframe_vector<int>();
    auto frontier = &vertexSet;                 // make a frontier of the last layer of nodes...
    {
      DerefScope scope;
      frontier->push_back(scope, ROOT_NODE_ID);   // the root node goes first...
    }
    for (int i = 1; i < numNodes; i += SCOPE_COUNTER) {
      DerefScope scope;
      for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
        if (j == ROOT_NODE_ID) {
          frontier->push_back(scope, 0);              // other nodes get put in their own index...
        } else {
          frontier->push_back(scope, j);              // except for the zero node, whose index is the root node's number.
        }
      }
    }
      
    auto vertexSet2 = manager->allocate_dataframe_vector<int>(); // for petsdirbyh
    auto newFrontier = &vertexSet2;               // make a frontier of the current layer of nodes
                                                  // this will remain empty for now.

    int visitedNodes = 1;
    // bool bottomActivated = false;
    int bottomOffset;
    while (visitedNodes < ((numNodes * 1.0) * numNodes / numEdges)) { // while the percentage of visited nodes is less than the inverse of the average degree...
      //std::cout << "Top-down." << std::endl;
      //traversedEdges += hybridStep(manager, distances, graph, visitedNodes - count, frontier, &count); // for hybridstep
      traversedEdges += topDownStep(manager, distances, graph, frontier, &count, newFrontier);
      {
        auto temp = frontier;   // swap the frontiers -- we're moving to the next layer now
        frontier = newFrontier;
        newFrontier = temp;
      }
      visitedNodes += count;                                          // update the number of visited odes
      if (count == 0) {
        break;  // to prevent infinite loops
      }
      {
        //auto temp = frontier;       // swap the frontiers -- we're moving to the next layer now
        //frontier = newFrontier;
        //newFrontier = temp;
      }
      iterator++;
    }
    //bottomOffset = visitedNodes;            // the remaining nodes are at the end of the data structure, starting at index visitedNodes...
    bottomCount = numNodes - visitedNodes;  // and having a size of numNodes - visitedNodes
    for (int i = 1; i < numNodes; i += SCOPE_COUNTER) {
      DerefScope scope;
      for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
        if (distances->at(scope, j) == -1) {
          bottomFrontier->push_back(scope, j);
        }
      }
    }
    while (visitedNodes < numNodes) {                               // while we've visited fewer nodes than we have...
      traversedEdges += bottomUpStep(manager, distances, graph, 0, bottomFrontier, &bottomCount, &iterator);
      visitedNodes = numNodes - bottomCount;    // update the number of visited nodes
      // std::cout << bottomCount << std::endl;
      iterator++;
    }
    gtepsDisplay(start, traversedEdges, "Hybrid");
    return;
  }

  void do_work(FarMemManager *manager) {
    std::cout << "Running " << __FILE__ "..." << std::endl;

    auto graph = manager->allocate_dataframe_vector<int>();       // allocate our vectors
    auto distances = manager->allocate_dataframe_vector<int>();
    auto distances2 = manager->allocate_dataframe_vector<int>();
    auto distances3 = manager->allocate_dataframe_vector<int>();
    auto usedNumbers = manager->allocate_dataframe_vector<int>();
    int numEdges;                                                 // don't forget our ints, allocate them too!
    int numNodes;

    if (USER_INPUT == true) { // input subsection.
      long unsigned int inputVar = 0;
      long int phaseI = 0;     // main phases
      long int phaseJ = 0;     // subphases
      long int inputNodes = 0; // # nodes
      long int inputEdges = 0; // # edges
      while (inputVar != -2) {
        bool refresh = false;
        DerefScope scope;
        while (inputVar != -2 && !refresh) {
          if (phaseJ % SCOPE_COUNTER == 0 && false) {
            switch (phaseI) {
              case 0: // constants
                switch (phaseJ) {
                  case 0: // graph header token
                    std::cout << "Recite the magic value" << std::endl; // to mimick the graph header token
                    std::cout << 0xDEADBEEF << std::endl;               // this converts to a decimal number in cout.
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
                std::cout << "Enter node #" << phaseJ << std::endl;
                break;
              case 2: // edge destinations, for each edge
                std::cout << "Enter edge #" << phaseJ << std::endl;
                break;
              default: // it should never get here.
                std::cout << "ERROR" << std::endl;
            }
          }
          {
            //std::cin >> inputVar; // input. replace with file reading on port.
            switch (phaseI) {
              case 0:
                switch (phaseJ) {
                  case 0:
                    inputVar = 0xDEADBEEF;
                    break;
                  case 1:
                    inputVar = INPUT_NUM_NODES;
                    break;
                  case 2:
                    inputVar = INPUT_NUM_EDGES;
                    break;
                }
                break;
              case 1:
                inputVar = INPUT_RATIO*phaseJ;
                inputVar -= (phaseJ % 100800);
                break;
              case 2:
                inputVar = (phaseJ + 3) % INPUT_NUM_NODES;
                break;
            }
          }
          switch (phaseI) {
            case 0: // constants
              switch (phaseJ) {
                case 0: // graph header token
                  if (inputVar != 0xDEADBEEF) {       // match the magic value...
                    std::cout << "ERROR" << std::endl;  // or else you'll get an error
                  }
                  phaseJ++;                           // there's another subphase...
                  break;
                case 1: // number of nodes
                  inputNodes = inputVar;        // we've inputted the number of nodes.
                  for (int k = 0; k < 2; k++) {
                    graph.push_back(scope, -1);   // expand the vector to accomodate the new constants
                  }
                  graphSet(&graph, 1, 0, -1, -1, inputNodes, &scope); // numNodes is at (1,0), set it to inputNodes
                  phaseJ++;                     // there's another subphase...
                  break;
                case 2: // number of edges
                  inputEdges = inputVar;          // we've inputted the number of edges.
                  graphSet(&graph, 0, 0, -1, -1, inputEdges, &scope); // numEdges is at (0,0), set it to inputEdges
                  phaseI++;   // move to phase 1, 
                  phaseJ = 0; // subphase 0.
                  refresh = true;
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
              graphSet(&graph, 2, phaseJ, inputNodes, inputEdges, inputVar, &scope); // put the node offset into the vector
              phaseJ++;                     // there's another subphase...
              if (phaseJ >= inputNodes) {   // unless there isn't...
                phaseI = 2;                   // then, move to phase 2...
                phaseJ = 0;                   // subphase 0.
                refresh = true;
              }
              if (phaseJ % SCOPE_COUNTER == 0) {
                refresh = true;
              }
              break;
            case 2: // edge offsets, for each edge
              for (int k = 0; k < 2; k++) {
                graph.push_back(scope, -1);   // expand the vector to accomodate the new edge
              }
              graphSet(&graph, 3, phaseJ, inputNodes, inputEdges, inputVar, &scope); // put the edge destination into the vector.
              phaseJ++;                       // there's another subphase...
              if (phaseJ >= inputEdges) {     // unless there isn't...
                inputVar = -2;                  // then exit the loop.
              }
              if (phaseJ % SCOPE_COUNTER == 0) {
                refresh = true;
              }
              break;
            default: // it should never get here.
              std::cout << "ERROR" << std::endl;
              inputVar = -2; // to exit infinite loops.
              break;
          }
        }
      } 
    }
    std::cout << "ZONECK 1 ";
    for (int repeats = 0; repeats < 1; repeats++) { // for testing purposes. set to repeats<1 for normal use.
      {

        //graph.clear();          // clear the vectors from the previous iteration
        distances.clear();
        distances2.clear();
        distances3.clear();
        usedNumbers.clear();
        int arraySize = INPUT_ARRAY_SIZE; // establish the size of the array
        std::random_device rd;  // set up the RNG
        std::mt19937_64 eng(rd());
        std::uniform_int_distribution<uint64_t> distr(0, arraySize - 1); // to select a random element of the array.

        { // - 
          //int edges[arraySize][arraySize];      // each entry indicates whether an edge exists between two nodes.
          //for (int i = 0; i < arraySize; i++) {
          //  edges[i][i] = 0;
          //  for (int j = 0; j < arraySize; j++) { // for each distinct pair of nodes:
          //    int randomNum = distr(eng);
          //    if (randomNum != 1 || j == i) {
          //      randomNum = 0;                      // if the random num wasn't 1, set it to 0.
          //    }
          //    if (j == i + 1) {
          //      randomNum = 1;                      // if the two nodes are consecutive, add the edge anyway.
          //    }
          //    edges[i][j] = randomNum;
          //    // edges[j][i] = randomNum;
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
          //graph.at_mut(scope, 1).push_back(scope, arraySize);  // graph[1].push_back(scope, arraySize);
          //int runningOffset = 0;
          //for (int i = 0; i < arraySize; i++) {
          //  graph.at_mut(scope, 2).push_back(scope, runningOffset); // graph[2].push_back(scope, runningOffset);
          //  for (int j = 0; j < arraySize; j++) {
          //    if (edges[i][j] == 1) {
          //      graph.at_mut(scope, 3).push_back(scope, j);         // graph[3].push_back(scope, j);
          //      runningOffset++;
          //    }
          //  }
          //}
        }

        if (USER_INPUT == false) { // graphgen 
          int runningOffset = 0;  // keep track of how many edges we've gone through so far
          int outDegree = 4;      // set the outdegree of each node
          for (int i = 0; i < (2*arraySize*(1+outDegree))+2; i += SCOPE_COUNTER) {
            DerefScope scope;
            for (int j = i; (j < (2*arraySize*(1+outDegree))+2) && (j < i+SCOPE_COUNTER); j++) {
              graph.push_back(scope, -1);   // expand the vector to accomodate the new edge
            }
          }
          graphSet(&graph, 1, 0, -1, -1, arraySize);  // add the graph's size (in nodes and in edges) as data to the graph
          graphSet(&graph, 0, 0, -1, -1, arraySize*outDegree);
          int graphGenNodes = arraySize;
          int graphGenEdges = arraySize*outDegree;
          std::cout << "ZONECK 2 ";
          for (int i = 0; i < arraySize; i += SCOPE_COUNTER) { // for each node...
            DerefScope scope;
            for (int j = i; (j < arraySize) && (j < i+SCOPE_COUNTER); j++) {
              if (DO_PRINT) {
                std::cout << "Node " << j << " running offset: " << runningOffset << "     ";
                std::cout << "Linked nodes:";
              }
              int oldOffset = runningOffset;        // keep an old value of runningOffset, for later
              if (j == arraySize - 1) {     // todo: stop messing with the graphgen
                int temp = distr(eng);                // add an edge from the current node...
                usedNumbers.push_back(scope, temp);   // to a random node if we're at the last node
              } else {
                usedNumbers.push_back(scope, j + 1);  // to the next node otherwise
              }
              int currentOutDegree = 1;
              while (currentOutDegree < outDegree) {  // while we still haven't gotten all edges...
                int randomNum = distr(eng);             // select a random node
                bool isUnique = true;
                for (int k = 0; k < currentOutDegree; k++) {
                  if (usedNumbers.at(scope, k) == randomNum) {
                    isUnique = false;                   // try again if the edge there already exists
                  }
                }
                if (isUnique && randomNum != j) {             // if the edge isn't redundant...
                  usedNumbers.push_back(scope, randomNum);      // add it
                  currentOutDegree++;                           // count it
                }
              }
              for (int k = 0; k < outDegree; k++) {   // for each edge from this node...
                //int randomNum = distr(eng);
                //if (randomNum != 1 || j == i) {
                //  randomNum = 0;           // if the random num wasn't 1, set it to 0.
                //}
                //if (j == i + 1) {
                //  randomNum = 1;           // if the two nodes are consecutive, add the edge anyway.
                //}
                //if (randomNum == 1) {
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
                  std::cout << " " << usedNumbers.at(scope, k);
                }
                int tempInt;
                {
                  tempInt = usedNumbers.at(scope, k);
                }
                graphSet(&graph, 3, runningOffset, graphGenNodes, graphGenEdges, tempInt, &scope); // add the edge to the graph
                runningOffset++;                                                      // update the running offset
              }
              graphSet(&graph, 2, j, graphGenNodes, graphGenEdges, oldOffset, &scope);  // add the offset to the graph
              if (DO_PRINT) {
                std::cout << std::endl;
              }
              usedNumbers.clear();
            }
          }
        }

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
        std::cout << "ZONECK 3 ";
        numEdges = graphAt(&graph, 0, 0, -1, -1);   // numEdges is at (0,0)
        numNodes = graphAt(&graph, 1, 0, -1, -1);   // numNodes is at (1,0)

        { // for incoming edges
          auto nodeCounts = manager->allocate_dataframe_vector<int>();   // for the indegree of each node
          auto nodeScatter = manager->allocate_dataframe_vector<int>();  // to adjust the offset of each node
          for (int i = 0; i < numNodes; i += SCOPE_COUNTER) {   // at each node...
            DerefScope scope;
            for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
              nodeScatter.push_back(scope, 0);       // prefill the scatter and counts to 0
              nodeCounts.push_back(scope, 0);
            }
          }
          std::cout << "ZONECK 3.1 ";
          for (int i = 0; i < numNodes; i += (SCOPE_COUNTER/4)) {   // for each node...
            DerefScope scope;
            for (int j = i; (j < numNodes) && (j < i+(SCOPE_COUNTER/4)); j++) {
              int startEdge = graphAt(&graph, 2, j, numNodes, numEdges, &scope); // establish the edge range: start
              int endEdge = (j == numNodes - 1) ? numEdges : graphAt(&graph, 2, j + 1, numNodes, numEdges, &scope); // establish the node range: end
              for (int k = startEdge; k < endEdge; k++) {
                int targetNode = graphAt(&graph, 3, k, numNodes, numEdges, &scope); // establish the neighboring nodes
                {
                  nodeCounts.at_mut(scope, targetNode)++;                             // update the indegree of the target node
                }
              }
            }
          }
          std::cout << "ZONECK 3.2 ";
          graphSet(&graph, 4, 0, numNodes, numEdges, 0);    // the zeroth node has offset zero
          for (int i = 1; i < numNodes; i += SCOPE_COUNTER) { // for each other node...
            DerefScope scope;
            for (int j = i; (j < numNodes) && (j < i+SCOPE_COUNTER); j++) {
              int tempInt;
              {
                tempInt = nodeCounts.at_mut(scope, j - 1);
              }
              graphSet(&graph, 4, j, numNodes, numEdges, graphAt(&graph, 4, j - 1, numNodes, numEdges, &scope) + tempInt, &scope);
              // graph[4].push_back(scope, graph[4].at(scope, j - 1) + nodeCounts.at(scope, j - 1));
              // the offset is the sum of the previous offset and the indegree of the previous node
            }
          }
          std::cout << "ZONECK 3.3 ";
          {
            int savedJ = 0;
            int edgeOffset = 0;
            for (int i = 0; i < numNodes; i = savedJ) { // for each node...
              DerefScope scope;
              int scopeUsage = 0;
              for (int j = i; (j < numNodes) && (scopeUsage < SCOPE_COUNTER); j++) {
                //std::cout << "ZONECK 3.3.1 ";
                int startEdge = graphAt(&graph, 2, j, numNodes, numEdges, &scope); // establish the edge range: start
                startEdge += edgeOffset;
                edgeOffset = 0;
                int endEdge = (j == numNodes - 1) ? numEdges : graphAt(&graph, 2, j + 1, numNodes, numEdges, &scope);
                //std::cout << "ZONECK 3.3.2 ";
                for (int k = startEdge; k < endEdge; k++) {
                  scopeUsage += 4;
                  if (scopeUsage >= SCOPE_COUNTER) {
                    edgeOffset += k - graphAt(&graph, 2, j, numNodes, numEdges, &scope);
                    savedJ--;
                    //std::cout << "{" << savedJ << "/" << edgeOffset << "} "; // 
                    break;
                  }
                  //std::cout << "ZONECK 3.3.3 ";
                  int targetNode = graphAt(&graph, 3, k, numNodes, numEdges, &scope); // establish the neighboring nodes
                  int desiredElement = graphAt(&graph, 4, targetNode, numNodes, numEdges, &scope);   // start with the offset for the target node
                  //std::cout << "ZONECK 3.3.4 ";
                  {
                    desiredElement += nodeScatter.at_mut(scope, targetNode);      // and add the scatter for the target node
                    nodeScatter.at_mut(scope, targetNode)++;                      // and increment the scatter, to not overwrite the newly added edge
                  }
                  //std::cout << "ZONECK 3.3.5-{" << scopeUsage << "} ";
                  graphSet(&graph, 5, desiredElement, numNodes, numEdges, j, &scope);                // use that as an index to add the edge
                  //std::cout << "ZONECK 3.3.6 ";
                }
                savedJ++;
              }
            }
          }
        }
        std::cout << "ZONECK 4 ";
        { // graph print 
          if (DO_PRINT) {
            for (int i = 0; i < numNodes; i++) {
              int starting = graphAt(&graph, 2, i, numNodes, numEdges);
              int ending = (i == numNodes - 1) ? (numEdges) : (graphAt(&graph, 2, i + 1, numNodes, numEdges));  // TODO: use this syntax for other ending nodes
              // int ending = (i == numNodes - 1) ? (numEdges) : (graphAt(graph, scope, 2, i + 1));
              std::cout << "Info ~ Node " << i << " [[ " << starting << " ]] is linked to nodes:";                      // print debug info: node # and offset
              for (int j = starting; j < ending; j++) {
                std::cout << " " << graphAt(&graph, 3, j, numNodes, numEdges);                                    // print debug info: list of neighbouring nodes
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



          //for (int i = 0; i < numNodes; i++) {
          //  std::cout << distances.at_mut(scope, i);
          //}
          //std::cout << std::endl;
        }
      }
        if (TOP_DOWN_TRAV) {
          topDownTraverse(manager, &distances, &graph);   // call each BFS algorithm
        }
        if (BOTTOM_UP_TRAV) {
          bottomUpTraverse(manager, &distances2, &graph);
        }
        if (HYBRID_TRAV) {
          hybridTraverse(manager, &distances3, &graph);
        }
      if (TOP_DOWN_TRAV && BOTTOM_UP_TRAV && HYBRID_TRAV) { // compare 
        if (DO_PRINT) {
          for (int i = 0; i < numNodes; i++) {
            DerefScope scope;
            std::cout << distances.at_mut(scope, i);  // print debug info: topdown distances
          }
          std::cout << std::endl;
          for (int i = 0; i < numNodes; i++) {
            DerefScope scope;
            std::cout << distances2.at_mut(scope, i); // print debug info: bottomup distances
          }
          std::cout << std::endl;
          for (int i = 0; i < numNodes; i++) {
            DerefScope scope;
            std::cout << distances3.at_mut(scope, i); // print debug info: bottomup distances
          }
          std::cout << std::endl;
        }
        bool isOkay = true;
        for (int i = 0; i < numNodes; i++) { // for each node:
          DerefScope scope;
          if (distances.at_mut(scope, i) != distances2.at_mut(scope, i) || distances2.at_mut(scope, i) != distances3.at_mut(scope, i)) {
            std::cout << "X";   // throw error if the three distances don't agree
            isOkay = false;
            break;
          } else {
            if (DO_PRINT) {
              std::cout << "-"; // this dash indicates nothing's wrong with your distances
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
