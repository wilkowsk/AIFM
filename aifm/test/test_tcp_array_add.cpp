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

#include <chrono>
using namespace std::chrono;

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
    int numNodes = graph->at_mut(*scope, 1); // numNodes is at (1,0)
    int numEdges = graph->at_mut(*scope, 0); // numEdges is at (0,0)
    switch(x) {
      case 5:
        loc += numNodes;   // total: loc = (2*numNodes) + numEdges + 2
      case 4:
        loc += numEdges;   // total: loc = numNodes + numEdges + 2
      case 3:
        loc += numNodes;   // total: loc = numNodes + 2
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

  int top_down_step(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  DataFrameVector<int>* frontier,
  int* count,
  DataFrameVector<int>* newFrontier,
  DerefScope* scope) {
    int numEdges = graphAt(graph, scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int newCount = 0;           // we're working on a new layer, so it's empty for now
    newFrontier->clear();  // therefore, clear the new frontier
    {
      for (int i = 0; i < *count; i++) {           // for each node in the frontier:
        int node = frontier->at(*scope, i);          // establish the node
        if (DO_PRINT) {
          std::cout << "Working on node " << node << std::endl;
        }
        int startEdge = graphAt(graph, scope, 2, node); // establish the edge range: start
        int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, scope, 2, node + 1);
        for (int neighbor = startEdge; neighbor < endEdge; neighbor++) { // for each edge:
          int outgoing = graphAt(graph, scope, 3, neighbor);               // establish the neighboring nodes
          traversedEdges++;
          if (distances->at(*scope, outgoing) == -1) {                       // when the node hasn't been visited yet,
            distances->at_mut(*scope, outgoing) = distances->at(*scope, node) + 1;  // write its distance
            newCount++;                                                           // increment the new count
            newFrontier->push_back(*scope, outgoing);                             // add the node to the new frontier
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
    *count = newCount;  // the counts are swapped as well, for the same reason.
    return traversedEdges;
  }

  void top_down_sort(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    DerefScope scope;

    int numEdges = graphAt(graph, &scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, &scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;
    
    for (int i = 0; i < numNodes; i++) { // initialise the distances of each node...
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);   // to 0 on the root node...
      } else {
        distances->push_back(scope, -1);  // and sentinel -1 on all else.
      }
    }

    int count = 1;    // the number of nodes on the previous layer.
    //int newCount = 0; // the number of nodes on this layer.

    auto vertexSet = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertexSet;                  // make a frontier of the last layer of nodes...
      frontier->push_back(scope, ROOT_NODE_ID);     // currently containing only the root
      
    auto vertexSet2 = manager->allocate_dataframe_vector<int>();
    auto newFrontier = &vertexSet2;               // make a frontier of the current layer of nodes
                                                    // this will remain empty for now.

    if (DO_PRINT) {
      std::cout << std::endl;
    }
    while (count != 0) {  // while there are still nodes to sort through...
      traversedEdges += top_down_step(manager, distances, graph, frontier, &count, newFrontier, &scope);
      {
        auto temp = frontier;       // swap the frontiers -- we're moving to the next layer now
        frontier = newFrontier;
        newFrontier = temp;
      }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "{{{{ Top-down traversal. }}}}" << std::endl;
    std::cout << "{{{{ " << duration.count() << " ns }}}}" << std::endl;
    std::cout << "{{{{ " << traversedEdges << " edges }}}}" << std::endl;
    std::cout << "{{{{ " << (traversedEdges * 1.0) / duration.count() << " gteps }}}}" << std::endl;
  }

  int bottom_up_step(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  int offset,
  DataFrameVector<int>* frontier,
  int* count,
  int* iterator,
  DerefScope* scope) {
    int numEdges = graphAt(graph, scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int newCount = *count;                           // TODO: test whether this line is unnecessary.
    for (int i = *count - 1; i >= 0; i--) {      // for each node...
      int node = frontier->at(*scope, i + offset);          // establish the node
      if (DO_PRINT) {
        std::cout << "Working on node " << node << std::endl;
      }
      int startEdge = graphAt(graph, scope, 4, node); // establish the edge range: start
      int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, scope, 4, node + 1);
      for (int neighbor = startEdge; neighbor < endEdge; neighbor++) { // for each edge:
        int outgoing = graphAt(graph, scope, 5, neighbor); // establish the neighboring nodes
        traversedEdges++;
        if (distances->at(*scope, outgoing) == *iterator) {  // when the node is on the frontier,
          distances->at_mut(*scope, node) = *iterator + 1;      // write its distance
          newCount--;                                         // decrement the new count
          {
            //auto temp = frontier->at(*scope, i);   // move the node to the end of the list
            frontier->at_mut(*scope, i + offset) = frontier->at(*scope, newCount + offset); // todo: make this more efficient
            //frontier->at_mut(*scope, newCount) = temp;
          }
          frontier->pop_back(*scope);              // and delete it
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
    if (newCount == *count) {
      //std::cout << "This graph is not connected. ERROR." << std::endl;
      //break;
      *count = 0;
    } else {
      *count = newCount; // update the count
    }
    return traversedEdges;
  }

  void bottom_up_sort(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    DerefScope scope;

    int numEdges = graphAt(graph, &scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, &scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;

    for (int i = 0; i < numNodes; i++) {           // initialise the distances of each node...
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);             // to 0 for the root node...
      } else {
        distances->push_back(scope, -1);            // and sentinel -1 for all else
      }
    }
    auto vertexSet = manager->allocate_dataframe_vector<int>();
    auto frontier = &vertexSet;                    // make a frontier of the unvisited nodes
                                                    // we'll be initialising this later

    int count = numNodes - 1;    // the number of unvisited nodes (updated after each iteration)
    //int newCount = count;         // the number of unvisited nodes (updated more frequently)
    int iterator = 0;             // the iteration count.

    for (int i = 0; i < numNodes; i++) {
      if (i != ROOT_NODE_ID) {                      // except for the root...
        frontier->push_back(scope, i);              // put every node into the "unvisited" frontier
      }
    }

    while (count != 0) {  // while there are still nodes to sort through...
      traversedEdges += bottom_up_step(manager, distances, graph, 0, frontier, &count, &iterator, &scope);
      iterator++;       // and the iterator
    }
    
    //auto vertexSet2 = manager->allocate_dataframe_vector<int>(); // TODO: are these unnecessary?
    //auto newFrontier = &vertexSet2;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "{{{{ Bottom-up traversal. }}}}" << std::endl;
    std::cout << "{{{{ " << duration.count() << " ns }}}}" << std::endl;
    std::cout << "{{{{ " << traversedEdges << " edges }}}}" << std::endl;
    std::cout << "{{{{ " << (traversedEdges * 1.0) / duration.count() << " gteps }}}}" << std::endl;
    return;
  }

  int hybrid_step(
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph,
  int offset,
  DataFrameVector<int>* frontier,
  int* count,
  DerefScope* scope) {
    int numEdges = graphAt(graph, scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;

    int newCount = 0;           // we're working on a new layer, so it's empty for now
    {
      for (int i = 0; i < *count; i++) {           // for each node in the frontier:
        int node = frontier->at(*scope, i + offset);          // establish the node
        if (DO_PRINT) {
          std::cout << "Working on node " << node << std::endl;
        }
        int startEdge = graphAt(graph, scope, 2, node); // establish the edge range: start
        int endEdge = (node == numNodes - 1) ? numEdges : graphAt(graph, scope, 2, node + 1);
        for (int neighbor = startEdge; neighbor < endEdge; neighbor++) { // for each edge:
          int outgoing = graphAt(graph, scope, 3, neighbor);               // establish the neighboring nodes
          traversedEdges++;
          if (distances->at(*scope, outgoing) == -1) {                       // when the node hasn't been visited yet,
            distances->at_mut(*scope, outgoing) = distances->at(*scope, node) + 1;  // write its distance
            {
                                                                        // now i'm gonna warn you, what's about to happen is a bit odd
              int backTraverser = outgoing;                             // make a traverser into the list of nodes...
              int frontTraverser = frontier->at(*scope, backTraverser); // starting at the index corresponding to the number of the target node.
              while (frontTraverser != outgoing) {                      // while we haven't arrived at the target node yet...
                //std::cout << backTraverser << " < ";
                backTraverser = frontTraverser;
                frontTraverser = frontier->at(*scope, backTraverser);     // advance the traverser into the index corresponding to the number of the current node.
              }                                                         // this while loop is guaranteed to halt, and runs in amortized constant time.
              //std::cout << backTraverser << std::endl;
              frontier->at_mut(*scope, backTraverser) = frontier->at(*scope, offset + *count + newCount); // copy from the first untraversed node in the vector to the location of the target node
              frontier->at_mut(*scope, offset + *count + newCount) = outgoing;  // copy the target node to replace the first untraversed node in the vector
              if (DO_PRINT) {
                for (int i = 0; i < numNodes; i++) {
                  if (i < 10) {
                    std::cout << 0;
                  }
                  std::cout << i;
                  if (i == offset + *count - 1 || i == offset + *count + newCount) {
                    std::cout << "|";
                  } else {
                    std::cout << " ";
                  }
                }
                std::cout << std::endl;
              
                for (int i = 0; i < numNodes; i++) {
                  if (frontier->at(*scope, i) < 10) {
                    std::cout << 0;
                  }
                  std::cout << frontier->at(*scope, i);
                  if (i == offset + *count - 1 || i == offset + *count + newCount) {
                    std::cout << "|";
                  } else {
                    std::cout << " ";
                  }
                }
                std::cout << std::endl;
              }
            }
            newCount++;                                                           // increment the new count
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
    *count = newCount;  // the counts are swapped as well, for the same reason.
    return traversedEdges;
  }

  void hybrid_sort (
  FarMemManager *manager,
  DataFrameVector<int>* distances,
  DataFrameVector<int>* graph) {
    auto start = high_resolution_clock::now();
    DerefScope scope;

    int numEdges = graphAt(graph, &scope, 0, 0);   // numEdges is at (0,0)
    int numNodes = graphAt(graph, &scope, 1, 0);   // numNodes is at (1,0)

    int traversedEdges = 0;

    for (int i = 0; i < numNodes; i++) {           // initialise the distances of each node...
      if (i == ROOT_NODE_ID) {
        distances->push_back(scope, 0);             // to 0 for the root node...
      } else {
        distances->push_back(scope, -1);            // and sentinel -1 for all else
      }
    }
    //auto bottomVertexSet = manager->allocate_dataframe_vector<int>();
    //auto bottomFrontier = &bottomVertexSet;                    // make a frontier of the unvisited nodes
                                                    // we'll be initialising this later

    int bottomCount = 0;    // the number of unvisited nodes (updated after each iteration)
    //int bottomNewCount = bottomCount;         // the number of unvisited nodes (updated more frequently)
    int iterator = 0;             // the iteration count.

    for (int i = 0; i < numNodes; i++) {
      if (i != ROOT_NODE_ID) {                      // except for the root...
                      // put every node into the "unvisited" frontier
      }
    }

    int count = 1;    // the number of nodes on the previous layer.
    //int newCount = 0; // the number of nodes on this layer.

    auto vertexSet = manager->allocate_dataframe_vector<int>();
      auto frontier = &vertexSet;                  // make a frontier of the last layer of nodes...
      frontier->push_back(scope, ROOT_NODE_ID);     // currently containing only the root
      for (int i = 1; i < numNodes; i++) {
        if (i == ROOT_NODE_ID) {
          frontier->push_back(scope, 0);
        } else {
          frontier->push_back(scope, i);
        }
      }
      
    //auto vertexSet2 = manager->allocate_dataframe_vector<int>();
    //auto newFrontier = &vertexSet2;               // make a frontier of the current layer of nodes
                                                    // this will remain empty for now.

    int visitedNodes = 1;
    bool bottomActivated = false;
    int bottomOffset;
    while (visitedNodes < numNodes) {
      //std::cout << "Level " << iterator << " forthcoming. Current nodes: " << visitedNodes << std::endl;
      if (visitedNodes < ((numNodes * 1.0) * numNodes / numEdges)) {
        //std::cout << "Top-down." << std::endl;
        traversedEdges += hybrid_step(manager, distances, graph, visitedNodes - count, frontier, &count, &scope);
        visitedNodes += count;
        {
          //auto temp = frontier;       // swap the frontiers -- we're moving to the next layer now
          //frontier = newFrontier;
          //newFrontier = temp;
        }
      } else {
        //std::cout << "Bottom-up." << std::endl;
        if (bottomActivated == false) {
          bottomActivated = true;
          //for (int i = 0; i < numNodes; i++) {
          //  if (distances->at_mut(scope, i) == -1) {
          //    bottomFrontier->push_back(scope, i);
          //    bottomCount++;
          //  }
          //}
          bottomOffset = visitedNodes;
          bottomCount = numNodes - visitedNodes;
        }
        traversedEdges += bottom_up_step(manager, distances, graph, bottomOffset, frontier, &bottomCount, &iterator, &scope);
        visitedNodes = numNodes - bottomCount;
      }
      iterator++;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "{{{{ Hybrid traversal. }}}}" << std::endl;
    std::cout << "{{{{ " << duration.count() << " ns }}}}" << std::endl;
    std::cout << "{{{{ " << traversedEdges << " edges }}}}" << std::endl;
    std::cout << "{{{{ " << (traversedEdges * 1.0) / duration.count() << " gteps }}}}" << std::endl;
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
/*
    { // input subsection.
      unsigned int inputVar;
      int phaseI = 0;     // main phases
      int phaseJ = 0;     // subphases
      int inputNodes = 0; // # nodes
      int inputEdges = 0; // # edges
      DerefScope scope;
      while (inputVar != -2) {
        switch (phaseI) {
          case 0: // constants
            switch (phaseJ) {
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
            std::cout << "Enter node #" << phaseJ << std::endl;
            break;
          case 2: // edge destinations, for each edge
            std::cout << "Enter edge #" << phaseJ << std::endl;
            break;
          default: // it should never get here.
            std::cout << "ERROR" << std::endl;
        }
        std::cin >> inputVar; // input. replace with file reading on port.
        switch (phaseI) {
          case 0: // constants
            switch (phaseJ) {
              case 0: // graph header token
                if (inputVar != 0xDEADBEEF) { // match the magic value...
                  std::cout << "ERROR" << std::endl; // or else you'll get an error
                }
                phaseJ++; // there's another subphase...
                break;
              case 1: // number of nodes
                inputNodes = inputVar;          // we've inputted the number of nodes.
                for (int k = 0; k < 2; k++) {
                  graph.push_back(scope, -1);   // expand the vector to accomodate the new constants
                }
                graphSet(&graph, &scope, 1, 0, inputNodes); // numNodes is at (1,0), set it to inputNodes
                phaseJ++; // there's another subphase...
                break;
              case 2: // number of edges
                inputEdges = inputVar;          // we've inputted the number of edges.
                graphSet(&graph, &scope, 0, 0, inputEdges); // numEdges is at (0,0), set it to inputEdges
                phaseI++;   // move to phase 1, 
                phaseJ = 0;         // subphase 0.
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
            graphSet(&graph, &scope, 2, phaseJ, inputVar); // put the node offset into the vector
            phaseJ++; // there's another subphase...
            if (phaseJ >= inputNodes) { // unless there isn't...
              phaseI = 2;  // then, move to phase 2...
              phaseJ = 0;                // subphase 0.
            }
            break;
          case 2: // edge offsets, for each edge
            for (int k = 0; k < 2; k++) {
              graph.push_back(scope, -1);   // expand the vector to accomodate the new edge
            }
            graphSet(&graph, &scope, 3, phaseJ, inputVar); // put the edge destination into the vector.
            phaseJ++; // there's another subphase...
            if (phaseJ >= inputEdges) { // unless there isn't...
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
    */
    for (int repeats = 0; repeats < 1; repeats++) { // for testing purposes. set to repeats<1 for normal use.
      {
        DerefScope scope;

        //graph.clear();          // clear the vectors from the previous iteration
        distances.clear();
        distances2.clear();
        distances3.clear();
        usedNumbers.clear();
        int arraySize = 30; // establish the size of the array
        std::random_device rd;  // set up the RNG
        std::mt19937_64 eng(rd());
        std::uniform_int_distribution<uint64_t> distr(0, arraySize - 1); // an edge has only a small chance to be present.

        { // - 
          //int edges[arraySize][arraySize];    // each entry indicates whether an edge exists between two nodes.
          //for (int i = 0; i < arraySize; i++) {
          //  edges[i][i] = 0;
          //  for (int j = 0; j < arraySize; j++) { // for each distinct pair of nodes:
          //    int randomNum = distr(eng);
          //    if (randomNum != 1 || j == i) {
          //      randomNum = 0;           // if the random num wasn't 1, set it to 0.
          //    }
          //    if (j == i + 1) {
          //      randomNum = 1;           // if the two nodes are consecutive, add the edge anyway.
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
              int randomNum = distr(eng);                     // select a random node
              bool isUnique = true;
              for (int j = 0; j < currentOutDegree; j++) {
                if (usedNumbers.at(scope, j) == randomNum) {
                  isUnique = false;                           // try again if the edge there already exists
                }
              }
              if (isUnique && randomNum != i) {             // if the edge isn't redundant...
                usedNumbers.push_back(scope, randomNum);      // add it
                currentOutDegree++;                           // count it
              }
            }
            for (int j = 0; j < outDegree; j++) {   // for each edge from this node...
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

        numEdges = graphAt(&graph, &scope, 0, 0);  // numEdges at 0,0
        numNodes = graphAt(&graph, &scope, 1, 0);  // numNodes at 1,0

        { // for incoming edges
          auto nodeCounts = manager->allocate_dataframe_vector<int>();   // for the indegree of each node
          auto nodeScatter = manager->allocate_dataframe_vector<int>();  // to adjust the offset of each node
          for (int i = 0; i < numNodes; i++) {   // at each node...
            nodeScatter.push_back(scope, 0);       // prefill the scatter and counts to 0
            nodeCounts.push_back(scope, 0);
          }
          for (int i = 0; i < numNodes; i++) {   // for each node...
            int startEdge = graphAt(&graph, &scope, 2, i); // establish the edge range: start
            int endEdge = (i == numNodes - 1) ? numEdges : graphAt(&graph, &scope, 2, i + 1);
            for (int j = startEdge; j < endEdge; j++) {
              int targetNode;
              {
                targetNode = graphAt(&graph, &scope, 3, j);
                // establish the neighboring nodes
              }
              nodeCounts.at_mut(scope, targetNode)++; // update the indegree of the target node
            }
          }
          graphSet(&graph, &scope, 4, 0, 0);    // the zeroth node has offset zero
          for (int i = 1; i < numNodes; i++) { // for each other node...
            graphSet(&graph, &scope, 4, i, graphAt(&graph, &scope, 4, i - 1) + nodeCounts.at_mut(scope, i - 1));
            // graph[4].push_back(scope, graph[4].at(scope, i - 1) + nodeCounts.at(scope, i - 1));
            // the offset is the sum of the previous offset and the indegree of the previous node
          }
          for (int i = 0; i < numNodes; i++) { // for each node...
            int startEdge = graphAt(&graph, &scope, 2, i); // establish the edge range: start
            int endEdge = (i == numNodes - 1) ? numEdges : graphAt(&graph, &scope, 2, i + 1);
            for (int j = startEdge; j < endEdge; j++) {
              int targetNode;
              {
                targetNode = graphAt(&graph, &scope, 3, j);
                // establish the neighboring nodes
              }
              int desiredElement = graphAt(&graph, &scope, 4, targetNode);   // start with the offset for the target node
              desiredElement += nodeScatter.at_mut(scope, targetNode);      // and add the scatter for the target node
              graphSet(&graph, &scope, 5, desiredElement, i);                 // use that as an index to add the edge
              nodeScatter.at_mut(scope, targetNode)++;                      // and increment the scatter, to not overwrite it
            }
          }
        }

        { // graph print 
          if (DO_PRINT) {
            for (int i = 0; i < numNodes; i++) {
              int starting = graphAt(&graph, &scope, 2, i);
              int ending = (i == numNodes - 1) ? (numEdges) : (graphAt(&graph, &scope, 2, i + 1));  // TODO: use this syntax for other ending nodes
              // int ending = (i == numNodes - 1) ? (numEdges) : (graphAt(graph, scope, 2, i + 1));
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



          //for (int i = 0; i < numNodes; i++) {
          //  std::cout << distances.at_mut(scope, i);
          //}
          //std::cout << std::endl;
        }
      }
        top_down_sort(manager, &distances, &graph);   // call each BFS algorithm
        bottom_up_sort(manager, &distances2, &graph);
        hybrid_sort(manager, &distances3, &graph);
      { // compare 
        DerefScope scope;
        if (DO_PRINT) {
          for (int i = 0; i < numNodes; i++) {
            std::cout << distances.at_mut(scope, i);  // print debug info: topdown distances
          }
          std::cout << std::endl;
          for (int i = 0; i < numNodes; i++) {
            std::cout << distances2.at_mut(scope, i); // print debug info: bottomup distances
          }
          std::cout << std::endl;
          for (int i = 0; i < numNodes; i++) {
            std::cout << distances3.at_mut(scope, i); // print debug info: bottomup distances
          }
          std::cout << std::endl;
        }
        bool isOkay = true;
        for (int i = 0; i < numNodes; i++) {
          if (distances.at_mut(scope, i) != distances2.at_mut(scope, i) || distances2.at_mut(scope, i) != distances3.at_mut(scope, i)) {
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
