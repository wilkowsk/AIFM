#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "CycleTimer.h"
#include "bfs.h"
#include "graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define BOTTOMUP_NOT_VISITED_MARKER    0




void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) { // replace with something more suitable for AIFM arrays
    list->alloc_count = count;
    list->present = (int*)malloc(sizeof(int) * list->alloc_count);
    vertex_set_clear(list);
}



void bottom_up_step(
    graph* g,
    vertex_set* frontier,    // auto frontier = manager->allocate_array<vertex_set, kNumEntries>();
    int* distances,  // auto distances = manager->allocate_array<vertex_set, kNumEntries>();
    int iteration /*, DerefScope scope*/) {
    for (int i = 0; i < g->num_nodes; i++) {
        if (frontier->present[i] == BOTTOMUP_NOT_VISITED_MARKER) {
            int start_edge = g->incoming_starts[i]; // int start_edge = g->incoming_starts.at(scope, i);
            int end_edge = (i == g->num_nodes-1)? g->num_edges : g->incoming_starts[i + 1]; // int end_edge = (i == g->num_nodes-1)? g->num_edges : g->incoming_starts.at(scope, i+1);
            for(int neighbor = start_edge; neighbor < end_edge; neighbor ++) {
                int incoming = g->incoming_edges[neighbor]; // int incoming = g->incoming_edges.at(scope, neighbor);
                if(frontier->present[incoming] == iteration) { // if(frontier->present.at(scope, incoming) == iteration) {
                    distances[i] = distances[incoming] + 1; // distances.at_mut(scope, i) = distances.at(scope, incoming) + 1;
                    frontier->count++;
                    frontier->present[i] = iteration + 1; // frontier->present.at_mut(scope, i) = iteration + 1;
                    break;
                }
            }
        }
    }

}

void bfs_bottom_up(graph* graph, solution* sol /*, DerefScope scope*/)
{
    // 15-418/618 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set list1;
    // vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    // vertex_set_init(&list2, graph->num_nodes);
    int iteration = 1;

    vertex_set* frontier = &list1; // auto frontier = manager->allocate_array<vertex_set, kNumEntries>();
    // vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    // for (int i=0; i<graph->num_nodes; i++)
        // sol->distances[i] = NOT_VISITED_MARKER;
    
        
    // setup frontier with the root node    
    // just like put the root into queue
    frontier->present[frontier->count++] = 1; // frontier->present.at_mut(scope, frontier->count++) = 1;

    // set the root distance with 0
    sol->distances[ROOT_NODE_ID] = 0; // sol->distances.at_mut(scope, ROOT_NODE_ID) = 0;
    

    // just like pop the queue
    while (frontier->count != 0) {
        frontier->count = 0;

        double start_time = CycleTimer::currentSeconds();
        

        bottom_up_step(graph, frontier, sol->distances, iteration);

        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);

        iteration++;

    }



}




// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    graph* g,
    vertex_set* frontier, // make correct data types
    vertex_set* new_frontier,
    int* distances /*, DerefScope scope*/) // add DerefScope scope
{

    for (int i=0; i<frontier->count; i++) {

        int node = frontier->present[i]; // int node = frontier->present.at(scope, i);

        int start_edge = g->outgoing_starts[node]; // int start_edge = g->outgoing_starts.at(scope, node);
        int end_edge = (node == g->num_nodes-1) ? g->num_edges : g->outgoing_starts[node+1]; // int end_edge = (node == g->num_nodes-1) ? g->num_edges : g->outgoing_starts.at(scope, node+1);

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor]; // int outgoing = g->outgoing_edges.at(scope, neighbor);

            if (distances[outgoing] == NOT_VISITED_MARKER) { // if (distances.at(scope, outgoing) == NOT_VISITED_MARKER) {
                distances[outgoing] = distances[node] + 1; // distances.at_mut(scope, outgoing) = distances.at(scope, node) + 1;
                int index = new_frontier->count++;
                new_frontier->present[index] = outgoing; // new_frontier->present.at_mut(scope, index) = outgoing;
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(graph* graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1; // auto frontier = manager->allocate_array<vertex_set, kNumEntries>();
    vertex_set* new_frontier = &list2; // auto new_frontier = manager->allocate_array<vertex_set, kNumEntries>();

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER; // sol->distances.at_mut(scope, i) = NOT_VISITED_MARKER;
    
        
    // setup frontier with the root node    
    // just like put the root into queue
    frontier->present[frontier->count++] = ROOT_NODE_ID; // frontier->present.at_mut(scope, frontier->count++) = ROOT_NODE_ID;

    // set the root distance with 0
    sol->distances[ROOT_NODE_ID] = 0; // sol->distances.at_mut(scope, ROOT_NODE_ID) = 0;

    // just like pop the queue
    while (frontier->count != 0) {

        double start_time = CycleTimer::currentSeconds();

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);

        // swap pointers. good thing these are pointers or we'd have to do something more complex
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}
