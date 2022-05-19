#ifndef __BFS_H__
#define __BFS_H__

struct graph;

struct solution
{
    int* distances; // make into proper data type
};

struct vertex_set {
    int  count;
    int  alloc_count;
    int* present;   // same here
};


void bfs_bottom_up(graph* graph, solution* sol);
void bfs_top_down(graph* graph, solution* sol);
void bfs_hybrid(graph* graph, solution* sol);

#endif
