#pragma once
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <cassert>
#include "usher_common.hpp"
/* #include <tbb/flow_graph.h> */
/* #include <tbb/reader_writer_lock.h> */
/* #include <tbb/scalable_allocator.h> */
/* #include <tbb/task_scheduler_init.h> */
/* #include <tbb/blocked_range.h> */
/* #include <tbb/task_group.h> */
/* #include <tbb/tbb.h> */
/* #include <tbb/mutex.h> */
/* #include "parsimony.pb.h" */
/* #include "Instrumentor.h" */

/* #if SAVE_PROFILE == 1 */
/* #  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__); */
/* #else */
/* #  define TIMEIT() */
/* #endif */

namespace MAT = Mutation_Annotated_Tree;

namespace Mutation_Annotated_DAG {

Mutation_Annotated_DAG::DAG* mat_to_dag(MAT::Tree* tree)

class Edge;
using EdgeVector = std::vector<Edge *>;

class Node {
  public:
    size_t id;
    EdgeVector rootward_edges;
    std::vector<EdgeVector *> clades;

    bool is_leaf();
    bool is_ua_node();
    std::vector<MAD::Node> children();

    bool add_edge(size_t clade_idx, MAD::Edge* edge);
    void add_parent_edge(MAD::Edge* edge);

    Node(std::vector<EdgeVector *> clades, );

};

class Edge {
  public:
    Node *parent;
    Node *child;
    std::vector<Mutation_Annotated_Tree::Mutation> mutations;

    void add_mutation(Mutation_Annotated_Tree::Mutation mut);
    void clear_mutations();

    Edge(MAD::Node * parent, MAD::Node * child, std::vector<MAT::Mutation> mutations);
};

class DAG {
  private:
    std::unordered_map<size_t, Node *> all_nodes;

  public:
    Tree() {
        root = NULL;
        curr_internal_node = 0;
        all_nodes.clear();
    }
    Node *root;

    void merge(Mutation_Annotated_DAG::DAG* newdag)
    Mutation_Annotated_DAG::DAG* sample()
    // tbb::concurrent_unordered_map<std::string, std::vector<std::string>>
    // condensed_nodes; tbb::concurrent_unordered_set<std::string>
    // condensed_leaves;
};
} // namespace Mutation_Annotated_DAG

