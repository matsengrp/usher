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
#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/task_group.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include "parsimony.pb.h"
#include "Instrumentor.h"

/* #if SAVE_PROFILE == 1 */
/* #  define TIMEIT() InstrumentationTimer timer##__LINE__(__PRETTY_FUNCTION__); */
/* #else */
/* #  define TIMEIT() */
/* #endif */

namespace MAT = Mutation_Annotated_Tree;

namespace Mutation_Annotated_DAG {


class Edge;
using EdgeVector = std::vector<Edge *>;

class Node {
  public:
    size_t id;
    EdgeVector rootward_edges;
    std::vector<EdgeVector *> clades;

    bool is_leaf();
    bool is_ua_node();
    std::vector<Node> children();

    void add_child_edge(Edge*);
    void add_child_edge(size_t clade_idx, Edge*);
    void add_parent_edge(Edge* edge);
    void remove_parent_edge(Edge* edge);
    std::set<MAT::Mutation> get_muts_rel_reference();
    std::string get_seq();

    Node();
    Node(std::vector<EdgeVector *> clades);
};

class Edge {
  public:
    Node *parent;
    Node *child;
    std::vector<Mutation_Annotated_Tree::Mutation> mutations;

    void add_mutation(Mutation_Annotated_Tree::Mutation mut);
    void clear_mutations();

    Edge();
    Edge(Node * parent, Node * child, std::vector<MAT::Mutation> mutations);
    ~Edge();
};

class DAG {
  public:
    //made this public to allow setting in mat_to_dag. I'm sure there's a
    //better way.
    std::unordered_map<size_t, Node *> all_nodes;
    std::string reference_sequence;
    Node *root;

    void merge(DAG* newdag);
    std::vector<Node *> postorder();
    DAG* sample();

    DAG(Node *);
    DAG(Node *, std::string);
};

DAG mat_to_dag(MAT::Tree&);
DAG load_mat_protobuf(const std::string&);
} // namespace Mutation_Annotated_DAG

