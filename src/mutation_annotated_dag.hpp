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

    Node();
    Node(std::vector<EdgeVector *> clades);

    bool is_leaf() const;
    bool is_ua_node() const;
    LazyIterator<Node *> children() const;

    void add_child_edge(Edge*);
    void add_child_edge(size_t clade_idx, Edge*);
    void add_parent_edge(Edge* edge);
    void remove_parent_edge(Edge* edge);
    std::set<MAT::Mutation> get_muts_rel_reference();
    std::string get_sequence() const;


};

class Edge {
  private:
    AbstractWeight dynamic_programming_data;

  public:
    Node *parent;
    Node *child;
    std::vector<Mutation_Annotated_Tree::Mutation> mutations;
    double probability;

    void add_mutation(Mutation_Annotated_Tree::Mutation mut);
    void clear_mutations();

    Edge();
    Edge(Node * parent, Node * child, std::vector<MAT::Mutation> mutations);
    ~Edge();
};

class DAG {
  private:
    void postorder_helper(Node*, std::unordered_set<Node *> &, std::vector<Node *> &);


  public:
    //made this public to allow setting in mat_to_dag. I'm sure there's a
    //better way.
    std::unordered_map<size_t, Node *> all_nodes;
    std::string reference_sequence;
    Node *root;

    void merge(DAG* newdag);
    std::vector<Node *> postorder();
    DAG* sample() const;

    DAG(Node *);
    DAG(Node *, std::string);

    ArbitraryPrecisionInteger count_histories() const;
    /* can be done as a call to postorder_history_weight_accumulation */

    void write_protobuf(std::string filename);
    /* ... in our own protobuf format like the draft that Russ suggested we use
     */

    DAG sample_history() const;
    /* Sample a history from the DAG, choosing edges below each clade according
     * to relative Edge.probability values
     */

    DAG get_history_by_index(ArbitraryPrecisionInteger) const;
    /* Use the indexing scheme developed by @clarisw to map integer indices
     * uniquely to histories in the DAG.
     */

    bool is_clade_tree() const;
    /* Return whether this history DAG is a history, i.e. contains exactly one
     * descendant edge from each node-clade pair. Equivalent to
     * count_histories() == 1 */

    void merge(std::vector<DAG*>);
    /* Merge the provided list of history DAGs into this one. Passing a list
     * allows reuse of the node mapping, instead of having to re-build it for
     * each merged DAG
     */

    void add_all_allowed_edges();
    /* Add all possible valid edges to the DAG. An edge is valid iff the target
     * node's clade union (union of child clades) is one of the child clades of
     * the parent node.
     */


    // ###### Abstract Dynamic Programming Methods #######
    // This is how I laid this out in Python, very flexible but perhaps not
    // ideal in c++. Here, if the passed functions require full sequences at
    // nodes, that may be a problem...
    // An AbstractWeight could be an int (like for parsimony), a double
    // (like for probabilities), or a container class like a (Counter) map.
    void postorder_history_weight_accumulation(
            AbstractWeight (*leaf_func)(Node&),
            AbstractWeight (*edge_weight_func)(Edge&),
            AbstractWeight (*accum_within_clade)(vector<AbstractWeight>),
            AbstractWeight (*accum_between_clade)(vector<AbstractWeight>)
    /* This simply implements the (I find, very common) idiom by which leaf nodes are assigned
     * a weight (stored in Node.dynamic_programming_data) via leaf_func, and
     * all other nodes are assigned a weight which is the result of
     * accum_between_clade on the weights associated to clades.
     * Each clade's weight is the result of accum_within_clade on weights
     * associated with all edges descending from that clade. An edge's weight
     * is calculated by calling edge_weight_func on the edge's parent and
     * target nodes, and adding the result to the child node's weight stored in
     * Node.dynamic_programming_data.
     *
     * The edge_weight_func could be non-symmetric, so it should always be
     * called with an edge's parent node as the first argument.
     *
     * the accum functions could take two AbstractWeights instead of a list,
     * but this way we can define what happens when accumulating an empty list
     * of weights (what's the additive/multiplicative identity of our weights)
     */

    CounterMap weight_count(
            Weight (*leaf_func)(Node&),
            Weight (*edge_weight_func)(Edge&),
            Weight (*accum_func)(vector<Weight>)
    );
    /* Counts the Weights of all histories in the hDAG, such as their parsimony
     * scores, and returns a map with observed weights as keys, and
     * counts of histories with each weight as values.
     *
     * This is just a call to postorder_history_weight_accumulation, where the
     * AbstractWeight object being accumulated is a CounterMap, where the keys are
     * Weights (in the sense of arguments to this method), and values are
     * counts. In the call to postorder_history_weight_accumulation, leaf_func
     * and edge_weight_func wrap the value from the leaf_func and edge_weight_func
     * passed to this function in a CounterMap with count 1. Accum_within_clade
     * is counter_sum, and accum_between_clades is counter_prod, using the
     * passed accum_func.
     */

    Weight optimal_weight_annotate(
            Weight (*leaf_func)(Node&),
            Weight (*edge_weight_func)(Edge&),
            Weight (*accum_func)(vector<Weight>),
            Weight (*optimal_func)(vector<Weight>)
    );
    /* Returns the optimal weight of a tree in the history DAG, annotating each
     * Node with the optimal weight of a sub-history beneath it (in
     * Node.dynamic_programming_data)
     *
     * This is just a call to postorder_history_weight_accumulation in which
     * accum_within_clade is optimal_func, and accum_between_clade is
     * accum_func.
     */


    void trim_optimal_weight(
            Weight (*leaf_func)(Node&),
            Weight (*edge_weight_func)(Edge&),
            Weight (*accum_func)(vector<Weight>),
            Weight (*optimal_func)(vector<Weight>)
    );
    /* First calls optimal_weight_annotate, then uses annotations to
     * trim the dag to only express optimal weight trees, by choosing the
     * optimal weight edges (comparing edge weight plus minimum weight of
     * a sub-history below that edge's target node).
     *
     * Or this could instead be done in a single postorder traversal, where
     * trimming happens at the same time as annotation.
     */
}
    
class CounterMap{
    /* Just something that counts items, like `collections.Counter` in Python.
     */
}

CounterMap counter_product(
        std::vector<CounterMap> counterlist,
        Weight (*accum_func)(vector<Weight>)
        )
/* A sort of cartesian product, which counts all the ways a Weight can
 * be achieved by applying accum_func to combinations of one key from each
 * passed counter in counterlist.
 */

CounterMap counter_sum(
        std::vector<CounterMap> counterlist
        )
/* Equivalent to unioning Counters: The returned counter has all the keys
 * from the passed list of counters, and when a key appears in multiple
 * counters, the resulting count is the sum of all the counts for that key
 * in counterlist.
 */


DAG mat_to_dag(MAT::Tree&);
DAG load_mat_protobuf(const std::string&);

DAG read_protobuf(const std::string&);
/* in our own protobuf format */
} // namespace Mutation_Annotated_DAG

