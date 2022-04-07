#include <stdexcept>
#include "mutation_annotated_dag.hpp"

namespace MAT = Mutation_Annotated_Tree;
namespace MAD = Mutation_Annotated_DAG;

int main(int argc, char** argv) {
    MAD::DAG history = MAD::load_mat_protobuf(argv[0]);
}

MAD::DAG MAD::mat_to_dag(MAT::Tree & tree) {
    // I understand that although I will be returning the DAG by value, since
    // it contains pointers to objects initialized in this scope, I must use
    // new for each such object if I want it to not be deleted on return?
    std::vector<Mutation_Annotated_Tree::Node*> tree_dfs = tree.depth_first_expansion();
    const size_t n_nodes = tree_dfs.size();
    MAD::Node* dag_dfs[n_nodes];
    // build dag nodes from tree nodes, starting with leaves
    for (auto node : tree_dfs){
        std::vector<MAD::EdgeVector *>* clades = new std::vector<MAD::EdgeVector *>;
        MAD::Node * p_newnode = new MAD::Node;
        // build child edges, one per clade
        for (auto child : node->children){
            // DFS implies dag_dfs will contain an already constructed node at
            // this index
            MAD::Edge * child_edge = new MAD::Edge(p_newnode, dag_dfs[child->dfs_idx], child->mutations);
            MAD::EdgeVector * clade = new MAD::EdgeVector;
            clade->push_back(child_edge);
            clades->push_back(clade);
        }
        p_newnode->clades = *clades;
        dag_dfs[node->dfs_idx] = p_newnode;
    }
    // TODO need to add a UA node, and compute the right mutations from
    // some passed root sequence.
    MAD::DAG newdag(dag_dfs[0]);
    return newdag;
}

MAD::DAG MAD::load_mat_protobuf(const std::string& protobuf_filename){
    MAT::Tree tree = MAT::load_mutation_annotated_tree(protobuf_filename);
    return MAD::mat_to_dag(tree);
}

MAD::Node::Node(std::vector<EdgeVector *> clades){
    for (auto edgevector : clades){
        if (edgevector->empty()){
            throw std::invalid_argument("Each clade must have at least one descendant edge.");
        }
    }
    this->clades = clades;
}

bool MAD::Node::is_leaf(){
    return clades.empty();
}

bool MAD::Node::is_ua_node(){
    return rootward_edges.empty();
}

void MAD::Node::add_child_edge(size_t clade_idx, MAD::Edge* edge){
    clades[clade_idx]->push_back(edge);
}

void MAD::Node::add_parent_edge(MAD::Edge* edge){
    rootward_edges.push_back(edge);
}

void MAD::Node::remove_parent_edge(MAD::Edge* edge){
    /* Remove a parent edge by pointer */
    // Slow and ugly, perhaps there's a better way?
    // e.g. using some sort of set would enforce there isn't another pointer
    // to an identical edge in the list...
    rootward_edges.erase(std::remove(rootward_edges.begin(), rootward_edges.end(), edge), rootward_edges.end());
}

std::set<MAT::Mutation> MAD::Node::get_muts_rel_reference(){
    std::set<MAT::Mutation> muts;
    MAD::Node *current_node = this;
    while (not current_node->is_ua_node()){
        /*current_node must then have a parent edge*/
        MAD::Edge *parent_edge = current_node->rootward_edges[0];
        for (auto mut : parent_edge->mutations){
            muts.insert(mut);
        }
        current_node = parent_edge->parent;
    }
    return muts;
}

std::string get_seq(){
    throw std::runtime_error("Not Implemented");
    return std::string();
}

MAD::Edge::Edge(MAD::Node * parent, MAD::Node * child, std::vector<MAT::Mutation> mutations){
    /* All edge methods will manage addition of parents, but require children to be explicitly specified. */
    this->parent = parent;
    this->child = child;
    child->add_parent_edge(this);
}

MAD::Edge::~Edge(){
    child->remove_parent_edge(this);
    // and then... member objects are destroyed automatically?
}

void MAD::Edge::add_mutation(MAT::Mutation mut){
    mutations.push_back(mut);
}

void MAD::Edge::clear_mutations(){
    mutations.clear();
}

MAD::DAG::DAG(Node *root){
    this->root = root;
}

MAD::DAG::DAG(Node *root, std::string reference_sequence){
    this->root = root;
    this->reference_sequence = reference_sequence;
}
