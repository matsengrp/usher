#include <stdexcept>
#include <iostream>
#include "mutation_annotated_dag.hpp"

using namespace Mutation_Annotated_DAG;

int main(int argc, char** argv) {
    std::cout << "loading protobuf" << argv[1] << std::endl;
    DAG history = load_mat_protobuf(argv[1]);
    std::cout << "loaded protobuf" << argv[1] << std::endl;
    std::vector<Node*> post_order = history.postorder();
    std::cout << "DAG has " << post_order.size() << " nodes" << std::endl;
    std::cout << "id of ua node is " << history.root->id << "" << std::endl;
    return 0;
}

DAG Mutation_Annotated_DAG::mat_to_dag(MAT::Tree & tree) {
    // I understand that although I will be returning the DAG by value, since
    // it contains pointers to objects initialized in this scope, I must use
    // new for each such object if I want it to not be deleted on return?
    auto dag_node_indices = new std::unordered_map<size_t, Node *>;
    std::vector<Mutation_Annotated_Tree::Node*> tree_dfs = tree.depth_first_expansion();
    std::cout << "tree has " << tree_dfs.size() << " nodes" << std::endl;
    const size_t n_nodes = tree_dfs.size();
    Node* dag_dfs[n_nodes];
    // build dag nodes from tree nodes, starting with leaves
    for (auto node = tree_dfs.rbegin(); node != tree_dfs.rend(); node++){
        std::vector<EdgeVector *>* clades = new std::vector<EdgeVector *>;
        Node * p_newnode = new Node;
        // build child edges, one per clade
        for (Mutation_Annotated_Tree::Node * child : (*node)->children){
            // DFS implies dag_dfs will contain an already constructed node at
            // this index
            Edge * child_edge = new Edge(p_newnode, dag_dfs[child->dfs_idx], child->mutations);
            EdgeVector * clade = new EdgeVector;
            clade->push_back(child_edge);
            clades->push_back(clade);
        }
        p_newnode->clades = *clades;
        size_t old_idx = (*node)->dfs_idx;
        dag_dfs[old_idx] = p_newnode;
        p_newnode->id = old_idx + 1;
        (*dag_node_indices)[old_idx + 1] = p_newnode;
    }
    Node * ua_node = new Node;
    ua_node->id = 0;
    (*dag_node_indices)[0] = ua_node;
    std::vector<MAT::Mutation> mutations;
    // TODO compute correct mutations rel ua node sequence
    Edge * ua_child_edge = new Edge(nullptr, dag_dfs[0], mutations);
    ua_node->add_child_edge(ua_child_edge);
    DAG newdag(ua_node);
    newdag.all_nodes = *dag_node_indices;
    return newdag;
}

std::vector<Node *> Mutation_Annotated_DAG::DAG::postorder(){
    std::vector<Node *> visit_order;
    std::unordered_set<Node *> visited;
    DAG::postorder_helper(root, visited, visit_order);
    return visit_order;
}

void Mutation_Annotated_DAG::DAG::postorder_helper(Node* current_node, std::unordered_set<Node*> & visited, std::vector<Node *> & visit_order){
    if (not visited.count(current_node)){
        for (auto child : current_node->children()){
            DAG::postorder_helper(child, visited, visit_order);
        }
        current_node->id = visit_order.size();
        visit_order.push_back(current_node);
    }
}

DAG Mutation_Annotated_DAG::load_mat_protobuf(const std::string& protobuf_filename){
    MAT::Tree tree = MAT::load_mutation_annotated_tree(protobuf_filename);
    return mat_to_dag(tree);
}

Mutation_Annotated_DAG::Node::Node(){
    id = 0;
}

Mutation_Annotated_DAG::Node::Node(std::vector<EdgeVector *> clades){
    for (auto edgevector : clades){
        if (edgevector->empty()){
            throw std::invalid_argument("Each clade must have at least one descendant edge.");
        }
    }
    this->clades = clades;
}

bool Mutation_Annotated_DAG::Node::is_leaf(){
    return clades.empty();
}

bool Mutation_Annotated_DAG::Node::is_ua_node(){
    return rootward_edges.empty();
}

std::vector<Node *> Mutation_Annotated_DAG::Node::children(){
    //could I return some sort of iterator instead of initializing a whole new
    //vector for this?
    std::vector<Node *> child_list;
    for (auto clade : clades){
        for (auto child_edge : *clade){
            child_list.push_back(child_edge->child);
        }
    }
    return child_list;
}


void Mutation_Annotated_DAG::Node::add_child_edge(Edge* edge){
    if (clades.size() > 1){
        throw std::invalid_argument("To add an edge to a node with more than one clade you must specify a clade index.");
    } else if (clades.empty()){
        clades.push_back(new EdgeVector);
    } // clades contains one EdgeVector now.
    clades[0]->push_back(edge);
}

void Mutation_Annotated_DAG::Node::add_child_edge(size_t clade_idx, Edge* edge){
    clades[clade_idx]->push_back(edge);
}

void Mutation_Annotated_DAG::Node::add_parent_edge(Edge* edge){
    rootward_edges.push_back(edge);
}

void Mutation_Annotated_DAG::Node::remove_parent_edge(Edge* edge){
    /* Remove a parent edge by pointer */
    // Slow and ugly, perhaps there's a better way?
    // e.g. using some sort of set would enforce there isn't another pointer
    // to an identical edge in the list...
    rootward_edges.erase(std::remove(rootward_edges.begin(), rootward_edges.end(), edge), rootward_edges.end());
}

std::set<MAT::Mutation> Mutation_Annotated_DAG::Node::get_muts_rel_reference(){
    std::set<MAT::Mutation> muts;
    Node *current_node = this;
    while (not current_node->is_ua_node()){
        /*current_node must then have a parent edge*/
        Edge *parent_edge = current_node->rootward_edges[0];
        for (auto mut : parent_edge->mutations){
            muts.insert(mut);
        }
        current_node = parent_edge->parent;
    }
    return muts;
}

std::string Mutation_Annotated_DAG::Node::get_seq(){
    throw std::runtime_error("Not Implemented");
    return std::string();
}

Mutation_Annotated_DAG::Edge::Edge(Node * parent, Node * child, std::vector<MAT::Mutation> mutations){
    /* All edge methods will manage addition of parents, but require children to be explicitly specified. */
    this->parent = parent;
    this->child = child;
    child->add_parent_edge(this);
}

Mutation_Annotated_DAG::Edge::~Edge(){
    child->remove_parent_edge(this);
    // and then... member objects are destroyed automatically?
}

void Mutation_Annotated_DAG::Edge::add_mutation(MAT::Mutation mut){
    mutations.push_back(mut);
}

void Mutation_Annotated_DAG::Edge::clear_mutations(){
    mutations.clear();
}

Mutation_Annotated_DAG::DAG::DAG(Node *root){
    this->root = root;
}

Mutation_Annotated_DAG::DAG::DAG(Node *root, std::string reference_sequence){
    this->root = root;
    this->reference_sequence = reference_sequence;
}
