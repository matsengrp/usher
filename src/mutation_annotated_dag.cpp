#include "mutation_annotated_dag.hpp"

namespace MAT = Mutation_Annotated_Tree;
namespace MAD = Mutation_Annotated_DAG;

MAD::DAG* MAD::mat_to_dag(MAT::Tree* tree) {
    
}

MAD::Node::Node(std::vector<EdgeVector *> clades, ){
    for (auto edgevector : clades){
        if (edgevector->empty()){
            throw std::invalid_argument("Each clade must have at least one descendant edge.")
        }
    }
    this->clades = clades;
}

void MAD::Node::is_leaf(){
    return clades.empty()
}

void MAD::Node::add_parent_edge(MAD::Edge* edge){
    rootward_edges.push_back(edge);
}

MAD::Edge::Edge(MAD::Node * parent, MAD::Node * child, std::vector<MAT::Mutation> mutations){
    /* All edge methods will manage addition of parents, but require children to be explicitly specified. */
    this->parent = parent;
    this->child = child;
    /* What happens when this object is destroyed? */
    child->add_parent_edge(this);
}

void MAD::Edge::add_mutation(MAT::Mutation mut){
    mutations.push_back(mut);
}

void MAD::Edge::clear_mutations(){
    mutations.clear()
}
