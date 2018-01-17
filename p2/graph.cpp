#include "graph.hpp" // always include corresponding header first

#include <ostream>
#include<iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "dsu.hpp"

namespace ED
{
void Node::set_lambda(double l){
    lambda = l;
} 

/////////////////////////////////////////////
//! \c Node definitions
/////////////////////////////////////////////
void Node::set_coordinates(double new_x, double new_y) {
     x = new_x;
     y = new_y;
}

void Node::add_neighbor(NodeId const id)
{
   _neighbors.push_back(id);
}

/////////////////////////////////////////////
//! \c Graph definitions
/////////////////////////////////////////////

Edge & Graph :: get_edge(size_type idx){
    return _edges[idx];
}

Edge & Graph :: get_backup_edge(size_type u, size_type v){
    return _backup_edges[get_edge_index(u,v)];

}

Graph::Graph(NodeId const num_nodes)
   :
   _nodes(num_nodes),
   _num_edges(0)
{}

Graph::~Graph()
{
    _nodes.clear();
}

Edge::Edge (Node *u, Node *v, double d) : node_u(u), node_v(v), dist(d)
{}

void Graph :: generate_edges(){
    size_type num_nodes = _nodes.size();
    for(size_type i=0 ; i < num_nodes;  i++){
        _nodes[i].set_id(i);
        for(size_type j =0 ; j < i; j++){
            add_edge(j, i);
        }
    }
}

void Graph :: reset_edges() {
    _edges = _backup_edges;
} 

void Graph::fix_forbidden_edges(std::vector<std::pair<size_type, size_type > > &F, std::vector<std::vector<size_type> > &R){
     
    //We change the cost to infinity of every edge in the 
    //list F
    for(size_type i =0 ; i < F.size(); i++){
        _edges[get_edge_index(F[i].first, F[i].second)].dist = 1e14;
    }

    //If there is a node with 2 required edges
    //then we forbid every other edge
    for(size_type i =0 ; i < R.size(); i++){
        assert(R[i].size() <= 2);
        if(R[i].size() == 2){
            for(size_type j = 0; j < R.size(); j++)
                if(R[i][0] != j && R[i][1] != j && i != j){
                    _edges[get_edge_index(i, j)].dist = 1e14;
                }
                    
        }
    }
}

void Graph::fix_lambdas_and_sort_edges(std::vector<double> &lambdas){
    
    for(size_type i =0 ; i < lambdas.size(); i++)
        _nodes[i].set_lambda(lambdas[i]);
    
    sort(_edges.begin(), _edges.end()); 
    /*std::cout << "##############\n";
    for(size_type i =0 ; i < _edges.size(); i++)
        std::cout << _edges[i].dist << " ids=" <<_edges[i].get_first()->get_id() << ", " << _edges[i].get_second()->get_id()<<" <-> " <<" lambdas" << _edges[i].get_first()->get_lambda() << " " << _edges[i].get_second()->get_lambda() << "\n";*/
}

void Graph::add_edge(NodeId node1_id, NodeId node2_id)
{
   if (node1_id == node2_id)
   {
      throw std::runtime_error("ED::Graph class does not support loops!");
   }

   // minimum redundancy :-), maybe a bit overkill...
   auto impl = [this](NodeId a, NodeId b)
   {
      Node & node = _nodes.at(a);
      node.add_neighbor(b);
   };

   impl(node1_id, node2_id);
   impl(node2_id, node1_id);

   _num_edges+=2;
   _backup_edges.push_back(ED :: Edge(&_nodes[node1_id], &_nodes[node2_id], get_distance(node1_id, node2_id)));
}

void Graph::set_coordinates(size_type node, double new_x, double new_y){
    _nodes[node].set_coordinates(new_x, new_y);
}

int Graph::distance(double x1, double y1, double x2, double y2)
{
    return std::lround(std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
}

int Graph::get_distance(size_type node_u, size_type node_v){

    return distance(_nodes[node_u].x, _nodes[node_u].y, _nodes[node_v].x, _nodes[node_v].y);

}

std::ostream & operator<<(std::ostream & str, Graph const & graph)
{
   str << "c This encodes a graph in DIMACS format\n"
       << "p edge " << graph.num_nodes() << " " << graph.num_edges() << "\n";

   for (NodeId node_id = 0; node_id < graph.num_nodes(); ++node_id)
   {
      auto const & node = graph.node(node_id);

      for (auto const & neighbor_id : node.neighbors())
      {
         // output each edge only once
         if (node_id < neighbor_id)
         {
            str << "e " << to_dimacs_id(node_id) << " " << to_dimacs_id(neighbor_id) << "\n";
         }
      }
   }

   str << std::flush;
   return str;
}


/////////////////////////////////////////////
//! global functions
/////////////////////////////////////////////

NodeId from_dimacs_id(DimacsId const dimacs_id)
{
   if (dimacs_id == 0)
   {
      throw std::runtime_error("Invalid (0) DIMACS id.");
   }

   return static_cast<NodeId>(dimacs_id - 1);
}

DimacsId to_dimacs_id(NodeId const node_id)
{
   if (node_id == std::numeric_limits<NodeId>::max())
   {
      throw std::runtime_error("Invalid (inf) node id.");
   }

   return static_cast<DimacsId>(node_id + 1);
}

} // namespace ED
