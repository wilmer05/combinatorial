#ifndef ALG_HPP
#define ALG_HPP

#include "graph.hpp"
#include "dsu.hpp" 
#include <queue>

/**
    @file edmonds.hpp

    @brief This file implements the Edmond's algorithm
**/
namespace ALG{ // for the implementation of Edmonds Algorithm

    using size_type = std::size_t;


    /** 
        @class Edmonds

        @brief Implements the edmonds algorithm, it has every operation in different methods, like shrinking, augment and extension
    **/
    class Edmonds{
    
        public:
           
           /** @brief Constructor, **/
           Edmonds(size_type);


           void reset_vectors();

           void grow_tree(size_type);

           void run();
           
//            ED::Graph extract_solution();

           void  match(size_type, size_type);

           void add_edge_to_pending_list(size_type, size_type);

           void add_neighbors_to_pending_list(size_type);


           bool exposed_vertex(size_type);

           void extend_tree(size_type, size_type);

           void shrink(size_type, size_type);

           void remove_not_contained_edges(std::vector<std::pair<size_type, size_type> > &edges);

            void augment_matching(size_type node, size_type parent);

        private:
            ED::Graph graph;
            DSU::Dsu dsu;
            std::vector<size_type> matched_to; 
            std::queue<std::pair<size_type, size_type> > pending_edges; 
            size_type num_original_nodes;
            std::vector<bool> on_tree;
            std::vector<bool> even_node;
            std::vector<bool> odd_node;
            std::vector<size_type> parent;


            //Do I need this? \/
            std::vector<std::pair<int,int> > cycle_edge;
            std::vector<bool> visited;
            std::vector<std::vector<std::pair< size_type, size_type> > > edges_in_cycle; 
            std::vector<size_type> node_closes_cycle;
            std::vector<size_type> nodes_in_tree;
            std::vector<size_type> frustrated;
    }; 


};


#endif
