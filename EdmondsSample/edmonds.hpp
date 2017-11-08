#ifndef ALG_HPP
#define ALG_HPP

#include "graph.hpp"
#include "dsu.hpp" 
#include <queue>
namespace ALG{

    using size_type = std::size_t;

    class Edmonds{
    
        public:
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

    }; 


};


#endif
