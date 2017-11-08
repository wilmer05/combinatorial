#include"edmonds.hpp"
#include <utility>
#include <vector>

namespace ALG{

    Edmonds :: Edmonds(size_type num_nodes) : graph(num_nodes), dsu(num_nodes){

        num_original_nodes = num_nodes;
        for(size_type i =static_cast<size_type>(0) ; i < num_nodes;i++)
            matched_to.push_back(i);
        
    }


    void Edmonds::reset_vectors(){
        while(pending_edges.size())
            pending_edges.pop();
        on_tree = std::vector<bool>(graph.num_nodes(), false);
        even_node = std::vector<bool>(graph.num_nodes(), false);
        odd_node = std::vector<bool>(graph.num_nodes(), false);
        parent = std::vector<size_type>(graph.num_nodes());
        
        for(size_type i =static_cast<size_type>(0) ; i < graph.num_nodes();i++)
            parent[i] = i;
    }

    bool Edmonds::exposed_vertex(size_type node){
        return matched_to[node] == node;
    }


    void Edmonds::extend_tree(size_type node, size_type even_parent){
        parent[node] = even_parent;
        size_type match = matched_to[node];
        
        odd_node[node] = true;
        parent[match] = node;
        even_node[match] = true;

        add_neighbors_to_pending_list( match );

    }

    void Edmonds::shrink(size_type node_u, size_type node_v){
    
    }

    void Edmonds::grow_tree(size_type root_node){
        on_tree[root_node] = true;
        add_neighbors_to_pending_list(root_node);
        even_node[root_node] = true;


        while(pending_edges.size()){
           std::pair<size_type, size_type> edge = pending_edges.front();
           pending_edges.pop();

           size_type node_x = edge.first;
           size_type node_y = edge.second;
          
           size_type repr_x = dsu.find(node_x);
           size_type repr_y = dsu.find(node_y);

           if(!even_node[repr_x] || odd_node[repr_y]){
                std::swap(repr_x, repr_y);
                if( !even_node[repr_x] || odd_node[repr_y]) 
                    continue;
           }

           //Condition 1: augment matching
           if(!on_tree[repr_y] && exposed_vertex(repr_y)){
                

           }
           //Extend tree
           else if(!on_tree[repr_y] && !exposed_vertex(repr_y)){
                extend_tree(repr_y, repr_x);     
           }
           //Shrink
           else{
                shrink(repr_x, repr_y);         
           }

        } 

        

    }
    //Edmonds Algrotihm implementation
    void Edmonds::run(){
       reset_vectors(); 
      
       for(size_type i =static_cast<size_type>(0) ; i < num_original_nodes; i++)
            if(matched_to[i] == i)
                grow_tree(i);

    }

    
    inline void  Edmonds::match(size_type node_u, size_type node_v){
          matched_to[node_u] = node_v;
          matched_to[node_v] = node_u;
    }

    inline void Edmonds::add_edge_to_pending_list(size_type node_u, size_type node_v){
        pending_edges.push(std::make_pair(node_u, node_v));
    }

    inline void Edmonds::add_neighbors_to_pending_list(size_type node){
        for(auto const &neighbour : graph.node( node ).neighbors()){
            size_type representative = dsu.find(neighbour);
            if(!on_tree[representative]){
                add_edge_to_pending_list(node, neighbour); 
                on_tree[representative] = true;
            }
        }
    }
};
