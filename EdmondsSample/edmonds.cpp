#include"edmonds.hpp"
#include <utility>
#include <vector>

namespace ALG{

    Edmonds :: Edmonds(size_type num_nodes) : graph(num_nodes), dsu(num_nodes), on_tree(num_nodes, false), visited(num_nodes, false), frustrated(num_nodes, false){

        num_original_nodes = num_nodes;
        for(size_type i =static_cast<size_type>(0) ; i < num_nodes;i++)
            matched_to.push_back(i);
         
        for(size_type i = static_cast<size_type>(0) ; i < graph.num_nodes();i++)
            parent[i] = i;

        nodes_in_tree.clear();
        edges_in_cycle.clear();
        
    }


    void Edmonds::reset_vectors(){
        while(pending_edges.size())
            pending_edges.pop();

        /*on_tree = std::vector<bool>(graph.num_nodes(), false);
        even_node = std::vector<bool>(graph.num_nodes(), false);
        odd_node = std::vector<bool>(graph.num_nodes(), false);
        */

        
        //Reset parent vectors
        for(size_type i = static_cast<size_type>(0); i < nodes_in_tree.size(); i++){
            parent[nodes_in_tree[i]] = nodes_in_tree[i];
            on_tree[nodes_in_tree[i]] = false;
            even_node[nodes_in_tree[i]] = false;
            odd_node[nodes_in_tree[i]] = false;
        }
        
        
        on_tree.resize(graph.num_nodes());
        even_node.resize(graph.num_nodes());
        odd_node.resize(graph.num_nodes());
        parent.resize(graph.num_nodes());
        visited.resize(graph.num_nodes());
        edges_in_cycle.clear();
        

        dsu.clean(nodes_in_tree);
        nodes_in_tree.clear();
        node_closes_cycle.clear();
        
    }

    bool Edmonds::exposed_vertex(size_type node){
        return matched_to[node] == node;
    }


    void Edmonds::extend_tree(size_type node, size_type even_parent){
        size_type match = matched_to[node];
        
        odd_node[node] = true;
        odd_node[match] = false;
        even_node[match] = true;
        even_node[node] = false;

        parent[match] = node;
        parent[node] = even_parent;

        add_neighbors_to_pending_list( match );
        nodes_in_tree.push_back(node);
        nodes_in_tree.push_back(match);


        //The nodes are in the tree
        on_tree[node] = true;
        on_tree[match] = true;

    }

    void Edmonds::remove_not_contained_edges(std::vector<std::pair<size_type, size_type> > &edges){
        int idx = -1;
        for(size_type i = 0 ; i < edges.size() && idx < 0; i++)  {
            if(!visited[edges[i].first || !visited[edges[i].second]])
                idx = i;
        }
        if(idx > 0)
            edges.erase(edges.begin() + idx, edges.end()); 
    
    }

    void Edmonds::shrink(size_type node_u, size_type node_v){
        
        cycle_edge.push_back(std::make_pair(node_u, node_v));


        size_type num_nodes = dsu.num_nodes();
        dsu.add_node();
        dsu.join(node_u, node_v);
        dsu.set_parent(node_u, num_nodes);
         
        size_type r_cycle_node, l_cycle_node;
        r_cycle_node = node_u;
        l_cycle_node = node_v;
        
        std::vector<std::pair<size_type,size_type> > l_path, r_path;
        l_path.push_back(std::make_pair(node_v, node_u));

        size_type common_predecessor ;

        //Cycle until we found a node reached by both paths
        while(true){
            //if(matched_to[r_cycle_node] != r_cycle_node){
                visited[r_cycle_node] = true;
            //}
            //if(matched_to[l_cycle_node] != l_cycle_node){
                visited[l_cycle_node] = true;
            //}
            size_type parent_r = dsu.find(parent[r_cycle_node]);
            size_type parent_l = dsu.find(parent[l_cycle_node]);
            
            if(!exposed_vertex(r_cycle_node))
                r_path.push_back(std::make_pair(r_cycle_node,parent_t));
            if(!exposed_vertex(l_cycle_node))
                l_path.push_back(std::make_pair(l_cycle_node,parent_l));

            r_cycle_node = parent_r;
            l_cycle_node = parent_t;;

            //Found common node
            if(visited[r_cycle_node] || visited[l_cycle_node] ||
                r_cycle_node==l_cycle_node){


                visited[r_cycle_node] = true;
                visited[l_cycle_node] = true;


                size_type tmp_node = common_predecessor = r_cycle_node;
                
                if(visited[l_cycle_node]) tmp_node = common_predecessor = l_cycle_node;

                //Set as univisted nodes not in the cycle
                while(visited[tmp_node]){
                    visited[tmp_node] = false;
                    tmp_node = dsu.find(parent[tmp_node]);
                }

                visited[common_predecessor] = true;
                 //Remove edges not contained in the cycle
                 //This gives a total time of O(|E(C)|) for the
                 //shrinking step
                 remove_not_contained_edges(l_path);  
                 remove_not_contained_edges(r_path);
                 

                 //ordering the edges in the cycle is just one
                 //of the paths reversed joined to the other

                 //Possible Improvement, just add to the left path
                 reverse(r_path.begin(), r_path.end());
                 std::vector<std::pair<size_type, size_type> > cycle;
                 cycle.reserve(l_path.size() + r_path.size());
                 cycle.insert(cycle.end(), l_path.begin(), l_path.end());
                 cycle.insert(cycle.end(), r_path.begin(), r_path.end());

                 edges_in_cycle.push_back(cycle);
                 break;
                  

            }

        }

        //Updates every vector for the tree
        on_tree.push_back(true);
        even_node.push_back(true);
        odd_node.push_back(false);
        
        visited.push_back(false);
        node_closes_cycle.push_back(common_predecessor);
        parent.push_back(parent[common_predecessor]);
    }

    void Edmonds::augment_matching(size_type node, size_type par){
        parent[node] = par;
    
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
               augment_matching(repr_y, repr_x); 
               reset_vectors();
               return;

           }
           //Extend tree
           else if(!on_tree[repr_y] && !exposed_vertex(repr_y)){
                extend_tree(repr_y, repr_x);     
           }
           //Shrink
           else{
                shrink(repr_x, repr_y);         
                size_type last_cycle = edges_in_cycle.size() - 1;
                for(size_type i =0 ; i < edges_in_cycle[last_cycle].size(); i++){
                    size_type node_a = edges_in_cycle[last_cycle][i].first;
                    size_type node_b = edges_in_cycle[last_cycle][i].second;

                    if(visited[node_a] && odd_node[node_a]){
                        add_neighbors_to_pending_list(node_a);        
                    }
                        
                    visited[node_a] = false;

                    if(visited[node_b] && odd_node[node_b]){
                        add_neighbors_to_pending_list(node_b);        
                    }
                    visited[node_b] = false;

                }
           }

        } 

       //Frustrated tree, mark nodes as frustrated
       
        for(size_type i = 0 ; i < nodes_in_tree.size(); i++)  {
            if(nodes_in_tree[i] < num_original_nodes)
                frustrated[nodes_in_tree[i]] = true;
        }
        reset_vectors();

    }
    //Edmonds Algrotihm implementation
    void Edmonds::run(){
       reset_vectors(); 
      
       for(size_type i =static_cast<size_type>(0) ; i < num_original_nodes; i++)
            if(exposed_vertex(i) && !frustrated[i])
                grow_tree(i);

    }

    
    //Function that matches node_u to node_v
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

            if(!on_tree[representative] && !frustrated[representative] && !frustrated[node]){
                add_edge_to_pending_list(node, neighbour); 
                on_tree[representative] = true;
            }
        }
    }
};
