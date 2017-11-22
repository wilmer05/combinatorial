#include"edmonds.hpp"
#include <utility>
#include <vector>

namespace ALG{

    Edmonds :: Edmonds(size_type num_nodes) : graph(num_nodes), dsu(num_nodes), on_tree(num_nodes, false), even_node(num_nodes, false),odd_node(num_nodes, false),  parent(num_nodes),  actual_node_to_parent(num_nodes), frustrated(num_nodes, false) {

        num_original_nodes = num_nodes;
        matched_to.resize(num_nodes);
        belongs_to_pseudonode.resize(num_nodes);
        for(size_type i =static_cast<size_type>(0) ; i < num_nodes;i++){
            matched_to[i] = i;
            belongs_to_pseudonode[i] = i;
            parent[i] = i;
            actual_node_to_parent[i] = i;
        }
         

        nodes_in_tree.clear();
        edges_in_cycle.clear();
        augmented.clear();
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
            actual_node_to_parent[nodes_in_tree[i]] = nodes_in_tree[i];
            belongs_to_pseudonode[nodes_in_tree[i]] = nodes_in_tree[i];
            on_tree[nodes_in_tree[i]] = false;
            even_node[nodes_in_tree[i]] = false;
            odd_node[nodes_in_tree[i]] = false;
        }
        
        
        on_tree.resize(graph.num_nodes());
        even_node.resize(graph.num_nodes());
        odd_node.resize(graph.num_nodes());
        parent.resize(graph.num_nodes());
        actual_node_to_parent.resize(graph.num_nodes());
        visited.resize(graph.num_nodes());
        belongs_to_pseudonode.resize(graph.num_nodes());
        edges_in_cycle.clear();
        pseudonodes_edges.clear();
        augmented.clear();
        

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

        actual_node_to_parent[node] = node;
        actual_node_to_parent[match] = match;

        add_neighbors_to_pending_list( match );
        if(!on_tree[node])
            nodes_in_tree.push_back(node);

        if(!on_tree[match])
            nodes_in_tree.push_back(match);

        //The nodes are in the tree
        on_tree[node] = true;
        on_tree[match] = true;

    }

    void Edmonds::remove_not_contained_edges(std::vector<std::pair<size_type, size_type> > &edges){
        int idx = -1;
        for(size_type i = 0 ; i < edges.size() && idx < 0; i++)  {
            size_type parent_a = dsu.find(edges[i].first);
            size_type parent_b = dsu.find(edges[i].second);
            if(!visited[parent_a] || !visited[parent_b])
                idx = i;
        }
        if(idx >= 0)
            edges.erase(edges.begin() + idx, edges.end()); 
    
    }

    void Edmonds::shrink(size_type node_u, size_type node_v){
        
        cycle_edge.push_back(std::make_pair(node_u, node_v));


        size_type num_nodes = dsu.num_nodes();
        dsu.add_node();

        //dsu.join(node_u, node_v);
        //dsu.set_parent(node_u, num_nodes);
         
        size_type r_cycle_node, l_cycle_node;
        r_cycle_node = dsu.find(node_u);
        l_cycle_node = dsu.find(node_v);
        
        visited[r_cycle_node] = true;
        visited[l_cycle_node] = true;

        std::vector<std::pair<size_type,size_type> > l_path, r_path, l_path_pseudonodes, r_path_pseudonodes;
        l_path.push_back(std::make_pair(node_u, node_v));
        l_path_pseudonodes.push_back(std::make_pair(r_cycle_node, l_cycle_node));
        

        size_type common_predecessor;

        //Cycle until we found a node reached by both paths
        while(true){
            size_type parent_r = dsu.find(parent[r_cycle_node]);
            size_type parent_l = dsu.find(parent[l_cycle_node]);
            
            if(!exposed_vertex(r_cycle_node)){
                r_path.push_back(std::make_pair(actual_node_to_parent[r_cycle_node], parent[r_cycle_node]));
               r_path_pseudonodes.push_back(std::make_pair(r_cycle_node, dsu.find(parent[r_cycle_node]))); 
            }
            if(!exposed_vertex(l_cycle_node)){
                l_path.push_back(std::make_pair(actual_node_to_parent[l_cycle_node], parent[l_cycle_node]));
               l_path_pseudonodes.push_back(std::make_pair(l_cycle_node, dsu.find(parent[l_cycle_node]))); 
            }
            r_cycle_node = parent_r;
            l_cycle_node = parent_l;

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
                 remove_not_contained_edges(l_path_pseudonodes);
                 remove_not_contained_edges(r_path_pseudonodes);
                 

                 //ordering the edges in the cycle is just one
                 //of the paths reversed joined to the other

                 //Possible Improvement, just add to the left path
                 reverse( r_path.begin(), r_path.end());
                 reverse( r_path_pseudonodes.begin(), r_path_pseudonodes.end());
                 std::vector<std::pair<size_type, size_type> > cycle, pseudonode_cycle;
                 cycle.reserve( l_path.size() + r_path.size());
                 pseudonode_cycle.reserve(l_path_pseudonodes.size() + r_path_pseudonodes.size());
                 cycle.insert(cycle.end(), l_path.begin(), l_path.end());
                 cycle.insert(cycle.end(), r_path.begin(), r_path.end());
                 pseudonode_cycle.insert(pseudonode_cycle.end(), l_path_pseudonodes.begin(), l_path_pseudonodes.end());
                 pseudonode_cycle.insert(pseudonode_cycle.end(), r_path_pseudonodes.begin(), r_path_pseudonodes.end());


                 //Join the partition set of vertices
                 for(size_type i = 0 ; i < pseudonode_cycle.size(); i++) {
                    size_type parent_a = pseudonode_cycle[i].first;
                    size_type parent_b = pseudonode_cycle[i].second;

                    belongs_to_pseudonode[parent_a] = belongs_to_pseudonode[parent_b] = num_nodes;
                    dsu.set_parent(parent_a, num_nodes);
                    dsu.set_parent(parent_b, num_nodes);
                    
                 }
                 edges_in_cycle.push_back(cycle);
                 break;

            }

        }

        //Updates every vector for the tree
        on_tree.push_back(true);
        even_node.push_back(true);
        odd_node.push_back(false);
        
        visited.push_back(false);
        belongs_to_pseudonode.push_back(num_nodes);
        node_closes_cycle.push_back(common_predecessor);

        parent.push_back(dsu.find(parent[common_predecessor]));
        actual_node_to_parent.push_back(actual_node_to_parent[dsu.find(common_predecessor)]);
        augmented.push_back(false);
    }

    bool Edmonds::is_pseudonode(size_type node){
        return dsu.find(node) >= num_original_nodes;
    }


    size_type Edmonds::get_next_node_to_match(size_type node, int incident_kind_of_edge){
    
        //if the node has to be matched, then this is the next node
        if(incident_kind_of_edge == 1) 
            return node;

        size_type pseudonode = dsu.find(node);
        return actual_node_to_parent[pseudonode];

    }

    bool Edmonds::node_belongs_to_edge(std::pair<size_type, size_type> edge, size_type node){
        return edge.first == node || edge.second == node; 
    }

    void Edmonds::augment_matching_on_cycle(size_type pseudonode, size_type ignored_node){
        int  root_edge_index= -1;
        int pseudo_index = pseudonode - num_original_nodes;

        if(augmented[pseudo_index])
            return;

        augmented[pseudo_index] = true;
        int cycle_size = pseudonodes_edges[pseudo_index].size();
        //First we find the ignored node
        for(int i = 0 ; i < cycle_size && root_edge_index < 0; i++){
           if(node_belongs_to_edge(pseudonodes_edges[pseudo_index][i], ignored_node)){
                root_edge_index = i; 
           } 
        }
        if(node_belongs_to_edge(pseudonodes_edges[pseudo_index][(root_edge_index + 1) % cycle_size], ignored_node)) 
            root_edge_index = (root_edge_index + 1) % cycle_size;
        int second_root_edge_index = root_edge_index - 1;
        if(second_root_edge_index < 0) 
            second_root_edge_index = cycle_size - 1;

        int current_edge_to_match = (root_edge_index + 1) % cycle_size;

        //Then match the correct nodes

        while(current_edge_to_match != root_edge_index && current_edge_to_match != second_root_edge_index){
            
            size_type node_u = edges_in_cycle[pseudo_index][current_edge_to_match].first; 
            size_type node_v = edges_in_cycle[pseudo_index][current_edge_to_match].second; 

            match(node_u, node_v);
            match_all_pseudonodes_of(node_u, 1);
            match_all_pseudonodes_of(node_v, 1);
            current_edge_to_match += 2;
            current_edge_to_match %= cycle_size;
        } 

    }





    void Edmonds::match_all_pseudonodes_of(size_type node, int next_kind_of_edge){
        size_type node_to_match = get_next_node_to_match(node, next_kind_of_edge);
        size_type last_pseudonode = dsu.find(node_to_match); 
        size_type ignored_node = node;
        while( ignored_node != last_pseudonode ){
            augment_matching_on_cycle(belongs_to_pseudonode[ignored_node], ignored_node);
            ignored_node = belongs_to_pseudonode[ignored_node];
        } 
    }

    void Edmonds::augment_matching(size_type node, size_type par){
        parent[node] = par;
        size_type root_pseudonode = dsu.find(current_root);

        //This variable indicates the next kind of edge
        //it's 0 for an unmatched edge
        //and 1 for a matched edge
        int next_kind_of_edge = 1;
        while(dsu.find(node) != root_pseudonode){

            //The node is inside a cycle
            if(node != dsu.find(node)) {
                match_all_pseudonodes_of(node, next_kind_of_edge);        
            } else if(next_kind_of_edge){
                matched_to[node] = parent[node];
                matched_to[parent[node]] = node;
            }
            
            next_kind_of_edge = 1 - next_kind_of_edge;
            node = parent[dsu.find(node)];
        }
        //augment matching on the cycle of the root
        augment_matching_on_cycle(root_pseudonode, current_root);
    
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

           if(repr_x == repr_y) 
                continue;

           if(!even_node[repr_x] || odd_node[repr_y]){
                std::swap(repr_x, repr_y);
                if( !even_node[repr_x] || odd_node[repr_y]) 
                    continue;
           }

           //Condition 1: augment matching
           if(!on_tree[repr_y] && exposed_vertex(repr_y)){
               augment_matching(node_y, node_x); 
               reset_vectors();
               return;

           }
           //Extend tree
           else if(!on_tree[repr_y] && !exposed_vertex(repr_y)){
                extend_tree(node_y, node_x);     
           }
           //Shrink
           else{
                shrink(node_x, node_y);         
                size_type last_cycle = edges_in_cycle.size() - 1;


                //We add edges from odd nodes
                //Every node is actually two times, so we care that
                //we don't add unnecesary edges
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
            if(exposed_vertex(i) && !frustrated[i]){
                current_root  = i;
                grow_tree(i);
            }

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
            size_type repr_neigh = dsu.find(neighbour);
            size_type repr_node = dsu.find(node);
                
            if(repr_neigh != repr_node && !frustrated[neighbour] && !frustrated[node]){
                add_edge_to_pending_list(node, neighbour); 
                on_tree[repr_neigh] = true;
            }
        }
    }
};
