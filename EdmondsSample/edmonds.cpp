#include"edmonds.hpp"
#include <utility>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace ALG{

    Edmonds :: Edmonds(size_type num_nodes) : graph(num_nodes), dsu(num_nodes), on_tree(num_nodes, false), even_node(num_nodes, false),odd_node(num_nodes, false),  parent(num_nodes),  actual_node_to_parent(num_nodes), frustrated(num_nodes, false) {

        //We initialize every variable in our class
        num_original_nodes = num_nodes;
        matched_to.resize(num_nodes);
        belongs_to_pseudonode.resize(num_nodes);

        //Initially each node is set to be matched to itself (this
        //actually means that it's not matched to other node
        //The parent of a node is also itself
        //And the node that connects a node to its parent is
        //itself
        for(size_type i = 0 ; i < num_nodes;i++){
            matched_to[i] = i;
            belongs_to_pseudonode[i] = i;
            parent[i] = i;
            actual_node_to_parent[i] = i;
        }
         

        nodes_in_tree.clear();
        edges_in_cycle.clear();
        augmented.clear();
    }

    Edmonds::~Edmonds(){
        //We clear every vector and structures inside our class
        dsu.clean();
        while(pending_edges.size())
            pending_edges.pop();

        on_tree.clear();
        odd_node.clear();
        even_node.clear();
        parent.clear();
        actual_node_to_parent.clear();
        matched_to.clear();
        visited.clear();
        edges_in_cycle.clear();
        pseudonodes_edges.clear();
        nodes_in_tree.clear();
        frustrated.clear();
        belongs_to_pseudonode.clear();
        augmented.clear();
    }

    void Edmonds::find_maximal_matching(){

        //We iterate over every node and if it is possible to match it
        //we match it
        for(size_type node = 0 ; node < graph.num_nodes(); node++){
            if(!exposed_vertex(node))
                continue;
            for(size_type neighbour : graph.node( node ).neighbors()){
                if(exposed_vertex(neighbour)){
                    match(node, neighbour);
                    break;
                }
            }
        } 
    }

    void Edmonds::reset_vectors(){
        while(pending_edges.size())
            pending_edges.pop();

        //Reset the vectors only for nodes in the current alternating
        //tree since the vector entries for other nodes have not been
        //changed
        for(size_type i = static_cast<size_type>(0); i < nodes_in_tree.size(); i++){
            parent[nodes_in_tree[i]] = nodes_in_tree[i];
            actual_node_to_parent[nodes_in_tree[i]] = nodes_in_tree[i];
            belongs_to_pseudonode[nodes_in_tree[i]] = nodes_in_tree[i];
            on_tree[nodes_in_tree[i]] = false;
            even_node[nodes_in_tree[i]] = false;
            odd_node[nodes_in_tree[i]] = false;
        }
        
        
        //We resize the vectors to the number of nodes in the 
        //original graph

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
        

        //We also clean the DSU structures only for the used nodes
        dsu.clean(nodes_in_tree);
        nodes_in_tree.clear();
       
    }

    bool Edmonds::exposed_vertex(size_type node){
        //A node is exposed if it's matched to itself
        return matched_to[node] == node;
    }


    void Edmonds::extend_tree(size_type node, size_type even_parent){
        size_type match = matched_to[node];
        
        //We extend the tree with node being odd and the match
        //as an even node
        odd_node[node] = true;
        odd_node[match] = false;
        even_node[match] = true;
        even_node[node] = false;

        parent[match] = node;
        parent[node] = even_parent;

        actual_node_to_parent[node] = node;
        actual_node_to_parent[match] = match;

        //We add the incident edges to the matched node
        add_incident_edges_to_pending_list( match );

        //If the nodes weren't in the tree, then we add them
        //To the vector of used nodes
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

        //The visited vector is marked wieh the nodes that are actually in the cycle
        //So we remove from the list of edges passed as arguent those edges
        //with a non-visited node

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

        //We add a new node to the DSU, because we will add a 
        //new pseudonode
        size_type num_nodes = dsu.num_nodes();
        dsu.add_node();
         
        /**
            r_cycle_node refers to the node in the right path 
            of the cycle. l_cycle_node is analogue but on the 
            left path
        **/
        size_type r_cycle_node, l_cycle_node;

        //We look for the pseudonodes corresponding to the nodes w
        //hich close the cycle

        r_cycle_node = dsu.find(node_u);
        l_cycle_node = dsu.find(node_v);
        
        //We mark as visited every vertex in the construction
        //of the left path and right path
        visited[r_cycle_node] = visited[l_cycle_node] = true;

        //vector of edges corresponding to the 
        //left path (l_path) and right path (r_path)
        //to construct the cycle (they connect nodes in the
        //original graph)
        std::vector<std::pair<size_type,size_type> > l_path, r_path;
        
        //vector of edges corresponding to the 
        //left path (l_path) and right path (r_path)
        //to construct the cycle (they connect maybe pseudonodes
        //and maybe nodes)
        std::vector<std::pair<size_type,size_type> > l_path_pseudonodes, r_path_pseudonodes;


        //The first edge in the left path is the one that closes
        //the cycle, after that l_path and r_path are disjoint
        l_path.push_back(std::make_pair(node_u, node_v));
        l_path_pseudonodes.push_back(std::make_pair(r_cycle_node, l_cycle_node));

        //Now we have to find the common predecessor
        //to close the cycle and build the pseudonode
        size_type common_predecessor;
        size_type root_repr = dsu.find(current_root);

        //Iterate until we found a node reached by both paths
        //which is the commmon_predecessor
        while(true){
            size_type parent_r = dsu.find(parent[r_cycle_node]);
            size_type parent_l = dsu.find(parent[l_cycle_node]);
            visited[r_cycle_node] = true; 
            visited[l_cycle_node] = true;
            

            //We add the edge on the right side of the cycle
            //if the edge is not from the root to itself
            if(dsu.find(parent[r_cycle_node]) != r_cycle_node){
                r_path.push_back(std::make_pair(actual_node_to_parent[r_cycle_node], parent[r_cycle_node]));
                r_path_pseudonodes.push_back(std::make_pair(r_cycle_node, parent_r)); 
            
            }

            //We add the edge on the left side of the cycle
            //if the edge is not from the root to itself
            if(dsu.find(parent[l_cycle_node]) != l_cycle_node){
                l_path.push_back(std::make_pair(actual_node_to_parent[l_cycle_node], parent[l_cycle_node]));
                l_path_pseudonodes.push_back(std::make_pair(l_cycle_node, parent_l)); 
            }

            r_cycle_node = parent_r;
            l_cycle_node = parent_l;

            //If we find a common predecessor which closes the cycle
            //it means that the paths that we are following has 
            //already visited the predecessor by one of the path
            //or both path are in the same node now
            if((visited[r_cycle_node] && r_cycle_node != root_repr) 
            || (visited[l_cycle_node] && l_cycle_node != root_repr)  
            || r_cycle_node == l_cycle_node){

                
                //We create a temportal node to set as 
                //unvisited nodes not in the cycle
                size_type tmp_node = common_predecessor = r_cycle_node;
                
                if(visited[l_cycle_node] && l_cycle_node != root_repr)
                    tmp_node = common_predecessor = l_cycle_node;

                visited[r_cycle_node] = true;
                visited[l_cycle_node] = true;

                //Set as univisted nodes not in the cycle
                while(visited[tmp_node]){
                    visited[tmp_node] = false;
                    tmp_node = dsu.find(parent[tmp_node]);
                }

                 //Visited means that the node is in the cycle
                 //of course the common_predecessor is in it
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
                 //So we reverse the right path
                 reverse( r_path.begin(), r_path.end());
                 reverse( r_path_pseudonodes.begin(), r_path_pseudonodes.end());
                 //Vector of edges in the cycle (new pseudonode)
                 std::vector<std::pair<size_type, size_type> > cycle, pseudonode_cycle;

                 cycle.reserve( l_path.size() + r_path.size());
                 pseudonode_cycle.reserve(l_path_pseudonodes.size() + r_path_pseudonodes.size());

                 //We join the paths in the cycle
                 cycle.insert(cycle.end(), l_path.begin(), l_path.end());
                 cycle.insert(cycle.end(), r_path.begin(), r_path.end());

                 //We join the paths in the cycle with pseudonodes
                 pseudonode_cycle.insert(pseudonode_cycle.end(), l_path_pseudonodes.begin(), l_path_pseudonodes.end());
                 pseudonode_cycle.insert(pseudonode_cycle.end(), r_path_pseudonodes.begin(), r_path_pseudonodes.end());

                 //Finally we join the partition sets in the dsu  
                 //of vertices in the new pseudonode
                 join_partition_sets(pseudonode_cycle, num_nodes);

                 //We add the vector of edges
                 //corresponding to this new pseudonode
                 edges_in_cycle.push_back(cycle);
                 pseudonodes_edges.push_back(pseudonode_cycle);
                 break;

            }

        }

        //Adds the new pseudonode to the arrays in our algorithm
        add_new_pseudonode(common_predecessor, num_nodes);
        
    }

    void Edmonds::join_partition_sets(std::vector<std::pair<size_type, size_type> > &cycle, size_type num_nodes){
        for(size_type i = 0 ; i < cycle.size(); i++) {
             size_type parent_a = cycle[i].first;
             size_type parent_b = cycle[i].second;

             belongs_to_pseudonode[parent_a] = belongs_to_pseudonode[parent_b] = num_nodes;
             dsu.set_parent(parent_a, num_nodes);
             dsu.set_parent(parent_b, num_nodes);
             
        }
    }


    void Edmonds::add_new_pseudonode(size_type common_predecessor, size_type num_nodes){
        //Updates every vector for the tree
        on_tree.push_back(true);
        even_node.push_back(true);
        odd_node.push_back(false);
        
        visited.push_back(false);

        //the new pseudonode belongs to itself
        belongs_to_pseudonode.push_back(num_nodes);


        if(dsu.find(common_predecessor) == dsu.find(current_root))
            parent.push_back(num_nodes);
        else
            parent.push_back(parent[common_predecessor]);

        actual_node_to_parent.push_back(actual_node_to_parent[common_predecessor]);
        augmented.push_back(false);


    }


    bool Edmonds::is_pseudonode(size_type node){
        return dsu.find(node) >= num_original_nodes;
    }



    size_type Edmonds::get_next_node_to_match(size_type node, int incident_kind_of_edge){
    
        //if the node has to be matched, then this is the next node
        if(incident_kind_of_edge != 1) 
            return node;

        size_type pseudonode = dsu.find(node);
        return actual_node_to_parent[pseudonode];

    }

    bool Edmonds::node_belongs_to_edge(std::pair<size_type, size_type> edge, size_type node){
        return edge.first == node || edge.second == node; 
    }


    void Edmonds::augment_matching_on_cycle(size_type pseudonode, size_type ignored_node){
        int pseudo_index = pseudonode - num_original_nodes;

        //It's not a pseudonode
        if(pseudo_index < 0)
            return;

        //If we already agumented this cycle, then we go back
        if(augmented[pseudo_index])
            return;

        //else we mark this cycle as augmented
        augmented[pseudo_index] = true;
        int cycle_size = pseudonodes_edges[pseudo_index].size();

        //First we find the right most edge where the ignored node
        //belongs to
        //Idex of the first edge containing the ignored node
        long  root_edge_index= -1;
        for(long  i = 0 ; i < cycle_size && root_edge_index < 0; i++){
           if(node_belongs_to_edge(pseudonodes_edges[pseudo_index][i], ignored_node)){
                root_edge_index = i; 
           } 
        }

        //If the next edge contains the ignored node then it is 
        //the rightmost 
        if(node_belongs_to_edge(pseudonodes_edges[pseudo_index][(root_edge_index + 1) % cycle_size], ignored_node)) 
            root_edge_index = (root_edge_index + 1) % cycle_size;

        //Index of the other edge containing the ignored node
        int second_root_edge_index = root_edge_index - 1;
        if(second_root_edge_index < 0) 
            second_root_edge_index = cycle_size - 1;

        //We now match every two edges in the cycle
        int current_edge_to_match = (root_edge_index + 1) % cycle_size;

        while(current_edge_to_match != root_edge_index && current_edge_to_match != second_root_edge_index){
            size_type node_u = edges_in_cycle[pseudo_index][current_edge_to_match].first; 
            size_type node_v = edges_in_cycle[pseudo_index][current_edge_to_match].second; 
            
            //We match the nodes in the edge that are part of the
            //original graph G (without pseudonodes)
            match(node_u, node_v);

            //We have also to find a perfect matching
            //In every pseudonode containing these nodes
            match_all_pseudonodes_of(node_u, 0);
            match_all_pseudonodes_of(node_v, 0);


            current_edge_to_match += 2;
            current_edge_to_match %= cycle_size;
        } 

    }



    void Edmonds::match_all_pseudonodes_of(size_type node, int next_kind_of_edge){

        //Given the kind of edge connecting the pseudonode where 
        //node is we have to match the node given as parameter
        //or the node connecting the pseudonode to it's parent
        size_type node_to_match = get_next_node_to_match(node, next_kind_of_edge);


        size_type last_pseudonode = dsu.find(node_to_match); 

        //If this is not a pseudonode then we return
        if(last_pseudonode < num_original_nodes) 
            return;

        //Otherwise we recurse in every pseudonode containing 
        //the node matched to find a perfect matching
        size_type ignored_node = node_to_match;
        while( ignored_node != last_pseudonode && !augmented[belongs_to_pseudonode[ignored_node] - num_original_nodes]){
            augment_matching_on_cycle(belongs_to_pseudonode[ignored_node], ignored_node);
            ignored_node = belongs_to_pseudonode[ignored_node];
        } 
    }

    void Edmonds::augment_matching(size_type node, size_type par){
        //We set the parent of the exposed node
        parent[node] = par;
        on_tree[node] = true;
        odd_node[node] = true;
        nodes_in_tree.push_back(node);

        size_type root_pseudonode = dsu.find(current_root);


        //This variable indicates the next kind of edge
        //it's 0 for an unmatched edge
        //and 1 for a matched edge
        int next_kind_of_edge = 1;

        size_type current_pseudonode, last_pseudonode;
        while((current_pseudonode = dsu.find(node)) != root_pseudonode){
            //If the node is inside a cycle
            //then we have to recurse and try to find
            //a perfect matching in every pseudonode where
            //it is, and match it if the kind of edge 
            //that must go to the parent is non-matched
            if(node != current_pseudonode) {
                match_all_pseudonodes_of(node, next_kind_of_edge);        
                if(!next_kind_of_edge){
                    match(node, actual_node_to_parent[last_pseudonode]);
                    match_all_pseudonodes_of(node, 0);
                    match_all_pseudonodes_of(actual_node_to_parent[last_pseudonode], 0);
                    //match(node, parent[node]);
                } else {
                    match(parent[current_pseudonode], actual_node_to_parent[current_pseudonode]);
                    match_all_pseudonodes_of(parent[current_pseudonode], 0);
                    match_all_pseudonodes_of(actual_node_to_parent[current_pseudonode], 0);
                }

            //Else the node is not inside a pseudonode
            //So we could match it to the last node or to its parent
            } else if(next_kind_of_edge){
                match(node, parent[node]);
                match_all_pseudonodes_of(parent[node], 0);
            } else{
                match(actual_node_to_parent[last_pseudonode], node);
                match_all_pseudonodes_of(actual_node_to_parent[last_pseudonode], 0);
            }



            last_pseudonode = current_pseudonode; 
            next_kind_of_edge = 1 - next_kind_of_edge;
            node = parent[current_pseudonode];
        }

        //augment matching on the cycle of the root
        match_all_pseudonodes_of(parent[last_pseudonode], 0);
            
    }


    void Edmonds::add_incident_edges_of_last_cycle_odd_nodes(){
         size_type last_cycle = edges_in_cycle.size() - 1;


         //We add edges from odd nodes
         //Every node is actually two times, so we care that
         //we don't add unnecesary edges
         for(size_type i =0 ; i < pseudonodes_edges[last_cycle].size(); i++){
             size_type node_a = pseudonodes_edges[last_cycle][i].first;
             size_type node_b = pseudonodes_edges[last_cycle][i].second;

             if(!visited[node_a] && odd_node[node_a]){
                 add_incident_edges_to_pending_list(node_a);        
             }
                 
             visited[node_a] = false;

             if(!visited[node_b] && odd_node[node_b]){
                 add_incident_edges_to_pending_list(node_b);        
             }
             visited[node_b] = false;

         }

    }


    void Edmonds::grow_tree(size_type root_node){
        on_tree[root_node] = true;
        add_incident_edges_to_pending_list(root_node);
        even_node[root_node] = true;
        nodes_in_tree.push_back(root_node);

        //While there are some edges that have not been considered yet
        //we stay in this while
        while(pending_edges.size()){
           std::pair<size_type, size_type> edge = pending_edges.front();
           pending_edges.pop();

           size_type node_x = edge.first;
           size_type node_y = edge.second;
          
           //We find the representatives of the endpoints of the edge
           size_type repr_x = dsu.find(node_x);
           size_type repr_y = dsu.find(node_y);

           //If they belong to the same pseudonode then we don't do
           //anything
           if(repr_x == repr_y) 
                continue;

           //If the edge is not going from a even node to
           //an exposed vertex or to other even node
           //we continue
           if(!even_node[repr_x] || odd_node[repr_y]){
                std::swap(repr_x, repr_y);
                if( !even_node[repr_x] || odd_node[repr_y]) 
                    continue;
           }

           //Condition 1: augment matching
           if(!on_tree[repr_y] && exposed_vertex(node_y)){
               augment_matching(node_y, node_x); 
               reset_vectors();
               return;

           }
           //Extend tree
           else if(!on_tree[repr_y] && !exposed_vertex(node_y)){
                extend_tree(node_y, node_x);     
           }
           //Shrink since the edge is between two even nodes
           else{
                shrink(node_x, node_y);         
                add_incident_edges_of_last_cycle_odd_nodes();
           }
        } 


        //We didn't find an exposed vertex fo the tree is Frustrated
        //we mark nodes the nodes in the current treeas frustrated
        for(size_type i = 0 ; i < nodes_in_tree.size(); i++)  {
            if(nodes_in_tree[i] < num_original_nodes)
                frustrated[nodes_in_tree[i]] = true;
        }
        reset_vectors();
    }
    //Edmonds Algrotihm implementation
    void Edmonds::run(){

       //first we find an initial maximal matching (even if we 
       //receive a hint, we try to augment it
       find_maximal_matching();
       reset_vectors(); 
      
       for(size_type i =static_cast<size_type>(0) ; i < num_original_nodes; i++){
            if(exposed_vertex(i) && !frustrated[i]){
                current_root  = i;
                grow_tree(i);
            }
       }

    }
    
    //Function that matches node_u to node_v
    void Edmonds::match(size_type node_u, size_type node_v){
          matched_to[node_u] = node_v;
          matched_to[node_v] = node_u;
    }

    inline void Edmonds::add_edge_to_pending_list(size_type node_u, size_type node_v){
        pending_edges.push(std::make_pair(node_u, node_v));
    }

    inline void Edmonds::add_incident_edges_to_pending_list(size_type node){
        for(auto const &neighbour : graph.node( node ).neighbors()){
            size_type repr_neigh = dsu.find(neighbour);
            size_type repr_node = dsu.find(node);
                
            if(repr_neigh != repr_node && !frustrated[neighbour] && !frustrated[node]){
                add_edge_to_pending_list(node, neighbour); 
            }
        }
    }

    void Edmonds::print_matching(){

        std::vector<std::pair<size_type, size_type> > matching;

        for(size_type i = 0  ; i < num_original_nodes; i++)
            if(matched_to[i] > i)
                matching.push_back(std::make_pair(i+1, matched_to[i] + 1));
        
        std::cout << "p edge " << num_original_nodes << " " << matching.size() << "\n";
        for(size_type i=0 ; i < matching.size(); i++){
            std::cout << "e " << matching[i].first << " " << matching[i].second;
            if(i + 1 < matching.size())
                std::cout <<"\n";

        }
    }
};
