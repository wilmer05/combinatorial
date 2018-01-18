
#include<iostream>
#include<algorithm>
#include<utility>
#include<queue>
#include<cmath>
#include<stack>
#include"graph.hpp"
#include"dsu.hpp"
#include"held_karp.hpp"
namespace ALGORITHM{ //Start of namespace ALGORITHM

    void SearchNode :: print(){
        std::cout << "Req -> \n";
        for(size_type i =0 ; i < R.size(); i++) {
            for(size_type j =0 ; j < R[i].size();j++){
                std::cout << i << " " << R[i][j] << "\n";
            }
        }

        std::cout << "########\n F -> \n";
        for(size_type i = 0 ; i < F.size();i++)
            std::cout << F[i].first << " " << F[i].second <<"\n";
        std::cout << "Last one tree\n";
        for(size_type i =0 ; i < last_1_tree.size(); i++)
            for(size_type j =0 ; j < last_1_tree[i].size(); j++)
                if(i < last_1_tree[i][j])
                    std::cout << i << " edge " << last_1_tree[i][j] << "\n";
        std::cout << "--------------------" << std::endl;
        std::cout.flush();


    }

    bool SearchNode :: can_require(std::vector<std::vector<size_type> > &R, size_type node_u, size_type node_v){
        return !std::count(R[node_u].begin(), R[node_u].end(), node_v) && R[node_v].size() < 2 && R[node_u].size() < 2;
    }

    bool SearchNode :: not_required(std::vector<size_type> &R, size_type node_v){
        return !std::count(R.begin(), R.end(), node_v);
    }

    bool SearchNode :: check_1_tree(){
	std::vector<bool> visited(last_1_tree.size(), false);
	std::priority_queue<size_type> q;
	q.push(0);
	visited[0] = true;
	size_type total_edges = 0;
	while(q.size()){
	    size_type nod = q.top();
	    q.pop();
	    for(size_type i =0 ; i < last_1_tree[nod].size(); i++){
		size_type neighour = last_1_tree[nod][i];
		total_edges++;
		if(!visited[neighour]){
		    visited[neighour] = true;
		    q.push(neighour);
		}
	    }
			
	}
	return !std::count(visited.begin(), visited.end(), false) && last_1_tree[0].size() == 2 && total_edges == 2 * visited.size();

    }

    bool SearchNode :: solution_is_2_regular(){
    
        bool regular_2 = true;
        for(size_type i = 0 ; regular_2 && i < last_1_tree.size(); i++){
            regular_2 = regular_2 && last_1_tree[i].size() == 2;
        } 

        return regular_2 && last_1_tree.size();
    }

    std::vector<SearchNode> SearchNode::get_children(){

        std::vector<SearchNode> children;

        //If we cannot get a proper tree with the forbidden edges
        //Then we return
        if(invalid_node)
            return children;
        //Tmp variable to generate every child
        SearchNode child;
	    child.last_total_cost = 0.0;
        child.R = R;
        child.F = F;
        child.lambda = lambda;
        child.invalid_node = false;
        child.root_node = false;

        bool node_with_required = false;

        //p node of the paper
        size_type node_p = 1e9;
        size_type max_counter_node = 1e9;
        std::vector<size_type> counter(lambda.size(), 0);
        for(size_type i = 0; i < F.size(); i ++){
            counter[F[i].first]++;
            counter[F[i].second]++;
        }

        //In case there is no node p incident to a required edge
        //then we have to select depending on the number of
        //smallest number of feasible incident edges
        //(Heuristic 2 in the paper)
        for(size_type i = 0 ; i < lambda.size() && node_p >= lambda.size(); i++){

            if(last_1_tree[i].size() > 2){

                if(!R[i].size() && (max_counter_node >= lambda.size()  || counter[i] > counter[max_counter_node]))
                    max_counter_node = i;

                if(R[i].size() == 1){
                    node_with_required = true;
                    node_p = i;
                    break;
                }
            }
        }

        if(!node_with_required){
            //Second heuristic
            //node with least feaseble incident edges
            node_p = max_counter_node;
        }

       if(node_p >= lambda.size())
            return children;

       for(size_type j =0 ; j < last_1_tree[node_p].size() && !children.size(); j++){

           size_type adj_node = last_1_tree[node_p][j];
           if(!can_require(R, node_p, adj_node))
               continue;
           
           //Child S3 in the paper, choosing 
           //e1 = {node_p, adj_node}
           child.F.push_back(std::make_pair(node_p, adj_node));
           children.push_back(child);
           child.F.pop_back();
           
           //We have to choose e2 now
           for(size_type k =0 ;  k < last_1_tree[node_p].size(); k++){
                size_type adj_node_2 = last_1_tree[node_p][k];

                //We just forbid this edge if it's not already 
                //required
                if(adj_node != adj_node_2 && not_required(R[node_p], adj_node_2)){
                    child.R[node_p].push_back(adj_node);
                    child.R[adj_node].push_back(node_p);

                    child.F.push_back(std::make_pair(node_p, adj_node_2));
                    children.push_back(child);
                    child.F.pop_back();

                    //It's possible to require e2 only if it's the
                    //second edge to require
                    if(1 == child.R[node_p].size() && can_require(R, node_p, adj_node_2)){
                        child.R[node_p].push_back(adj_node_2);
                        child.R[adj_node_2].push_back(node_p);
                        children.push_back(child);
                    }
                    break;
                }
           }
       } 
            
        

        /*if(children.size() < 2 || children.size() >=4){
            std::cout << node_p << "\n";
            print();
            std::cout << children.size() << "\n";
        }*/
        assert(children.size() >= 2 && children.size() < 4);
        return children;
    }

    /**
        @class HeldKarp
        @brief This class implements the algoritm

    **/
    void HeldKarp :: run_HK(SearchNode &node) {
        

        std::vector<std::pair<size_type, size_type> > &F = node.get_F();
        graph.reset_edges();
        graph.fix_forbidden_edges(F, node.get_R());
        size_type num_nodes = graph.num_nodes();

        //We initialize the Required edges to a vector of empty
        //vectors
        if(node.is_root_node())
            node.R = std::vector<std::vector<size_type> >(graph.num_nodes(), std::vector<size_type>());

        //N changes depending on wether this node is the root
        //or not
        double N = (num_nodes + 3) / 4 + 5;
        if(node.is_root_node()){
            N = (num_nodes * num_nodes + 49) / 50 + num_nodes + 15;
	    //N = 0.3 * N;
	}
	//if(node.is_root_node()) std :: cout << " " << N << "\n";
       
        //We restart the last_1_tree
        node.last_1_tree = std::vector<std::vector<size_type> >(num_nodes, std::vector<size_type>()); 

        double t0;
        double tstep;
        double delta0;
        double delta_delta;
        double deltai;

        for(size_type step = 0 ; step < N ; step++){
            //We sort the edges based on the current lambda
            graph.fix_lambdas_and_sort_edges(node.get_lambda());
            //compute_1_tree(node);
            prims(node);
            //node.print();
            if(node.check_1_tree() && node.solution_is_2_regular())
                break;

            if(!step){
                if(node.is_root_node())
                    t0 = ((double) node.last_total_cost) / (2.0 * num_nodes);
                else
                    t0 = sum_lambda_root; 
                
                delta0 = (3.0 * t0) / (2.0 * N);
                delta_delta = t0 / (N * N - N);
                deltai = delta0;
                tstep = t0;
            }
	    /*if(node.is_root_node()){
	    std:: cout << "tstep= " << tstep << " step=" << step;
	    std:: cout << deltai << " " << delta_delta << "\n";
	    for(size_type i =0 ; i < node.lambda.size() ; i++)
		std::cout << node.lambda[i] << "\n";}*/
            update_lambda_function(node, tstep);
            tstep -= deltai;
            deltai -= delta_delta;
        }

        if(node.is_root_node()){
            std::vector<double> &lambda = node.get_lambda();
            sum_lambda_root = 0.0;
            for(size_type i =0 ;i  < lambda.size(); i++)
                sum_lambda_root += fabs(lambda[i]);
            sum_lambda_root *= (0.5 / num_nodes);
        }

        node.last_total_cost *= (1.0 - epsilon);
    }

    void HeldKarp :: update_lambda_function(SearchNode &node, double tstep) {

        //Update lambda function based on the nodes degrees and the t's
        std::vector<double> &lambda = node.get_lambda();

        for(size_type i =0 ; i < lambda.size(); i++){
             lambda[i] += tstep * d_parameter * (1.0 * node.last_1_tree[i].size() - 2.0);

             //In the first step there is not a Tree from the 
             //-1 iteration
             lambda[i] += tstep * ((1 - d_parameter) * 1.0 * node.last_second_1_tree[i].size() - 2.0);
        } 
    
    }


    void HeldKarp :: output_ans() {
    

       //The solution should be a 2-regular graph
       for(size_type i =0 ; i < best_solution.size(); i++){
            assert(best_solution[i].size() == 2);
       }

       
       assert(best_solution.size() == graph.num_nodes());
       size_type last,current;
       last = 0;
       current = best_solution[0][0];

       double sol = graph.get_distance(last, current);

       //We run through the 2-regular graph until we get
       //back to the first node
       if(stdout != out){
            fprintf(out, "TYPE : TOUR\n");  
            fprintf(out, "DIMENSION : %lu\n", graph.num_nodes());  
            fprintf(out, "TOUR_SECTION\n");  
            fprintf(out, "%lu\n", 1UL);
       }
       while(current!=0){
            if(stdout != out)
                fprintf(out,"%lu\n", current + 1);
            size_type tmp = last;
            last = current;
            if(best_solution[current][0] == tmp)
                current = best_solution[current][1];
            else current = best_solution[current][0];

            sol += graph.get_distance(last, current);
       }
        

       if(stdout != out){
            fprintf(out,"-1\nEOF");
       }
       std::cout << sol << "\n";
    }

    double HeldKarp :: get_lambda_edge_cost(SearchNode &node, size_type u, size_type v){
        return graph.get_edge(u, v).get_dist() + node.lambda[u] + node.lambda[v];
    }

    void HeldKarp :: prims(SearchNode &node){
        node.last_second_1_tree = node.last_1_tree;
        size_type num_nodes = node.R.size();
        
        std::vector<size_type> best_edge_to_add(num_nodes, 1);
        std::vector<bool> in_tree(num_nodes, false);
        node.last_1_tree = std::vector<std::vector<size_type> > (num_nodes);
        std::vector<std::vector<size_type> > &R = node.get_R();
        in_tree[1] = true;
        node.last_total_cost = 0.0;
        for(;;){
            size_type connect_node = 1e9;
            double best_found = infinity;
            //we find the best cost edge to add this to the 1
            for(size_type j = 1; j < num_nodes; j++){
                if(in_tree[j]) continue;

                size_type neigh = best_edge_to_add[j];

                if(!node.not_required(R[j], neigh)){
                    connect_node = j;
                    break;
                }
                double cost = get_lambda_edge_cost(node, j, neigh) ; //graph.get_edge(j, neigh).get_dist() + node.lambda(j) + node.lambda(neigh); 
                if(cost < best_found){
                    best_found = cost;
                    connect_node = j;
                }
            }

            if(connect_node >= num_nodes){
                break;
            }

            size_type neigh = best_edge_to_add[connect_node];
            node.last_1_tree[neigh].push_back(connect_node);
            node.last_1_tree[connect_node].push_back(neigh);
            in_tree[connect_node] = true;
            //node.last_total_cost += best_found;
            node.last_total_cost += graph.get_edge(connect_node, neigh).get_dist();

            //Now we find if connecting other nodes to the new node connect node is better 
            for(size_type i =2 ; i < num_nodes; i++){
                if(in_tree[i]) continue;
                size_type neigh = best_edge_to_add[i];
                double cost = get_lambda_edge_cost(node, connect_node, i); //graph.get_edge(connect_node, i).get_dist() + node.lambda(i) + node.lambda(connect_node);
                double current_cost = get_lambda_edge_cost(node, neigh, i); //graph.get_edge(neigh, i).get_dist() + node.lambda(i) + node.lambda(neigh);

                if( !node.not_required(R[i], connect_node) || cost < current_cost )
                    best_edge_to_add[i] = connect_node;
            }
            
        }

        //We have to add now the best two edges from node 0
        //These are the index of the best two neighbours
        size_type idx0 = 1,idx1 = 2; 
        if(get_lambda_edge_cost(node, 0, 1) >= get_lambda_edge_cost(node, 0, 2) || !node.not_required(R[0], 2))
            std::swap(idx0, idx1);
        for(size_type i =3 ; i < R.size(); i++){
            if(!node.not_required(R[0], i)){ 
                idx1 = i;
                std::swap(idx1,idx0);
            }else {
                 double cost = get_lambda_edge_cost(node, 0, idx1); 
                 double cost_0 = get_lambda_edge_cost(node, 0, idx0); 
                 double new_cost = get_lambda_edge_cost(node, 0, i);
                 if(node.not_required(R[0], idx1) && cost > new_cost)
                    idx1 = i;
                 double tmp_cost = get_lambda_edge_cost(node, 0, idx1);
                 if(tmp_cost < cost_0 && node.not_required(R[0], idx0))
                    std::swap(idx1, idx0);
                 
            }
        }
        node.last_1_tree[0].push_back(idx0);
        node.last_1_tree[0].push_back(idx0);
        node.last_1_tree[idx1].push_back(0);
        node.last_1_tree[idx0].push_back(0);
        node.last_total_cost += graph.get_edge(0, idx0).get_dist();
        node.last_total_cost += graph.get_edge(0, idx1).get_dist();
        
        std::vector<double> &lambda = node.get_lambda();
        for(size_type i =0 ; i < lambda.size();i++){
            node.last_total_cost += lambda[i] * (node.last_1_tree[i].size() - 2.0);
        }

        if(!node.last_second_1_tree.size())
            node.last_second_1_tree = node.last_1_tree;
        //node.print();
    }

    void HeldKarp :: compute_1_tree(SearchNode &node){
        node.last_second_1_tree = node.last_1_tree;
        DSU :: Dsu dsu(graph.num_nodes());
        std::vector<std::vector<size_type> > &R = node.get_R();
        std::vector<double> &lambda = node.get_lambda();

        //The current 1_tree has the required edges
        node.last_1_tree = R;
        node.last_total_cost = 0 ;
        size_type joins = 0;

        //The sorted edges will be used to produce the 1-tree
        std::vector<ED :: Edge> &edges = graph.get_edges();

        //We update the nodes that are joined by the required edges
        for(size_type i =0 ; i < R.size(); i++) {
            for(size_type j =0 ; j < R[i].size();j++){
                ED :: Edge &e = graph.get_backup_edge(i, R[i][j]);
                if(i >= R[i][j])
                    continue;
                node.last_total_cost += e.get_dist();
                if(i && R[i][j]){
                    joins++;
                    dsu.join(i, R[i][j]);
                }
            }
        }

        for(size_type i = 0 ; i < edges.size() && edges[i].get_dist() < infinity && joins + 1 < graph.num_nodes(); i++){
            size_type u = edges[i].get_first()->get_id();
            size_type v = edges[i].get_second()->get_id();
            //skip if the edge contains the first node as endpoint
            if(!u || !v)
                continue;
    
            //We join the nodes in the edge and add the edge to the 
            //1-tree
            if(dsu.find(u)!=dsu.find(v)) {
                joins++;
                dsu.join(u,v);
                node.last_total_cost += edges[i].get_dist();
                node.last_1_tree[u].push_back(v);
                node.last_1_tree[v].push_back(u);
            }
        }
    
        //We look for the best two edges with least cost
        //to add the first node
        int cnt = node.last_1_tree[0].size();
        for(size_type i =0 ; i < edges.size() && cnt < 2; i++){
            size_type u = edges[i].get_first()->get_id();
            size_type v = edges[i].get_second()->get_id();
            if(std::count(R[0].begin(), R[0].end(), u) || std::count(R[0].begin(), R[0].end(), v))
                continue;

            if(!u || !v) {
                cnt++;
                node.last_1_tree[u].push_back(v);
                node.last_1_tree[v].push_back(u);

                node.last_total_cost += edges[i].get_dist();
            }
        }
	    assert(node.last_1_tree[0].size() == 2);
	    node.actual_cost = node.last_total_cost;
        for(size_type i =0 ; i < lambda.size();i++){
            node.last_total_cost += lambda[i] * (node.last_1_tree[i].size() - 2.0);
        }

        if(!node.last_second_1_tree.size())
            node.last_second_1_tree = node.last_1_tree;

        node.num_joins = joins;
    }

    void HeldKarp :: set_upper_bound() {

	U =0.0;
	best_solution = std::vector<std::vector<size_type> >(graph.num_nodes());
	for(size_type i =0 ; i < graph.num_nodes(); i++){
	    U += graph.get_distance(i, (i+1) % graph.num_nodes());
	    best_solution[i].push_back((i+1) % graph.num_nodes());
	    best_solution[(1+i) % graph.num_nodes()].push_back(i);
	}

    }

    void HeldKarp :: branch_and_bound(){

        //We generate here n^2 edges
	//std::cout << "BLA1\n";
        graph.generate_edges();
	    set_upper_bound();
        //Best bound heuristic
        std::priority_queue<SearchNode> q;
        //std::stack<SearchNode> q;
        SearchNode root = SearchNode(graph.num_nodes());
        root.root_node = true;
        root.invalid_node = false;
	    root.last_1_tree.clear();

        run_HK(root); 
        q.push(root);
        std::cout << U << "=U\n";

        //While there is a node in the search space
        while(!q.empty()){
            SearchNode node = q.top();
            q.pop();
            //We run HeldKarp algorithm

            //If the tree found has a worse cost than the best
            //found so far then we discard the node
            if(node.last_total_cost >= U || node.invalid_node)
                continue;

            
            //If the tree found is 2-regular, then we
            //update our solution
            if(node.check_1_tree() && node.solution_is_2_regular()){
                U = node.last_total_cost;
                std::cout << U << "=U\n";
                best_solution = node.last_1_tree;
            }
            //Otherwise we look for the children
            else {
                std::vector<SearchNode> children = node.get_children();
                for(size_type i=0 ; i < children.size(); i++){
                    run_HK(children[i]);
                    if(children[i].last_total_cost < 1e9){
                        q.push(children[i]);
                    }
                }
            }
        }

        //End of the search space, we print the best solution
        output_ans();
    
    }

};

