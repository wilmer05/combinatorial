
#include<iostream>
#include<algorithm>
#include<utility>
#include<queue>
#include<cmath>
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
        for(size_type i = 0 ; i < lambda.size(); i++){

            /*we check if we can put an edge as required
            for(size_type j = 0 ; j < last_1_tree[i].size(); j++)
                if(R[last_1_tree[i][j]].size() < 2)
                    cnt_possible++;
            */
            if(last_1_tree[i].size() > 2){

                if(R[i].size() < 2 && (max_counter_node > lambda.size()  || counter[i] > counter[max_counter_node]))
                    max_counter_node = i;

                if(R[i].size() == 1){
                    node_with_required = true;
                    node_p = i;
                    break;
                }
            }
        }

        if(!node_with_required){
            node_p = max_counter_node;
        }

       if(node_p > lambda.size())
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
            
        

        if(children.size() <2 || children.size() >=4){
            std::cout << node_p << "\n";
            print();
            std::cout << children.size() << "\n";
        }
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
        size_type N = (num_nodes + 3) / 4 + 5;
        if(node.is_root_node())
            N = (num_nodes * num_nodes + 49) / 50 + num_nodes + 15;
       
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
            compute_1_tree(node);
            if(node.solution_is_2_regular())
                break;
            if(node.num_joins != num_nodes - 2){
                node.invalid_node = true;
                return;
            }

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
            update_lambda_function(node, tstep, step);
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

    }

    void HeldKarp :: update_lambda_function(SearchNode &node, double tstep, size_type step) {

        //Update lambda function based on the nodes degrees and the t's
        std::vector<double> &lambda = node.get_lambda();
//        std::cout << tstep << "\n";
//        std::cout << "Lambda before\n";

//        for(size_type i =0 ; i < lambda.size(); i++)
//            std::cout << lambda[i] << "\n";
        for(size_type i =0 ; i < lambda.size(); i++){
             lambda[i] += tstep * d_parameter * (1.0 * node.last_1_tree[i].size() - 2.0);

             //In the first step there is not a Tree from the 
             //-1 iteration
             if(step) lambda[i] += tstep * ((1 - d_parameter) * 1.0 * node.last_second_1_tree[i].size() - 2.0);
        } 
    
//        for(size_type i =0 ; i < lambda.size(); i++)
//            std::cout << lambda[i] << "\n";
    }


    void HeldKarp :: output_ans() {
    
       fprintf(out, "TYPE : TOUR\n");  
       fprintf(out, "DIMENSION : %lu\n", graph.num_nodes());  
       fprintf(out, "TOUR_SECTION\n");  

       //The solution should be a 2-regular graph
       for(size_type i =0 ; i < best_solution.size(); i++){
            assert(best_solution[i].size() == 2);
       }

       
       assert(best_solution.size() == graph.num_nodes());
       size_type last,current;
       last = 0;
       current = best_solution[0][0];

       //We run through the 2-regular graph until we get
       //back to the first node
       fprintf(out, "%lu\n", 1UL);
       while(current!=0){
            fprintf(out,"%lu\n", current + 1);
            size_type tmp = last;
            last = current;
            if(best_solution[current][0] == tmp)
                current = best_solution[current][1];
            else current = best_solution[current][0];
       }
        

       fprintf(out,"-1\nEOF");
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
            size_type u = edges[i].get_first().get_id();
            size_type v = edges[i].get_second().get_id();
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
            size_type u = edges[i].get_first().get_id();
            size_type v = edges[i].get_second().get_id();
            if(std::count(R[0].begin(), R[0].end(), u) || std::count(R[0].begin(), R[0].end(), v))
                continue;

            if(!u || !v) {
                cnt++;
                node.last_1_tree[u].push_back(v);
                node.last_1_tree[v].push_back(u);

                node.last_total_cost += edges[i].get_dist();
            }
        }

        for(size_type i =0 ; i < lambda.size();i++){
            node.last_total_cost += lambda[i] * (node.last_1_tree[i].size() - 2.0);
        }
        node.num_joins = joins;
    }

    void HeldKarp :: branch_and_bound(){
        U = 1e18;

        //We generate here n^2 edges
        graph.generate_edges();

        //Best bound heuristic
        std::queue<SearchNode> q;
        SearchNode root = SearchNode(graph.num_nodes());
        root.root_node = true;
        root.invalid_node = false;
        q.push(root);

        //While there is a node in the search space
        while(!q.empty()){
            SearchNode node = q.front();
            q.pop();

            //We run HeldKarp algorithm
            run_HK(node); 

            //If the tree found has a worse cost than the best
            //found so far then we discard the node
            if(node.last_total_cost >= U || node.invalid_node)
                continue;

            
            //If the tree found is 2-regular, then we
            //update our solution
            if(node.solution_is_2_regular()){
                U = node.last_total_cost;
                best_solution = node.last_1_tree;
            }
            //Otherwise we look for the children
            else {
                std::vector<SearchNode> children = node.get_children();
                for(size_type i=0 ; i < children.size(); i++)
                    if(children[i].last_total_cost < 1e9){
                        q.push(children[i]);
                    }
            }
        }

        //End of the search space, we print the best solution
        output_ans();
    
    }

};

