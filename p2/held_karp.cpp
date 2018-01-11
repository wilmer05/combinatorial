
#include<iostream>
#include<algorithm>
#include<utility>
#include<queue>
#include<cmath>
#include"graph.hpp"
#include"dsu.hpp"
#include"held_karp.hpp"
namespace ALGORITHM{ //Start of namespace ALGORITHM

    bool SearchNode :: solution_is_2_regular(){
    
        bool regular_2 = true;
        for(size_type i = 0 ; regular_2 && i < last_1_tree.size(); i++){
            regular_2 = regular_2 && last_1_tree[i].size() == 2;
        } 

        return regular_2;
    }

    std::vector<SearchNode> SearchNode::get_children(){

        std::vector<std::vector<size_type> > adj_tree(lambda.size());
        std::vector<SearchNode> children;


        //Tmp variable to generate every child
        SearchNode child;
        child.R = R;
        child.F = F;
        child.lambda = lambda;

        for(size_type i =0 ; i < lambda.size(); i++){
            if(last_1_tree[i].size() > 2){
                for(size_type j =0 ; j < last_1_tree[i].size(); j++){

                    size_type adj_node = last_1_tree[i][j];
                    
                    //We forbid the edge from i to adj_node
                    if( i < adj_node){
                        child.F.push_back(std::make_pair(i, adj_node));
                        //Add the child to the list of children
                        children.push_back(child);
                        child.F.pop_back();
                    }

                    //If we can add a required edge incident to i and
                    //the adj_node
                    if(R[i].size() < 2 && child.R[adj_node].size() < 2){

                        //We required the edge
                        child.R[i].push_back(adj_node);
                        child.R[adj_node].push_back(i);

                        for(size_type k = 0; k < last_1_tree[i].size(); k++){
                            //We look for other one incident to i to
                            //forbid
                            if(j == k) continue;
                            size_type forb_node = last_1_tree[i][k];

                            //Forbid the edge
                            child.F.push_back(std::make_pair(i,forb_node));
                            //then we add the child to the list
                            children.push_back(child);
                            child.F.pop_back();

                            //Then we see if we can required both edges
                            //treated before in other child
                            if(!R[i].size() && child.R[forb_node].size() < 2 && k > j){
                                child.R[i].push_back(forb_node);
                                child.R[forb_node].push_back(i);
                                children.push_back(child);

                                child.R[i].pop_back();
                                child.R[forb_node].pop_back();
                                
                            }
                        }
                        child.R[i].pop_back();
                        child.R[adj_node].pop_back();
                    }
                } 
            }        
        }
        

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
        for(size_type i =0 ; i < lambda.size(); i++){
             lambda[i] += tstep * d_parameter * (1.0 * node.last_1_tree[i].size() - 2.0);

             //In the first step there is not a Tree from the 
             //-1 iteration
             if(step) lambda[i] += tstep * ((1 - d_parameter) * 1.0 * node.last_second_1_tree[i].size() - 2.0);
        } 
    
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
            last = current;
            if(best_solution[current][0] == last)
                current = best_solution[current][1];
            else current = best_solution[current][0];
       }
        

       fprintf(out,"-1\nEOF");
    }

    void HeldKarp :: compute_1_tree(SearchNode &node){
        node.last_second_1_tree = node.last_1_tree;
        DSU :: Dsu dsu(graph.num_nodes());
        std::vector<std::vector<size_type> > &R = node.get_R();

        //The current 1_tree has the required edges
        node.last_1_tree = R;
        node.last_total_cost = 0 ;
        size_type joins = 0;

        //We update the nodes that are joined by the required edges
        for(size_type i =0 ; i < R.size(); i++) {
            for(size_type j =0 ; j < R[i].size();j++){
                ED :: Edge &e = graph.get_backup_edge(i, R[i][j]);
                dsu.join(i, R[i][j]);
                std::cout << "Req " << i << " " << R[i][j] << "\n";
                if(i && R[i][j])
                    joins++;
                node.last_total_cost += e.get_dist();
            }
        }

        //The sorted edges will be used to produce the 1-tree
        std::vector<ED :: Edge> &edges = graph.get_edges();
        for(size_type i = 0 ; i < edges.size() && joins + 2 < graph.num_nodes(); i++){
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
            if(!u || !v) {
                cnt++;
                node.last_1_tree[u].push_back(v);
                node.last_1_tree[v].push_back(u);

                node.last_total_cost += edges[i].get_dist();
            }
        }
        assert(node.last_1_tree[0].size() == 2);
    }

    void HeldKarp :: branch_and_bound(){
        U = 1e18;

        //We generate here n^2 edges
        graph.generate_edges();

        //Best bound heuristic
        std::priority_queue<SearchNode> q;
        SearchNode root = SearchNode(graph.num_nodes());
        root.root_node = true;
        q.push(root);

        //While there is a node in the search space
        while(!q.empty()){
            SearchNode node = q.top();
            q.pop();

            //We run HeldKarp algorithm
            run_HK(node); 

            //If the tree found has a worse cost than the best
            //found so far then we discard the node
            if(node.last_total_cost >= U)
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
                    q.push(children[i]);
            }
        }

        //End of the search space, we print the best solution
        output_ans();
    
    }


};

