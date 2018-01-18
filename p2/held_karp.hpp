#ifndef ALG
#define ALG

#include<iostream>
#include<algorithm>
#include<utility>
#include"graph.hpp"


namespace ALGORITHM{ //Start of namespace ALGORITHM

    const double d_parameter = 0.6;
    const double infinity = 1e9;
    const double epsilon = 1e-6;
    using size_type = std::size_t;

    /**
        @class Search node
        @brief This class implements a search node for the branch 
        and bound

    **/

    class SearchNode{
        public:

            //
            // Constructors
            //
            SearchNode() = default;

            SearchNode(size_type num_nodes) : root_node(false), lambda(std::vector<double> ( num_nodes, 0.0)){}


            //Getters
            std::vector<std::vector<size_type> > &get_R(){
                return R;
            }

            std::vector<std::pair<size_type, size_type> > &get_F(){
                return F;
            }


            std::vector<double> &get_lambda(){
                return lambda;
            }

            //Setters
            void set_lambda(std::vector<double> l){
                lambda = l;
            }


            /**
                @brief Function that indicates if a node is root
            **/
            bool is_root_node(){
                return root_node;
            }

            /**
                @brief Checks that the graph is connected
                and computes if the last solution is
                2-regular or not
            **/
            bool solution_is_2_regular();

            /**
                @brief redefining operator to sort nodes in the 
                priority queue
            **/
            bool operator<(const SearchNode &other) const {
                return last_total_cost > other.last_total_cost;
            }

            /**
                @brief Function that computes the children
                of a search node, based on the required and forbidden
                edges, and adding to them an edge, or a pair of them
            **/

            /**
                @brief Prints info related to the node just for debugging
            **/
            void print();

            std::vector<SearchNode> get_children();

            //Last two 1-trees found with the required and forbidden
            //Edges of this node
            std::vector<std::vector<size_type> > last_1_tree, last_second_1_tree;


            /**
                @brief indicated if an edge to a given node 
                belongs to the required nodes
            **/
            bool can_require(std::vector<std::vector<size_type> >&R, size_type node_u, size_type node_v);

            bool not_required(std::vector<size_type> &R, size_type node_v);



            /**
                @brief Check if a solution found is actually a 1-tree
                just to make sure that there are not cycles on other 
                length than n
            **/
	        bool check_1_tree();

            //Cost of the last solution found
            double last_total_cost;

            //Actual cost of the last 1_tree (without including 
            //the lambdas)
	        double actual_cost;

            //Boolean that indicates if it's the first node
            //in the search space
            bool root_node;

            /**
                lists of required and forbidden edges
            **/
            std::vector<std::pair<size_type, size_type> > F;
            std::vector<std::vector<size_type> > R;

            /**
                @brief a node is not valid when it has a cycle
                with the required edges that is not a tour
            **/ 
            bool invalid_node;

            /**
                lambda computed for this node so far
            **/
            std::vector<double> lambda;
    };

    /**
        @class HeldKarp
        @brief This class implements the algoritm

    **/
    class HeldKarp{
        private:
            //Graph of the problem
            ED :: Graph graph;

            //File to write or stdout
            FILE *out;

            //Variable for the cost of the best solution found so far
            double U;

            //Best solution found so far
            std::vector<std::vector<size_type> > best_solution;

            //Lambda computed in the root
            std::vector<double> lambda_root;

            //Sum to initialize the t0 on nodes that are not the 
            //root
            double sum_lambda_root;

        public:
            HeldKarp() : graph(ED :: Graph(0)){}
            HeldKarp(ED :: Graph g, FILE *output) : graph(g), out(output) {
            } 

            void branch_and_bound();

            /**
                @brief Runs the HK  algorithm  to compute a tour
            **/
            void run_HK(SearchNode &node);

            /**
                @brief prints the result of the algorithm
            **/
            void output_ans();


            /**
                @brief get the cost of an edge between u and v
                including the lambdas of u and v
            **/
            double get_lambda_edge_cost(SearchNode &node, size_type u, size_type v);

            /**
                @brief computes the min weight 1-tree based on 
                the required and forbidden edges in the SearchNode
            **/
            void prims(SearchNode &node);

            /**
                @brief we update the lambda function of the current node based on the time step and the ti of the current step
            **/
            void update_lambda_function(SearchNode &node, double tstep);

            /**
                @brief function to compute a first upper bound
            **/
	        void set_upper_bound();

    };


};

#endif //held_karp.hpp
