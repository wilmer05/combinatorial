#ifndef ALG
#define ALG

#include<iostream>
#include<algorithm>
#include<utility>
#include"graph.hpp"


namespace ALGORITHM{ //Start of namespace ALGORITHM

    const double d_parameter = 0.6;
    using size_type = std::size_t;

    /**
        @class Search node
        @brief This class implements a search node for the branch 
        and bound

    **/

    class SearchNode{
        public:
            SearchNode() = default;

            SearchNode(size_type num_nodes) : lambda(std::vector<double> ( num_nodes, 0.0)), root_node(false){}

            std::vector<std::vector<size_type> > &get_R(){
                return R;
            }

            std::vector<std::pair<size_type, size_type> > &get_F(){
                return F;
            }


            std::vector<double> &get_lambda(){
                return lambda;
            }

            void set_lambda(std::vector<double> l){
                lambda = l;
            }

            bool is_root_node(){
                return root_node;
            }

            void compute_children();

            bool solution_is_2_regular();


            bool operator<(const SearchNode &other) const {
                return last_total_cost > other.last_total_cost;
            }

            std::vector<SearchNode> get_children();

            std::vector<std::vector<size_type> > last_1_tree, last_second_1_tree;
            double last_total_cost;
            std::vector<SearchNode> children;
        private:
            std::vector<double> lambda;
            std::vector<std::pair<size_type, size_type> > F;
            std::vector<std::vector<size_type> > R;
            bool root_node;
    };

    /**
        @class HeldKarp
        @brief This class implements the algoritm

    **/
    class HeldKarp{
        private:
            ED :: Graph graph;
            FILE *out;
            double U;
            std::vector<std::vector<size_type> > best_solution;
            std::vector<double> lambda_root;
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

            void compute_1_tree(SearchNode &node);

            void update_lambda_function(SearchNode &node, double tstep, size_type step);

    };


};

#endif //held_karp.hpp
