#ifndef ALG_HPP
#define ALG_HPP

#include "graph.hpp"
#include "dsu.hpp" 
#include <queue>

/**
    @file edmonds.hpp

    @brief This file implements the Edmond's algorithm
**/
namespace ALG{ // for the implementation of Edmonds Algorithm

    using size_type = std::size_t;


    /** 
        @class Edmonds

        @brief Implements the edmonds algorithm, it has every operation in different methods, like shrinking, augment and extension
    **/
    class Edmonds{
    
        public:
           
           /** 
               @brief Constructor, it initializes every vector 
               inside the class
           **/
           Edmonds(size_type);

           ~Edmonds();

            /**
                @brief This function resets every vector inside the 
                class. It's used after an augmentation or after 
                a frustrated tree is found. Of course it doesn't 
                change the vector that says if a node is frustrated or
                not.
            **/
           void reset_vectors();

           /**
                @brief Function that finds an initial matching which
                is maximal
           **/
           void find_maximal_matching();

           /**
                @brief Given an exposed node, this functions sets 
                the node as the current root in the alternating tree 
                and starts the algorithm. This funtion calls the
                shrink, extend_tree and augment matching functions 
                depending on the case.

           **/
           void grow_tree(size_type root);


           /**
                @brief This is the main routine to run the algorithm, 
                it only looks for a exposed and non frustrated node 
                to call grow_tree in order to match the selected node
                if it's possible.
           **/
           void run();

           /**
                @brief Given two nodes, this function matches them
                and updates the vector matched_to to indicate that
                they are now matched.
           **/           
           void match(size_type, size_type);

           /**
                @brief Given two nodes that are connected by and edge,
                adds the edge to the pending list of edges to consider
                in order to grow the tree or find an exposed node.
                
           **/
           void add_edge_to_pending_list(size_type, size_type);

           /**
                @brief Given a node, adds the list of edges incident
                to it to the pending list of edges to consider
                to grow the alternating tree.
           **/
           void add_incident_edges_to_pending_list(size_type);

           /**
                @brief Given a node, returns if it is a pseudonode
                or not
           **/
           bool is_pseudonode(size_type);

           /**
                @brief Given a node, returns if it is an expsed vertex
                or not
           **/
           bool exposed_vertex(size_type);

           /**
                @brief Given two nodes node_u and node_v (node_v
                already belongs to the alternating tree) adds node_u 
                to the tree as well as the node that is matched to it.
           **/
           void extend_tree(size_type node_u, size_type node_v);

           /**
               @brief Given two nodes that are even (possibly 
               pseudonodes) in the alternating tree that close an odd
               cycle in the tree, it shrinks this odd cycle 
               by looking which is the common predecessor of them
               and then shrinking the cycle. It runs in 
               O(|E(C)|) where C is the cycle we are shrinking.
           **/
           void shrink(size_type, size_type);

           /**
                @brief Given an edge which is a pair of two nodes
                identifiers, and a node, returns if the node belongs
                to the edge
           **/
           bool node_belongs_to_edge(std::pair<size_type, size_type> edge, size_type node);

           /**
                @brief In order to find the common predecessor in the
                shrink step, sometimes we visit not contained
                edges or nodes in the cycle (at most |C|). So given
                a list of edges that could belong or not to the cycle
                it removes the edges that were not contained in the 
                cycle.
           **/
           void remove_not_contained_edges(std::vector<std::pair<size_type, size_type> > &edges);


            /**
                @brief Given a exposed node and it's parent,
                adds the node to the alternating tree and proceeds
                from the exposed vertex until the root matching every
                the nodes every two edges, and if one of the nodes
                in a new matching edge belongs to a pseudonode it
                calls the function match_all_pseudonodes_of in order
                to apply the proposition 22 to find a perfect matching
                in this pseudonode without the node matched in the
                highest level (in the current M-alternating tree).
                
            **/
            void augment_matching(size_type node, size_type parent);


            /**
                @brief Given a pseudonode and a node to ignore in 
                this pseudonode, applies the proposition 22
                (maybe recursively calling augment_matching_on_cycle 
                if there are some pseudonodes in the current 
                pseudonode) and matches the edges
                properly (the edges adjacent to the ignored node
                in the cycle of the pseudonode have to be unmatched).
            **/

            void augment_matching_on_cycle(size_type pseudonode, size_type ignored_node);


            /**
                @brief This function receives a node in the original 
                graph and a type of edge (matched for 1 and 0 
                otherwise) incident to the parent in the highest
                level of the alternating tree (the one with the
                biggest pseudonodes).
                
                If the incident edge to the parent has to be a 
                non-matched edge then this function matches the node
                and modifies the matching in every pseudonode where
                it belongs to, finding a perfect matching in those
                pseudonodes without this node. 
                
                Otherwise, the edge that has to be 
                mathed corresponds to the edge going from the
                node in the original graph connecting the pseudonode
                of the input node in the highest level of the 
                alternating tree to its parent.

                Note: The node to be matched in this step is given by
                the function get_next_node_to-match.
            **/
            void match_all_pseudonodes_of(size_type node, int next_kind_of_edge);


            /**
                @brief Given a node in the original graph and a type
                of incident edge to the parent of the pseudonode where
                this node belongs to, returns if we have to match
                the given node or the node that connectes the 
                pseudonode of this node to the parent.
            **/
            size_type get_next_node_to_match(size_type node, int incident_kind_of_edge);

            /**
                Prints the current matching of the algorithm
            **/
            void print_matching();

            /**
                @brief Function that adds the incident edges in the
                the last added pseudonode (cycle) to the list
                of pending edges to consider to grow the alternating
                tree.
            **/
            void add_incident_edges_of_last_cycle_odd_nodes();


            /**
                @brief Adds a new pseudonode to the tree and
                updates every vector
            **/
            void add_new_pseudonode(size_type predecessor, size_type pseudonode_number);

            /**
                @brief Join partition sets in a cycle (newpseudonode)
            **/
            void join_partition_sets(std::vector<std::pair<size_type, size_type> > &cycle, size_type pseudonode_id);

            
            /**
                @brief variable corresponding to the graph given
                as input
            **/
            ED::Graph graph;
        private:

            /**
                @brief Disjoint Set Union data structure to keep track
                of the current tree partitions
            **/
            DSU::Dsu dsu;

            /**
                @brief Queue (list) of the edges to be considered
                to grow the tree
            **/
            std::queue<std::pair<size_type, size_type> > pending_edges; 
            /**
                @brief Number of nodes in the original graph,
                (without the pseudonodes) 
            **/
            size_type num_original_nodes;

            /**
                @brief Current root of the alternating tree
            **/
            size_type current_root;

            /**
                @brief Vector that given a node tells if this
                node is already on the tree or not.
            **/
            std::vector<bool> on_tree;

            /**
                @brief Boolean vector that is true in a node
                is even in the current alternating tree or not
            **/
            std::vector<bool> even_node;

            /**
                @brief Boolean vector that is true in a node
                is odd in the current alternating tree or not
            **/
            std::vector<bool> odd_node;

            /**
                @brief Vector that indicates for a given node (or
                pseudonode) which is the parent of this node.

                Note: the parent is always a node in the original 
                graph. This helps to keep track to who's the actual
                nodes connecting two pseudonodes (or a pseudonode
                and a node) in the current alternating tree.
            **/
            std::vector<size_type> parent;

            /**
                @brief Vector that given a node or a pseudonode
                indicate which is the node in the original graph
                that connects this pseudonode to its parent, i.e
                for nodes in the original graph is the node itself
                but for pseudonodes is the node in the original graph
                connecting this pseudonode to the parent.
            **/
            std::vector<size_type> actual_node_to_parent;

            /**
                @brief Vector to keep track what is the node matched
                to a given node.
            **/
            std::vector<size_type> matched_to; 

            /**
                @brief Vector used in the shrinking step to keep track                which nodes are actually in the cycle to remove
                then edges that are not in the cycle.
            **/
            std::vector<bool> visited;

            /**
                @brief Vector of edges belonging two pseudonodes.
                These are two lists: 
                    edges_in_cycle are list of edges for each 
                    pseudonode whom incident vertices are nodes
                    in the original graph

                    peudonodes_edges is the list of edges for each
                    pseudonode whom incident vertices are nodes
                    in the alternating tree at a given moment
                    so the vertices in the edges thes lists 
                    could be pseudonodes or nodes in the original 
                    graph.
            **/
            std::vector<std::vector<std::pair< size_type, size_type> > > edges_in_cycle, pseudonodes_edges; 

            /**
                @brief Vector to keep track which nodes are currently
                in the tree used to clean everything in the 
                reset_vectors function.In this way we don't have 
                to reset every vector for all nodes in the graph
                but just for the ones we used during the creation
                of an alternating tree.
            **/
            std::vector<size_type> nodes_in_tree;

            /**
                @brief Boolean vector that indicates if a given node
                is frustrated or not
            **/
            std::vector<bool> frustrated;


            /**
                @brief Vector that indicates which is the smallest 
                pseudonode where given node or pseudonode belongs to.
            **/
            std::vector<size_type> belongs_to_pseudonode;

            /**
                @brief Boolean vector used in the augmentation 
                function to keep track which pseudonodes the 
                proposition 22 has been applied
            **/
            std::vector<bool> augmented;
    }; 


}; // namespace ALG


#endif
