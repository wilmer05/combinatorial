#ifndef DSU_HPP
#define DSU_HPP
#include<vector>
#include <cassert>

namespace DSU{

    using size_type = std::size_t;


    /**
       @class Dsu
       @brief A @c Dsu stores an array of identifiers to implement the disjoint set union datastructure
    **/

    class Dsu{
        public:


            /** @brief Creates an empty tree **/
            Dsu() = default;

            /** @brief Creates a tree with a fixed number of nodes**/
            Dsu(size_type const size);

            /** 
                @brief adds a node to the current tree. 
                This function is used when we create new nodes 
                due to odd cycles 
            **/
            void add_node();

            /** 
                @brief adds a node to the current tree with a fiex 
                parent
            **/
            void add_node(size_type parent);

            /** 
                @brief Clean the current tree, leaving an 
                empty tree 
            **/
            void clean();

            /** @brief Clean the current tree, leaving the same 
            tree or the original size (before the added composed 
            nodes) and reset the parents of the nodes in the vector*/
            void clean(std::vector<size_type> const &nodes);

            /** 
                @brief Finds the root of the tree where node 
                belongs to. It means the representative
            **/
            size_type find(size_type node);


            /** 
                @brief Join the trees where the node_u and 
                node_v belong
            **/
            void join(size_type node_u, size_type node_v);

            /** 
                @brief Returns the number of nodes in the 
                current tree 
            **/
            size_type num_nodes() const;
            
            /** @brief Returns true if node_u and node_v are in the same tree**/ 
            bool have_same_tree(size_type node_u, size_type node_v);

            /** 
                @brief Sets the root of the tree where node_v 
                belongs as the parent of the root of the tree 
                where node_u belongs
            **/
            void set_parent(size_type node_u, size_type node_v);

        private:
            std::vector<size_type> _parent; 
            size_type original_size;
              

    }; // Class Dsu

inline void Dsu::clean() {

    _parent.clear();

}
inline void Dsu::clean(std::vector<size_type> const &nodes){
    for(size_type i = 0 ; i < nodes.size(); i++){
        _parent[nodes[i]] = nodes[i];
    }
    _parent.resize(original_size);
    for(size_type i =0 ; i < _parent.size(); i++) assert(_parent[i] == i);
}

inline size_type Dsu::num_nodes() const{

    return _parent.size();

} 


inline void Dsu::set_parent(size_type node_u, size_type node_v){
    size_type repr_v = find(node_v);
    size_type repr_u = find(node_u);

    _parent[repr_u] = repr_v;

} 

};



#endif
