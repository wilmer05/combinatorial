#ifndef DSU_HPP
#define DSU_HPP
#include<vector>

namespace DSU{

    using size_type = std::size_t;


    /**
       @class Dsu
       @brief A @c Dsu stores an array of identifiers to implement the disjoint set union
              datastructure
    **/

    class Dsu{
        public:
            
            Dsu() = default;

            Dsu(size_type const);

            void add_node(size_type);
            void add_node();
            void clean();

            size_type find(size_type);
            void join(size_type, size_type);
            size_type num_nodes() const;
            bool have_same_tree(size_type, size_type) ;

        private:
            std::vector<size_type> _parent; 
              

    };

inline void Dsu::clean() {

    _parent.clear();

}

inline size_type Dsu::num_nodes() const{

    return _parent.size();

} 

};



#endif
