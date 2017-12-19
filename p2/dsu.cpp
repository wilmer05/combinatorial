#include "dsu.hpp"

namespace DSU{


    ////////////////////////////////////////////
    //! \c Dsu definition
    ////////////////////////////////////////////

    Dsu::Dsu(size_type size) : _parent(size){
        original_size = size;
        for(size_type i = 0 ; i < size ; i++)
            _parent[i] = i;
    }

    void Dsu::add_node(size_type parent){
        _parent.push_back(parent); 
    }


    void Dsu::add_node(){
        int size = _parent.size();
        _parent.push_back(size); 
    }

    size_type Dsu::find(size_type node_id){
   
        return _parent[node_id] = (node_id == _parent[node_id] ? node_id : find(_parent[node_id])); 

    }

    void Dsu::join(size_type node_u, size_type node_v){
        size_type parent_u = find(node_u); 
        size_type parent_v = find(node_v); 
        
        _parent[parent_u] = parent_v;
    }

    bool Dsu::have_same_tree(size_type node_u, size_type node_v){
        return find(node_u) == find(node_v); 
    }

};
