#ifndef ALG_HPP
#define ALG_HPP

#include "graph.hpp"
#include "dsu.hpp" 
#include <vector>
namespace ALG{

    using size_type = std::size_t;

    class Edmonds{
    
        public:
            Edmonds(size_type);

           void run(){
           
           }
//            ED::Graph extract_solution();


        private:
            ED::Graph graph;
            DSU::Dsu dsu;
            std::vector<bool> matched_nodes; 

    }; 


};


#endif
