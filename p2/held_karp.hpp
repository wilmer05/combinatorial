#ifndef ALG
#define ALG

#include<iostream>
#include<algorithm>
#include<utility>
#include"graph.hpp"


namespace ALGORITHM{ //Start of namespace ALGORITHM

    /**
        @class HeldKarp
        @brief This class implements the algoritm

    **/
    class HeldKarp{
        private:
            ED :: Graph graph;
            FILE *out;
        public:
            HeldKarp() : graph(ED :: Graph(0)){}
            HeldKarp(ED :: Graph g, FILE *output) : graph(g), out(output) {
            } 

            void run() {
               fprintf(out, "TYPE : TOUR\n");  
               fprintf(out, "DIMENSION : %lu\n", graph.num_nodes());  
               fprintf(out, "TOUR_SECTION\n");  
            }
    };


};

#endif //held_karp.hpp
