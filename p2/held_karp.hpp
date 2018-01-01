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

            /**
                @brief Runs the algorithm  to compute a tour
            **/
            void run();
    };


};

#endif //held_karp.hpp
