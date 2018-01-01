
#include<iostream>
#include<algorithm>
#include<utility>
#include"graph.hpp"
#include"held_karp.hpp"
namespace ALGORITHM{ //Start of namespace ALGORITHM

    /**
        @class HeldKarp
        @brief This class implements the algoritm

    **/
    void HeldKarp :: run() {
       fprintf(out, "TYPE : TOUR\n");  
       fprintf(out, "DIMENSION : %lu\n", graph.num_nodes());  
       fprintf(out, "TOUR_SECTION\n");  
    }


};

