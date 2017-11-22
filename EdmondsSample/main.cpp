#include <iostream>
#include <cstdlib>
#include <cstring>
#include "graph.hpp"
#include "edmonds.hpp"

const int max_string_size = 50000;
char input[max_string_size];

//! return graph consisting of path with \c num_nodes many vertices
/*static ED::Graph create_path(ED::NodeId num_nodes)
{
   ED::Graph result{num_nodes};

   for (ED::NodeId node_id = 0; node_id + 1 < num_nodes; ++node_id)
   {
      result.add_edge(node_id, node_id + 1);
   }
   // will be moved
   return result;
}*/

int usage(){
        std::cout << "Usage: ./main --graph file1.mdx [--hint file2.dmx]\n";
        return EXIT_FAILURE; 
}

void read_file(char *file_name, ED::Graph &graph){
    FILE * file;
    file =fopen(file_name, "r");
    int m_edges;
    int n_vertices;
    int u,v;

    //Reading file while it is not the end of file

    while(fscanf(file, "%[^\n]\n", input) !=EOF){
        //Skip comments
        if(input[0] == 'c') continue;

        //Getting the size of the graph
        if(input[0] == 'p') {
            sscanf(input, "p edge %d %d", &n_vertices, &m_edges); 
            graph = ED::Graph(n_vertices);
        }


        //Reading edge
        if(input[0] == 'e'){
            sscanf(input, "e %d %d", &u, &v);      
            u--; v--;
            graph.add_edge(u,v);
        }
         
    }
    fclose(file);
}

void solve_edmonds(ED :: Graph &g, ED :: Graph &hint){
    ALG::Edmonds algorithm(g.num_nodes());

    for(ED::size_type node = 0 ; node < g.num_nodes(); node++){
        for(ED ::size_type neighbour : g.node( node ).neighbors()){
            algorithm.graph.add_edge(node, neighbour);
        } 
    }
    for(ED::size_type node = 0 ; node < hint.num_nodes(); node++){
        for(ED::size_type neighbour : hint.node( node ).neighbors()){
            algorithm.match(node, neighbour);
        } 
    }
    for(ED::size_type node = 0 ; node < g.num_nodes(); node++){
        for(ED ::size_type neighbour : g.node( node ).neighbors()){
            if(algorithm.exposed_vertex(node) && algorithm.exposed_vertex(neighbour)){
                algorithm.match(node, neighbour);
                break;
             }
        } 
    }
    algorithm.run();


      
}

int main(int argc, char**argv)
{
    ED :: Graph graph(0), initial_matching(0);
    if(argc >= 3){
        if(argc == 4 || argc >= 6 || strcmp(argv[1], "--graph")) 
            return usage();
       
        if(argc > 3 && strcmp(argv[3], "--hint"))
            return usage();

        read_file(argv[2], graph);
        if(argc > 3 ){
            if(strcmp(argv[3], "--hint")) 
                return usage();
            std::cout << "Reading hint\n";
            read_file(argv[4], initial_matching); 
        }

        solve_edmonds(graph, initial_matching); 

    } else
        return usage(); 

   return EXIT_SUCCESS;
}
