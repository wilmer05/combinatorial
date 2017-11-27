#include <iostream>
#include <cstdlib>
#include <cstring>
#include "graph.hpp"
#include "edmonds.hpp"

const int max_string_size = 50000;
char input[max_string_size];

ALG::Edmonds algorithm(0);

int usage(){
        std::cout << "Usage: ./main --graph file1.mdx [--hint file2.dmx]\nThese files should exist.\n";
        exit(EXIT_FAILURE); 
}

void read_file(char *file_name, ED::Graph &graph){
    FILE * file;
    file =fopen(file_name, "r");
    int m_edges;
    int n_vertices;

    //Variables for nodes to read
    int node_u,node_v;
    
    if(!file)
        usage();
    
    //Reading file while it is not the end of file

    while(fscanf(file, "%[^\n]\n", input) !=EOF){

        //Skip comments
        if(input[0] == 'c') 
            continue;

        //Getting the size of the graph
        if(input[0] == 'p') {
            sscanf(input, "p edge %d %d", &n_vertices, &m_edges); 
            graph = ED::Graph(n_vertices);
        }

        //Reading edge
        if(input[0] == 'e'){
            sscanf(input, "e %d %d", &node_u, &node_v);      

            //Nodes are 0-indexed
            node_u--; node_v--;
            graph.add_edge(node_u,node_v);
        }
         
    }
    fclose(file);
}

void solve_edmonds(ED :: Graph &g, ED :: Graph &hint){
    algorithm = ALG::Edmonds(g.num_nodes());
    algorithm.graph = g;


    /*
        We add the matching hint to our match in the algorithm
    */
    for(ALG::size_type node = 0 ; node < hint.num_nodes(); node++){
        for(ED::size_type neighbour : hint.node( node ).neighbors()){
            algorithm.match(node, neighbour);
        } 
    }


    /*
        We run our algorithm and print the result
    */
    algorithm.run();
    algorithm.print_matching();
      
}

int main(int argc, char**argv)
{
    ED :: Graph graph(0), initial_matching(0);

    /*
        Parse of the input parameters, we check if there is some 
        hint or not, and if the number of parameters is correct
    */
    if(argc >= 3){
        if(argc == 4 || argc >= 6 || strcmp(argv[1], "--graph")) 
            return usage();
       
        if(argc > 3 && strcmp(argv[3], "--hint"))
            return usage();

        read_file(argv[2], graph);
        if(argc > 3 ){
            if(strcmp(argv[3], "--hint")) 
                return usage();
            read_file(argv[4], initial_matching); 
        }

        solve_edmonds(graph, initial_matching); 

    } else
        return usage(); 

   return EXIT_SUCCESS;
}
