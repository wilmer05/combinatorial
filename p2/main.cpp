#include <iostream>
#include <cstdlib>
#include <cstring>
#include "graph.hpp"
#include "held_karp.hpp"

const int max_string_size = 50000;
char input[max_string_size];
int n_vertices;
const char DIMENSION_STR[] = "DIMENSION";
const char NODE_COORD_STR[] = "NODE_COORD_SECTION";
const char EOF_STR[] = "EOF";

int usage(){
        std::cout << "Usage: ./main --instance file.tsp [--solution file.opt.tour]\nThese files should exist.\n";
        exit(EXIT_FAILURE); 
}

void ignore_input_until(char c, char c2){
   int len = strlen(input);
   int last_idx = 0;
   for(int i = 0 ; i < len; i++){
       if(input[i] == c){
           if(c2 != 0 && i + 1 < len && input[1 + i] == c2)
                continue;
           else if(!c2)
                continue;
       }
       input[last_idx++] = input[i];
   }
   input[last_idx] = 0;

}

bool input_starts_with_str(const char *str, int len_input){

    int len_str = strlen(str);
    if(len_str > len_input) 
        return false;

    for(int i = 0; i < len_str; i++)
        if(input[i] != str[i])
            return false;
    
    return true;
}

void read_file(char *file_name, ED::Graph &graph){
    FILE * file;
    file =fopen(file_name, "r");
    if(!file)
        usage();
    int n_vertices, node;
    double coord_x, coord_y;
    
    //Reading file while it is not the end of file
    bool dimension_readed = false;
    bool node_coord_readed = false;
    while(fscanf(file, "%[^\n]\n", input) !=EOF){
        
        ignore_input_until(':', 0);
        ignore_input_until(' ', ' ');
        int len = strlen(input);
       
        if(input_starts_with_str(EOF_STR, len))
            break;
        //Getting the size of the graph
        if(node_coord_readed){
            sscanf(input, "%d %lf %lf", &node, &coord_x, &coord_y); 
            node--;
            graph.set_coordinates(node, coord_x, coord_y);
        }
        
        if(input_starts_with_str(DIMENSION_STR, len)) {
            sscanf(input, "DIMENSION %d", &n_vertices); 
            graph = ED::Graph(n_vertices);
            dimension_readed = true;
        }

        if(input_starts_with_str(NODE_COORD_STR, len))
            node_coord_readed = true;
         
    }
    fclose(file);
}

void solve(ED :: Graph &g, FILE *output){
    ALGORITHM :: HeldKarp algorithm(g, output);
    /*
        We run our algorithm and print the result
    */
    algorithm.run();
}

int main(int argc, char**argv)
{
    ED :: Graph graph(0);
    FILE *output = stdout;

    /*
        Parse of the input parameters, we check if there is some 
        hint or not, and if the number of parameters is correct
    */
    if(argc >= 3){
        if(argc == 4 || argc >= 6 || strcmp(argv[1], "--instance")) 
            return usage();
       
        if(argc > 3 && strcmp(argv[3], "--solution"))
            return usage();

        read_file(argv[2], graph);
        if(argc > 3 ){
            if(strcmp(argv[3], "--solution")) 
                return usage();
            output = fopen(argv[4], "w");
        }

        
        solve(graph, output); 
        fclose(output);
    } else
        return usage(); 

   return EXIT_SUCCESS;
}
