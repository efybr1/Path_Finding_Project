#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *input_file, *output_file;
    char line[256];
    int num_lines, i;
    float x, y;

    // Open input node file
    input_file = fopen("C:/A_Project/A_TestMeshes/Test/Square.1.txt", "r");
    // Create and open output file
    output_file = fopen("C:/A_Project/A_TestMeshes/Test/Output.txt", "w");


    if (input_file == NULL) {
        printf("Error opening input file.\n");
        return 1;
    }


    if (output_file == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }

    // Get the number of nodes from the first line
    fgets(line, sizeof(line), input_file);
    sscanf(line, "%d", &num_lines);

    // loop that number of nodes
    for (i = 0; i < num_lines; i++) {
        // Reads line and splits it by the spaces
        fgets(line, sizeof(line), input_file);
        sscanf(line, "%*f %f %f", &x, &y);

        // print the output in the format (x,y) to output file
        fprintf(output_file, "(%f,%f)\n", x, y);
    }


    fclose(input_file);
    fclose(output_file);

    return 0;
}




