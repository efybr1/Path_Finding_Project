//Ben Richards

//---------------------------------Description-------------------------------------------//
//Program to add node and ele file together for viewing in "veeMeshViewer" software.
//Requires change to first line of the node file.
//And also the auto-generated comment in the node file to be deleted for some reason.
//Node and ele file then written to a .tri file and this is the output of the program.

//----------------------------------Includes---------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

//------------------------------------Structs--------------------------------------------//
struct Triangle
{
    unsigned int number;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> attributes;
};
struct Vertex
{
    unsigned int index;
    double x, y;
    unsigned int attributes;
    unsigned int boundaryMarker;
};

//--------------------------------Print Functions----------------------------------------//
void printTriangle(Triangle& triangle,int numAttributes)
 {
    std::cout << "Number: " << triangle.number << ", Nodes: ";
    for (unsigned int j = 0; j < triangle.nodes.size(); ++j) {
        std::cout << triangle.nodes[j] << " ";
    }
    if(numAttributes>0)
    {
        std::cout << "\nAttributes: ";
        for (unsigned int k = 0; k < triangle.attributes.size(); ++k) {
        std::cout << triangle.attributes[k] << " ";
    }
    }
    std::cout << std::endl;
}

void printVertex(Vertex& vertex)
{
    std::cout << "Vertex Number: " << vertex.index << ": Coordinates: (" << vertex.x << ", " << vertex.y << ")" << std::endl;
}

int main() {

//----------------------------------Open files-------------------------------------------//
    std::ifstream inputFile("C:/A_Project/triangle/A.1.ele");
    std::ifstream NodeInputFile("C:/A_Project/triangle/A.1.node");

    // Check if the file is opened successfully
    if (!inputFile.is_open()) {
        std::cout << "Error opening the file." << std::endl;
        return 1;}
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;}

//-------------------------------Read initial data---------------------------------------//
    // Read the first line of ele file and display ele stats
    unsigned int numTriangles, nodesPerTriangle, numAttributes;
    inputFile >> numTriangles >> nodesPerTriangle >> numAttributes;
    std::cout << "Ele file stats: " << std::endl;
    std::cout << "Number of triangles: " << numTriangles << "\nNodes per triangle: " << nodesPerTriangle << "\nNumber of attributes: " << numAttributes << "\n" << std::endl;

    // Read the first line of node file and display node stats
    unsigned int numVertices, dimension, numNodeAttributes, numBoundaryMarkers;
    NodeInputFile >> numVertices >> dimension >> numNodeAttributes >> numBoundaryMarkers;
    std::cout << "Node file stats: " << std::endl;
    std::cout << "Number of vertices: " << numVertices << "\nDimension: " << dimension << "\nNumber of attributes: " << numAttributes << "\nNumber of boundary markers: " << numBoundaryMarkers << "\n"<< std::endl;

    //-----------------------------Read Triangle data------------------------------------//
    // Read the remaining lines of ele file into the triangles vector
    std::vector<Triangle> triangles;
    for (int i = 0; i < numTriangles; ++i)
    {
        Triangle triangle; //Create triangle
        inputFile >> triangle.number; //Read triangle number

        // Read nodes
        for (int j = 0; j < nodesPerTriangle; ++j) {
            unsigned int node;
            inputFile >> node;
            triangle.nodes.push_back(node);}

        // Read attributes
        for (int k = 0; k < numAttributes; ++k) {
            unsigned int attribute;
            inputFile >> attribute;
            triangle.attributes.push_back(attribute);}

        //Push the now populated triangle onto the triangles array
        triangles.push_back(triangle);
    }
    // Close the file
    inputFile.close();

    //-------------------------------Read Vertex data------------------------------------//
    std::vector<Vertex> vertices;
    std::string line;
    std::getline(NodeInputFile, line); //Iterates over first line

    while (std::getline(NodeInputFile, line))
    {
        std::istringstream iss(line);
        Vertex vertex;
        iss >> vertex.index >> vertex.x >> vertex.y >> vertex.attributes >> vertex.boundaryMarker;
        vertices.push_back(vertex);
    }

    NodeInputFile.close();


    //For triangle 0
    //printTriangle(triangles[0],numAttributes);
    std::cout<< "Triangle 1 Nodes: " << std::endl;
    for(int i=0;i < triangles[0].nodes.size();i++)
    {
        std::cout << "Node "<< triangles[0].nodes[i] << ": (" << vertices[triangles[0].nodes[i]-1].x << ", " <<vertices[triangles[0].nodes[i]-1].y << ")" <<"" << std::endl;
    }
    return 0;



}

