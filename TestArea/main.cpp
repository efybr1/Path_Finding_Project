#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

//Class to encapsulate the concept of a node
class Node {
public:
    unsigned int number; //Node index
    double x, y; //x,y position
    std::vector<Node*> connectedNodes; //Connected node pointer array

    //Constructors
    Node(unsigned int num, double x_coord, double y_coord):number(num),x(x_coord),y(y_coord){} //Full constructor to initialise index and coordinates

    Node(unsigned int num):number(num){} //Constructor to initialise with only number

    //Sets
    void setConnection(Node* connectedNode)
    {
        connectedNodes.push_back(connectedNode);
    }
    void setCoordinates(double x_coord,double y_coord)
    {
        x=x_coord;
        y=y_coord;
    }

    //Prints
    void print()
    {
        std::cout << "Node number: " << number << ", x:  " << x << ", y: " << y << std::endl;
    }
    void printConnectedNodes()
    {
        if (connectedNodes.size()>0)
        {
            for(size_t i=0;i<connectedNodes.size();i++)
            {
                std::cout << "Connected Node: " <<connectedNodes[i]->number << std::endl;
            }
        }
    }
};

//Node needs to have a vector of the edges


int main() {

    std::ifstream inputEleFile("C:/A_Project/triangle/A.1.ele");
    std::ifstream NodeInputFile("C:/A_Project/triangle/A.1.node");
    if (!inputEleFile.is_open())
    {
        std::cout << "Error opening the ele file." << std::endl;
        return 1;
    }
    if (!NodeInputFile.is_open())
    {
        std::cerr << "Error opening node file!" << std::endl;
        return 1;
    }
    unsigned int numTriangles, nodesPerTriangle, numAttributes;
    inputEleFile >> numTriangles >> nodesPerTriangle >> numAttributes;
    std::cout << "Ele file stats: " << std::endl;
    std::cout << "Number of triangles: " << numTriangles << "\nNodes per triangle: " << nodesPerTriangle << "\nNumber of attributes: " << numAttributes << "\n" << std::endl;
    // Read the first line of node file and display node stats
    unsigned int numVertices, dimension, numNodeAttributes, numBoundaryMarkers;
    NodeInputFile >> numVertices >> dimension >> numNodeAttributes >> numBoundaryMarkers;
    std::cout << "Node file stats: " << std::endl;
    std::cout << "Number of vertices: " << numVertices << "\nDimension: " << dimension << "\nNumber of attributes: " << numAttributes << "\nNumber of boundary markers: " << numBoundaryMarkers << "\n"<< std::endl;


    std::vector<Node> nodes;
    std::string line;

    std::getline(NodeInputFile, line); // Skips the first line

    while (std::getline(NodeInputFile, line))
    {
        std::istringstream iss(line);
        unsigned int index;
        double x, y;
        // Add attributes and boundaryMarker if needed
        iss >> index >> x >> y;
        // Create a Node object and add it to the vector
        Node newNode(index, x, y);
        nodes.push_back(newNode);
    }

    for(int i=0;i<numVertices;i++)
    {
        nodes[i].print();
    }




    return 0;
}
