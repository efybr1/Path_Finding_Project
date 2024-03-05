// Ben Richards

//---------------------------------Description-------------------------------------------//
// Testing Dijkstra's naive implementation of Dijkstra's

//----------------------------------Includes---------------------------------------------//

#include <fstream>
#include <sstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>   // Include for std::sort, std::unique
#include <ctime>
#include <queue>

#include "wrapper.h"

//Make copy and assignment overloads for wrapper
//Same as normal but simply swap node pointers
//Then the Node's need to update the wrapper they point to as well!

/*// Overloaded specialization of std::swap for Wrapper class. Calls member swap function of Wrapper class
namespace std {
    template<>void swap(Wrapper& a, Wrapper& b) {
        a.swap(b);
        std::cout << "Swap" << std::endl;
        //Update Node with what wrapper now holds it
    }
}*/

//----------------------------------Functions--------------------------------------------//

//---------------------------------------------------------------------------------------//
//----------------------------shortestPathComparators-------------------------------------//
//---------------------------------------------------------------------------------------//

// Function used to compare two nodes. Returns boolean of whether a is < b
auto shortestPathComparator(Node* a, Node* b) {
    return a->getShortestPathToNode() > b->getShortestPathToNode();
};

// Define lambda function for comparing wrappers
auto compareByShortestPath = [](const Wrapper &a, const Wrapper &b) {
    return shortestPathComparator(a.getNodePtr(), b.getNodePtr());
};

//---------------------------------------------------------------------------------------//
//-----------------------------------readNodes-------------------------------------------//
//---------------------------------------------------------------------------------------//

// Function to read .node and .ele file to construct weighted graph of nodes.
// The .node file contains the position of each node, the .ele file the connectivity between them
void readNodes(std::vector<Node>& nodes)
{
    //Open the .node file for vertices
    std::ifstream NodeInputFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/01.1.node");
    if (!NodeInputFile.is_open())
    {
        std::cerr << "Error opening node file!" << std::endl;
        exit(1);
    }

    //Open the .ele file for segments
    std::ifstream inputEleFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/01.1.ele");
    if (!inputEleFile.is_open())
    {
        std::cout << "Error opening the ele file." << std::endl;
        exit(1);
    }

    //Read first line of .node file for parameters
    unsigned int numVertices, dimension, numNodeAttributes, numBoundaryMarkers;
    NodeInputFile >> numVertices >> dimension >> numNodeAttributes >> numBoundaryMarkers;

    std::cout << "Node file stats: " << std::endl;
    std::cout << "Number of vertices: " << numVertices << "\nDimension: " << dimension << "\nNumber of attributes: " << numNodeAttributes << "\nNumber of boundary markers: " << numBoundaryMarkers << "\n" << std::endl;

    //Read rest of node file into node objects
    std::string line;
    std::getline(NodeInputFile, line); // Skips the first line

    while (std::getline(NodeInputFile, line))
    {
        std::istringstream iss(line);
        unsigned int index;
        double x, y;
        iss >> index >> x >> y;
        Node newNode(index, x, y);
        nodes.push_back(newNode);
    }
    nodes.pop_back(); // Remove last line as it is read but does not contain node information (file format has text on the last line)

    //Read first line of .ele file for segment/edge parameters
    unsigned int numTriangles, nodesPerTriangle, numAttributes;
    inputEleFile >> numTriangles >> nodesPerTriangle >> numAttributes;
    std::cout << "Ele file stats: " << std::endl;
    std::cout << "Number of triangles: " << numTriangles << "\nNodes per triangle: " << nodesPerTriangle << "\nNumber of attributes: " << numAttributes << "\n" << std::endl;


    //Read the segments/edges from the .ele file into STD::pair objects. For efficiency use node number not node objects
    std::vector<std::pair<unsigned int, unsigned int>> nodePairs;
    for (unsigned int i = 0; i < numTriangles; i++)
    {
        unsigned int triangleNumber, node1, node2, node3;
        inputEleFile >> triangleNumber >> node1 >> node2 >> node3;

        std::pair<unsigned int, unsigned int> pair1(std::min(node1, node2), std::max(node1, node2));
        std::pair<unsigned int, unsigned int> pair2(std::min(node2, node3), std::max(node2, node3));
        std::pair<unsigned int, unsigned int> pair3(std::min(node3, node1), std::max(node3, node1));

        nodePairs.push_back(pair1);
        nodePairs.push_back(pair2);
        nodePairs.push_back(pair3);
    }

    std::sort(nodePairs.begin(), nodePairs.end());
    auto uniqueEnd = std::unique(nodePairs.begin(), nodePairs.end());
    nodePairs.erase(uniqueEnd, nodePairs.end());

    unsigned int number = 1;
    for (unsigned int i = 0; i < nodePairs.size(); i++)
    {
        if (nodePairs[i].first == number)
        {
            nodes[nodePairs[i].first - 1].setConnection(&nodes[nodePairs[i].second - 1]);
            nodes[nodePairs[i].second - 1].setConnection(&nodes[nodePairs[i].first - 1]);
        }
        else
        {
            number++;
            nodes[nodePairs[i].first - 1].setConnection(&nodes[nodePairs[i].second - 1]);
            nodes[nodePairs[i].second - 1].setConnection(&nodes[nodePairs[i].first - 1]);
        }
    }

    for (unsigned int i = 0; i < nodes.size(); i++)
    {
        nodes[i].calculateEdges();
    }
}



int main()
{

    //Create vector to hold all node objects and import nodes from mesh
    std::vector<Node> nodes;
    readNodes(nodes);

    //------------------------------------Dijkstra's-------------------------------------//

    //Start and End node
    int startNodeNumber = 1;
    int endNodeNumber = 17;
    nodes[startNodeNumber-1].setShortestPathToNode(0);

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i].printNode();
    }

    // Create a vector of Wrapper objects, one for each node
    std::vector<Wrapper> pq;
    for (unsigned int i=0; i<nodes.size();i++)
    {
        pq.push_back(Wrapper(i+1, &nodes[i]));
    }

    // Convert vector into a priority queue using make_heap
    std::make_heap(pq.begin(), pq.end(), compareByShortestPath);

    nodes[17].setShortestPathToNode(-1);
    //The +4 here would really be found by the number of the wrapper of that node, here we know it is pos 3
    //std::push_heap(pq.begin(), std::distance(pq.begin()), compareByShortestPath); //The +4 is to step the iterator to position 3

    // Display information about nodes and their associated wrappers through the wrappers
    std::cout << "\n--------------------------------------After--------------------------------------" << std::endl;
    for (unsigned int i = 0; i < pq.size(); ++i)
    {
        std::cout<< "WN through wrapper: " << pq[i].getNumber() << ". Node's WN through Node pointer: " << pq[i].getNodePtr()->getWrapper()->getNumber() << ". Node number: " << pq[i].getNodePtr()->getNumber() << std::endl;
    }

    std::cout << "\n-----------------------------Popping top element------------------------------\n" << std::endl;


    std::pop_heap(pq.begin(), pq.end(), compareByShortestPath);
    Wrapper topElement = pq.back();
    topElement.printWrapper();
    topElement.getNodePtr()->printNode();
    topElement.getNodePtr()->printConnectedNodes();



    //topElement.getNodePtr()->updateConnectedNodeLengths();
    pq.pop_back();

    /*for(unsigned int i=0;i<5;i++)
    {
        std::pop_heap(pq.begin(), pq.end(), compareByShortestPath);
        Wrapper topElement = pq.back();
        std::cout << "Popped smallest wrapper information: ";
        topElement.printWrapper();
        std::cout << "Popped smallest node information: ";
        topElement.getNodePtr()->printNode();
        std::cout << std::endl;
        pq.pop_back();
    }*/
    return 0;
}
