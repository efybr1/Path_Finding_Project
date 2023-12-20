//Ben Richards

//---------------------------------Description-------------------------------------------//
//Testing Dijkstra's naive implementation of Dijkstra's

//----------------------------------Includes---------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <utility>
#include <algorithm>  // Include for std::sort, std::unique
#include <math.h>

#include "node.h"

auto shortestPathComparator = [](const Node* a, const Node* b) -> bool
{
    return a->shortestPathToNode < b->shortestPathToNode;
};

int main() {
    //------------------------------------Open files-------------------------------------//
    std::ifstream inputEleFile("C:/A_Project/triangle/A.1.ele");
    std::ifstream NodeInputFile("C:/A_Project/triangle/A.1.node");
    //std::ifstream inputEleFile("C:/A_Project/Testing/File Manipulation/Star.1.ele");
    //std::ifstream NodeInputFile("C:/A_Project/Testing/File Manipulation/Star.1.node");
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

    //-----------------------------Read initial data-------------------------------------//
    // Read the first line of the .ele file and display node stats
    unsigned int numTriangles, nodesPerTriangle, numAttributes;
    inputEleFile >> numTriangles >> nodesPerTriangle >> numAttributes;
    std::cout << "Ele file stats: " << std::endl;
    std::cout << "Number of triangles: " << numTriangles << "\nNodes per triangle: " << nodesPerTriangle << "\nNumber of attributes: " << numAttributes << "\n" << std::endl;

    // Read the first line of the .node file and display node stats
    unsigned int numVertices, dimension, numNodeAttributes, numBoundaryMarkers;
    NodeInputFile >> numVertices >> dimension >> numNodeAttributes >> numBoundaryMarkers;
    std::cout << "Node file stats: " << std::endl;
    std::cout << "Number of vertices: " << numVertices << "\nDimension: " << dimension << "\nNumber of attributes: " << numAttributes << "\nNumber of boundary markers: " << numBoundaryMarkers << "\n"<< std::endl;

    //----------------------------Read Nodes from .node file-----------------------------//
    std::vector<Node> nodes;
    std::string line;
    std::getline(NodeInputFile, line); // Skips the first line

    while (std::getline(NodeInputFile, line)) //use getline to iterate through all lines in the file
    {
        //Read the line
        std::istringstream iss(line);
        unsigned int index;
        double x, y;

        //Read data from line into variables
        iss >> index >> x >> y;

        // Create a Node object with constructor to populate data and add it to the vector
        Node newNode(index, x, y);
        nodes.push_back(newNode);
    }
    nodes.pop_back(); //Remove last line as it is read but does not contain node information (file format)

    //------------------------Read Edges from .ele file as pairs-------------------------//
    //Create a vector of pairs of unsigned ints
    std::vector<std::pair<unsigned int, unsigned int> > nodePairs;

    for (unsigned int i = 0; i < numTriangles; i++) //For all the triangles
    {
        // Read triangle information
        unsigned int triangleNumber, node1, node2, node3;
        inputEleFile >> triangleNumber >> node1 >> node2 >> node3;

        //Read the three edges into pairs
        std::pair<unsigned int, unsigned int> pair1(std::min(node1, node2), std::max(node1, node2));
        std::pair<unsigned int, unsigned int> pair2(std::min(node2, node3), std::max(node2, node3));
        std::pair<unsigned int, unsigned int> pair3(std::min(node3, node1), std::max(node3, node1));

        // Add pairs to the vector
        nodePairs.push_back(pair1);
        nodePairs.push_back(pair2);
        nodePairs.push_back(pair3);
    }

    //----------------------Sort Edges and Remove non-unique edges-----------------------//
    std::sort(nodePairs.begin(), nodePairs.end());
    // Remove non-unique pairs using std::unique
    //std::vector<std::pair<unsigned int, unsigned int> >::iterator is replaced by auto in later c++ versions
    auto uniqueEnd = std::unique(nodePairs.begin(), nodePairs.end());
    nodePairs.erase(uniqueEnd, nodePairs.end());

    //-----------------Add the adjacencies to the nodes using the pairs------------------//

    //As the nodes are in order by number in the nodes array, we can use this to be able to add the adjacencies
    //Although vectors start indexing from zero and we have 1,2,3... so using number-1

    unsigned int number = 1;
    for (unsigned int i = 0; i < nodePairs.size(); i++)
    {

        if(nodePairs[i].first==number)
        {
            nodes[nodePairs[i].first-1].setConnection(&nodes[nodePairs[i].second-1]);
            nodes[nodePairs[i].second-1].setConnection(&nodes[nodePairs[i].first-1]);
        }
        else
        {
            number++;
            nodes[nodePairs[i].first-1].setConnection(&nodes[nodePairs[i].second-1]);
            nodes[nodePairs[i].second-1].setConnection(&nodes[nodePairs[i].first-1]);
        }
    }
    //----------------------------Calculate Node Path Lengths----------------------------//

    for(unsigned int i=0;i<nodes.size();i++)
    {
        nodes[i].calculateEdges();
        /*std::cout << "Node: " << i+1 << std::endl;
        nodes[i].printNode();
        nodes[i].printConnectedNodes();
        nodes[i].printDistancesToConnections();
        printf("\n");*/
    }

    //------------------------------------Dijkstra's-------------------------------------//

    //Start and End node
    int startNodeNumber = 1;
    int endNodeNumber = 6;

    //Get pointers to all the nodes in the unvisitedNodes array
    std::vector<Node *> unvisitedNodes;
    for(size_t i=0;i<nodes.size();i++)
    {
        unvisitedNodes.push_back(&nodes[i]);
    }

    //Set the current node to the Start node and set it's path length to 0
    Node* currentNode;
    currentNode = &nodes[startNodeNumber-1];
    currentNode->setShortestPathToNode(0);

    while(!unvisitedNodes.empty())
    {
        auto minNodeIt = std::min_element(unvisitedNodes.begin(), unvisitedNodes.end(), shortestPathComparator);
        currentNode = *minNodeIt;
        currentNode->setVisited();
        currentNode->updateConnectedNodeLengths();
        unvisitedNodes.erase(minNodeIt);
    }

    // Find and print the path to end node
    std::cout << "\nShortest path to Node "<< endNodeNumber <<":\n";
    currentNode = &nodes[endNodeNumber - 1];
    while (currentNode != nullptr) {
        currentNode->printNode();
        currentNode = currentNode->prevNode;
    }

    return 0;
}
