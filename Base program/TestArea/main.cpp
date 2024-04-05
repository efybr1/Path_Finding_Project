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

class Node {
private:

    unsigned int number; //Node index
    double x,y, shortestPathToNode; //x, y position, also shortest path distance
    bool visited; //bool to make sure node isn't visited twice, this should be taken out for optimization
    Node* prevNode; //Previous node variable to be able to trace the shortest path back for output.

    std::vector<Node*> connectedNodes; //Connected node pointer array
    std::vector<double> connectedNodeLengths; //Connected node lengths

    double pythagoras(const Node* connectedNode) const {
        return std::sqrt(std::pow(connectedNode->x - x, 2) + std::pow(connectedNode->y - y, 2));
    }

public:

    //Constructor
    Node(unsigned int num, double x_coord, double y_coord): number(num), x(x_coord), y(y_coord), shortestPathToNode(std::numeric_limits<double>::max()), visited(false), prevNode(nullptr) {};
    // Copy constructor assignment operator overload
    Node& operator=(const Node& existingNode){
        if (this != &existingNode) {  // If trying to = itself, don't.
            number = existingNode.number;
            x = existingNode.x;
            y = existingNode.y;
            shortestPathToNode = existingNode.shortestPathToNode;
            visited = existingNode.visited;
            prevNode = existingNode.prevNode;

            // Copy connectedNodes vector
            connectedNodes = existingNode.connectedNodes;

            // Copy connectedNodeLengths vector
            connectedNodeLengths = existingNode.connectedNodeLengths;
        }
        return *this;
    }

    //Sets
    void setConnection(Node* connectedNode){
        connectedNodes.push_back(connectedNode);
    }
    void setShortestPathToNode(double length){
        shortestPathToNode = length;
    }
    void setVisited(){
        visited = true;
    }
    //Gets
    unsigned int getNumber(){
        return number;
    }
    double getShortestPathToNode(){
        return shortestPathToNode;
    }
    Node* getPrevNode(){
    return prevNode;
    }

    //Utility Functions
    void calculateEdges() {
        for (size_t i = 0; i < connectedNodes.size(); i++) {
            double distance = pythagoras(connectedNodes[i]);
            connectedNodeLengths.push_back(distance);
        }
    }
     void updateConnectedNodeLengths() {
        for (size_t i = 0; i < connectedNodes.size(); ++i) {
            Node* connectedNode = connectedNodes[i];
            double newLength = shortestPathToNode + connectedNodeLengths[i]; //Calculate length to connected node through this node.

            // Only update if the new length is smaller than the newly calculated value - to maintain the shortest path
            if (newLength < connectedNode->shortestPathToNode) {
                connectedNode->shortestPathToNode = newLength;
                connectedNode->prevNode = this;
            }
        }
    }


    //Prints
    void printNode() {
    std::cout << "Node number: " << number << ", x:  " << x << ", y: " << y << ", shortest path length: " << shortestPathToNode << ", visited:  " << visited << std::endl;
}
    void printConnectedNodes(){
    if (connectedNodes.size() > 0) {
        for (size_t i = 0; i < connectedNodes.size(); i++) {
            std::cout << "Connected Node: " << connectedNodes[i]->number << std::endl;
        }
    }
}
    void printDistancesToConnections(){
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        double distance = connectedNodeLengths[i];
        std::cout << "Distance from node " << number << " to node " << connectedNodes[i]->number << ": " << distance << std::endl;
    }
}
};




//----------------------------------Functions--------------------------------------------//

//---------------------------------------------------------------------------------------//
//----------------------------shortestPathComparator-------------------------------------//
//---------------------------------------------------------------------------------------//

// Function used to compare two nodes. Returns boolean of whether a is < b
auto shortestPathComparator = [](Node* a, Node* b) -> bool
{
    return a->getShortestPathToNode() < b->getShortestPathToNode();
};

//---------------------------------------------------------------------------------------//
//-----------------------------------readNodes-------------------------------------------//
//---------------------------------------------------------------------------------------//

// Function to read .node and .ele file to construct weighted graph of nodes.
// The .node file contains the position of each node, the .ele file the connectivity between them
void readNodes(std::vector<Node>& nodes)
{
    //Open the .node file for vertices
    std::ifstream NodeInputFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/08.node");
    if (!NodeInputFile.is_open())
    {
        std::cerr << "Error opening node file!" << std::endl;
        exit(1);
    }

    //Open the .ele file for segments
    std::ifstream inputEleFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/08.ele");
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

    std::clock_t start = std::clock();
    while(!unvisitedNodes.empty())
    {
        auto minNodeIt = std::min_element(unvisitedNodes.begin(), unvisitedNodes.end(), shortestPathComparator);
        currentNode = *minNodeIt;
        currentNode->setVisited();
        currentNode->updateConnectedNodeLengths();
        unvisitedNodes.erase(minNodeIt);
    }
    std::clock_t end = std::clock();

    // Find and print the path to end node

    std::cout << "\nShortest path to Node "<< endNodeNumber <<":\n";
    currentNode = &nodes[endNodeNumber - 1];
    while (currentNode != nullptr) {
    currentNode->printNode();
        currentNode = currentNode->getPrevNode();
    }

    std::cout << "Shortest path length: " << nodes[endNodeNumber - 1].getShortestPathToNode() << std::endl;
    std::cout << "Time: " << (double)(((end - start) / (double)CLOCKS_PER_SEC)) << "s" << '\n';

    std::ofstream outputFile("testout.txt");
    outputFile << (double)(((end - start) / (double)CLOCKS_PER_SEC)) << std::endl;
    outputFile.close();

    /*for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i].printNode();
    }*/

    return 0;
}
