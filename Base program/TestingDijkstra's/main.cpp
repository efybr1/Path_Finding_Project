// Ben Richards

//---------------------------------Description-------------------------------------------//
// Base program 'naive' implementation of Dijkstra's shortest path algorithm

//----------------------------------Includes---------------------------------------------//

#include <fstream> //Used for file inputs and outputs
#include <sstream>
#include <math.h> //Used for calculations of path lengths
#include <iostream> //Used from output to command window
#include <vector> //Used for data structures
#include <utility> //Utility functions
#include <algorithm>   // Include for std::sort, std::unique
#include <ctime> //Used for timing


// Class to encapsulate the concept of a Node for creating a weighted map for Dijkstra's to
// solve. The class contains the unique identifier, coordinates, the Node's current shortest
// path from the designated start node, the location of the Node before it in this shortest
// path, and all the associated standard and utility functions.

class Node {
private:
    unsigned int number; //Node index
    double x,y, shortestPathToNode; //x, y position, also shortest path distance
    bool visited; //bool to make sure node isn't visited twice, this should be taken out for optimization
    Node* prevNode; //Previous node variable to be able to trace the shortest path back for output.

    std::vector<Node*> connectedNodes; //Connected node pointer array
    std::vector<double> connectedNodeLengths; //Connected node lengths

    //Private function as this is not used outside of the class
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
    unsigned int getNumber() const{
        return number;
    }
    double getX() const{
        return x;
    }
    double getY() const{
        return y;
    }
    double getShortestPathToNode(){
        return shortestPathToNode;
    }
    Node* getPrevNode(){
    return prevNode;
    }
    const std::vector<Node*>& getConnectedNodes() const {
        return connectedNodes;
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
    std::ifstream NodeInputFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/MainProgramTestMeshes/04.node");
    if (!NodeInputFile.is_open())
    {
        std::cerr << "Error opening node file!" << std::endl;
        exit(1);
    }

    //Open the .ele file for segments
    std::ifstream inputEleFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/MainProgramTestMeshes/04.ele");
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

//Function to output a file showing the shortest path found by the program on the input mesh
void outputFile(std::vector<Node>& nodes, unsigned int endNodeNumberIn)
{
    //Open output file
    std::ofstream outputFile("Checking.plsg");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    //Output nodes
    //First line
    outputFile << nodes.size() << " " << 2 << " " << 0 << std::endl;
    //All nodes
    for (const auto& node : nodes)
    {
        outputFile << node.getNumber() << " " << node.getX() << " " << node.getY() << std::endl;
    }

    //Output edges/segments
    //First line
    std::vector<std::pair<unsigned int, unsigned int>> segments;

    //Find all segments
    for (Node node : nodes)
    {
        unsigned int currentNodeNumber = node.getNumber();
        std::vector<Node*> connectedNodes = node.getConnectedNodes();

        //Iterate through connected nodes
        for (Node* connectedNode : connectedNodes)
        {
            unsigned int connectedNodeNumber = connectedNode->getNumber();

            // Add all segments to the vector
            if (currentNodeNumber < connectedNodeNumber)
            {
                segments.emplace_back(std::min(currentNodeNumber, connectedNodeNumber), std::max(currentNodeNumber, connectedNodeNumber));
            }
        }
    }

    //Sort segments vector
    std::sort(segments.begin(), segments.end());

    //Remove duplicates with unique and erase, similar to file read in
    segments.erase(std::unique(segments.begin(), segments.end()), segments.end());

    //Output file dimensions as per PSLG file format: Number of segments, dimensions (2) and characteristics (1, colour)
    outputFile << segments.size() << " " << 2 << " " << 1 << std::endl;

    //Output unique segments with identifying number
    unsigned int seg_number = 1;
    Node* currentNode = &nodes[endNodeNumberIn - 1];

    //For each semgent
    for (auto& segment : segments)
    {
        bool isSegmentOnShortestPath = false;

        // Check if both nodes of the segment lie on the shortest path
        Node* node1 = &nodes[segment.first - 1];
        Node* node2 = &nodes[segment.second - 1];
        currentNode = &nodes[endNodeNumberIn - 1];
        while (currentNode != nullptr) {
            if ((currentNode == node1 && currentNode->getPrevNode() == node2) || (currentNode == node2 && currentNode->getPrevNode() == node1)) {
                isSegmentOnShortestPath = true;
                break;
            }
            //Else
            currentNode = currentNode->getPrevNode();
        }
        outputFile << seg_number << " " << segment.first << " " << segment.second << " " << isSegmentOnShortestPath << std::endl;
        seg_number++;
    }

    outputFile.close();
    std::cout << "PSLG segments have been written to output file" << std::endl;
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

    //Perform Dijkstra's shortest path algorithm
    std::clock_t start = std::clock();
    while(!unvisitedNodes.empty())
    {
        //Find unvisited node with shortest path length
        auto minNodeIt = std::min_element(unvisitedNodes.begin(), unvisitedNodes.end(), shortestPathComparator);
        currentNode = *minNodeIt;
        currentNode->setVisited();
        //Update the path lengths of the current node's connected nodes
        currentNode->updateConnectedNodeLengths();
        //Remove the current node from the unvisited nodes vector
        unvisitedNodes.erase(minNodeIt);
    }
    std::clock_t end = std::clock();

    // Find and print the path to end node. This debug code is included as it may be useful for future developers
    /*std::cout << "\nShortest path to Node "<< endNodeNumber <<":\n";
    currentNode = &nodes[endNodeNumber - 1];
    while (currentNode != nullptr) {
        currentNode->printNode();
        currentNode = currentNode->getPrevNode();
    }*/


    std::cout << "Shortest path length: " << nodes[endNodeNumber - 1].getShortestPathToNode() << std::endl;
    std::cout << "Time: " << (double)(((end - start) / (double)CLOCKS_PER_SEC)) << "s" << '\n';

    // Output PSLG to a file
    outputFile(nodes,endNodeNumber);

    return 0;
}
