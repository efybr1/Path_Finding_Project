//Ben Richards

//---------------------------------Description-------------------------------------------//
//Current version of most recent project code


//----------------------------------Includes---------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <utility>
#include <algorithm>  // Include for std::sort, std::unique
#include <math.h>


//------------------------Classes------------------------------

//-------------------------------------------------------------
// Node: Node to encapsulate the data associated with the concept of a node.
//-------------------------------------------------------------
class Node {
public:
    unsigned int number; //Node index
    double x, y, shortestPathToNode; //x,y position, also shortest path distance
    bool visited; //bool to make sure node isn't visited twice, this should be taken out for optimisation
    Node* prevNode; //Previous node variable to be able to trace shortest path back for output.

    std::vector<Node*> connectedNodes; //Connected node pointer array
    std::vector<double> connectedNodeLengths;//Connected node lengths

    //Constructor
    Node(unsigned int num, double x_coord, double y_coord):number(num),x(x_coord),y(y_coord),shortestPathToNode(std::numeric_limits<double>::max()),visited(0),prevNode(nullptr){} //Full constructor to initialise index and coordinates

    //Sets/gets
    void setConnection(Node* connectedNode){connectedNodes.push_back(connectedNode);}
    void setShortestPathToNode(double length){shortestPathToNode=length;}
    void setVisited(){visited=1;}

    unsigned int getNumber(){return(number);}

    // Assignment operator overload
    Node& operator=(const Node& existingNode) {
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

    //Utility Functions
    void calculateEdges() {
        for (size_t i = 0; i < connectedNodes.size(); i++) {
            double distance = pythagoras(connectedNodes[i]);
            connectedNodeLengths.push_back(distance);
        }
    }
    double pythagoras(const Node* connectedNode) const {
        return std::sqrt(std::pow(connectedNode->x - x, 2) + std::pow(connectedNode->y - y, 2));
    }
    void updateConnectedNodeLengths() {
        for (size_t i = 0; i < connectedNodes.size(); ++i) {
            Node* connectedNode = connectedNodes[i];
            double newLength = shortestPathToNode + connectedNodeLengths[i]; //Calculate length to connected node through this node.

            // Only update if the new length is smaller than the newly calculated value - to maintain shortest path
            if (newLength < connectedNode->shortestPathToNode) {
                connectedNode->shortestPathToNode = newLength;
                connectedNode->prevNode = this;
            }
        }
    }

    //Prints
    void printNode(){
        std::cout << "Node number: " << number << ", x:  " << x << ", y: " << y << ", shortest path length: " << shortestPathToNode << ", visited:  "<< visited << std::endl;
    }
    void printConnectedNodes(){
        if (connectedNodes.size()>0)
        {
            for(size_t i=0;i<connectedNodes.size();i++)
            {
                std::cout << "Connected Node: " <<connectedNodes[i]->number << std::endl;
            }
        }
    }
    void printDistancesToConnections(){
    for (size_t i = 0; i < connectedNodes.size(); i++)
        {
        double distance = connectedNodeLengths[i];
        std::cout << "Distance from node " << number << " to node " << connectedNodes[i]->number << ": " << distance << std::endl;
        }
    }
};

    auto shortestPathComparator = [](const Node* a, const Node* b) -> bool {
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
        std::cout << "Node: " << i+1 << std::endl;
        nodes[i].printNode();
        nodes[i].printConnectedNodes();
        nodes[i].printDistancesToConnections();
        printf("\n");
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
