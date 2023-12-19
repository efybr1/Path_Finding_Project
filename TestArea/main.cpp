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

    int startnode = 1;
    int endnode = 31;
    std::vector<Node *> visitedNodes;

    Node* endNode = &nodes[endnode-1];
    Node* currentNode = &nodes[startnode-1];

    currentNode->printNode();
    endNode->printNode();

    currentNode->setVisited();
    currentNode->setShortestPathToNode(0);
    currentNode->printNode();
    currentNode->printConnectedNodes();

    //
    nodes[1].printNode();
    currentNode->updateConnectedNodeLengths();
    nodes[1].printNode();

    std::cout<< "Dijkstra's: "<<std::endl;

    //Psuedocode

    //While != visited nodes aka still nodes unvisited
        //Find unvisited with shortestpathtonode at this point
        //Update neighbours
        //Iterate round

    //Test case should include example where the first example would have failed but true Dijkstra's wont



    //When you push on heap, it will work out where to put the element in its data structure such that the heap is maintained
    //"Is heap" should test if it is good
    //In an ideal world we would just continually pop to get the element
    //However we not only pop the element, we need to update nodes which are still in the heap.
    //Calling makeheap again is appaling
        //Right way is to remove that element from the heap and re-insert it - but no command to remove from the middle of the heap, and we dont know
        //Where it is. Each node should contain a pointer into the heap to know where it is in the heap.
        //A node should know how to remove itself from the heap
        //Now need to go find some computer science slides to explain how pop and push heap work. Going to have to interact/modify them
        //Go into the standard STL header codes

    //Runtime is primarily on number of nodes.
    //So take 1 example and subdivide triangles over and over to prove this
    //Then take very different example and do the same, hopefully see the same scaling with the nodes

    //Take measurements for n^2 version to have a baseline.
    //
    //Profiler GPROF - on third year moodle page - profiler which will tell you how many ms used in each line of code




    // Update until at end node
    while (currentNode != endNode) {
        currentNode->printNode();
        currentNode->updateConnectedNodeLengths();

        // Mark the current node as visited
        currentNode->setVisited();
        visitedNodes.push_back(currentNode);

        // Find the next unvisited node with the shortest path - this section is what likely needs changed to shortest path from the big list of unvisited
        double minPath = std::numeric_limits<double>::max();
        Node* nextNode = nullptr;

        //Update all the connected nodes shortest paths and select next as shortest path
        for (Node* connectedNode : currentNode->connectedNodes)
            {
            if (!connectedNode->visited && connectedNode->shortestPathToNode < minPath)
                {
                minPath = connectedNode->shortestPathToNode;
                nextNode = connectedNode;
                }
            }

        if (nextNode == nullptr) {
            break;
        }
        currentNode = nextNode;
    }

    std::cout<< "Nodes: "<<std::endl;
    for(size_t i=0;i<nodes.size();i++)
    {
        nodes[i].printNode();
    }

    // Print the shortest path
    std::cout << "Shortest Path:" << std::endl;
    const Node* currentPathNode = endNode;
    while (currentPathNode != nullptr) {
        std::cout << "Node: " << currentPathNode->number << ", Shortest Path Length: " << currentPathNode->shortestPathToNode << std::endl;
        currentPathNode = currentPathNode->prevNode;
    }


    //Encapsulate get next node - initially naive way. But write the code around next node function so that we only
    //need to change the function.
    //If we want a rigerous search. We have a container of nodes. Our task is to access the one with the smallest value repeatedly.
    //We have n mandetory. Then linear search is n. If we store sorted - it would be nlogn. However we dont need to entirely sort.
    //Not doing this one off - doing it in an iterative sense - we want to know the smallest and then pop. And then we come back for the next
    //Smallest, remove and so on...

    //Apparently a standard algorithm container from this

    //We have values in the container we must find, however some values change each loop - we will need to modify the container to be
    //able to remove these values and re-add them in a way so we can always find the smallest.

    //Even on input routines - get into habit of getting order n analysis of them - showing consistent set of optimisation strats im applying
    //to my code.
    return 0;
}
