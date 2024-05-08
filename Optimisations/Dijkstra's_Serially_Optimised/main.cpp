// Ben Richards

//---------------------------------Description-------------------------------------------//
// Testing Dijkstra's naive implementation of Dijkstra's

//----------------------------------Includes---------------------------------------------//

#include <fstream> //Used for file inputs and outputs
#include <sstream>
#include <math.h> //Used for calculations of path lengths
#include <iostream> //Used from output to command window
#include <vector> //Used for data structures
#include <utility> //Utility functions
#include <algorithm>  // Include for std::sort, std::unique
#include <ctime> // Include for std::sort, std::unique
#include <queue> //Included for queue class

#include "wrapper.h" //Wrapper class

//----------------------------------Functions--------------------------------------------//

//---------------------------------------------------------------------------------------//
//----------------------------shortestPathComparators-------------------------------------//
//---------------------------------------------------------------------------------------//

// Function used to compare two nodes. Returns boolean of whether a is < b
auto shortestPathComparator(Node* a, Node* b) {
    return a->getShortestPathToNode() > b->getShortestPathToNode();
};

// Lambda function for comparing wrappers -  just use the node one
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
    std::ifstream NodeInputFile("C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/MainProgramTestMeshes/04.node");
    if (!NodeInputFile.is_open())
    {
        std::cout<< "Error opening node file!" << std::endl;
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
    // Open output file
    std::ofstream outputFile("FinalTest.plsg");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Output nodes
    // First line
    outputFile << nodes.size() << " " << 2 << " " << 0 << std::endl;
    // All nodes
    for (const auto& node : nodes)
    {
        outputFile << node.getNumber() << " " << node.getX() << " " << node.getY() << std::endl;
    }

    // Output edges/segments
    // First line
    std::vector<std::pair<unsigned int, unsigned int>> segments;

    // Find all segments
    for (Node node : nodes)
    {
        unsigned int currentNodeNumber = node.getNumber();
        std::vector<Node*> connectedNodes = node.getConnectedNodes();

        // Iterate through connected nodes
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

    // Sort segments vector
    std::sort(segments.begin(), segments.end());

    // Remove duplicates with unique and erase, similar to file read in
    segments.erase(std::unique(segments.begin(), segments.end()), segments.end());

    // Output file dimensions as per PSLG file format: Number of segments, dimensions (2) and characteristics (1, colour)
    outputFile << segments.size() << " " << 2 << " " << 1 << std::endl;

    // Output unique segments with identifying number
    unsigned int seg_number = 1;

    // Get the end node
    Node* endNode = &nodes[endNodeNumberIn - 1];

    // For each segment
    for (auto& segment : segments)
    {
        bool isSegmentOnShortestPath = false;

        // Iterate over each node on the shortest path
        Node* currentNode = endNode;
        while (currentNode->getPrevNode() != nullptr)
        {
            Node* nextNode = currentNode->getPrevNode();

            // Check if the segment formed by currentNode and nextNode matches the current segment being evaluated
            if ((currentNode->getNumber() == segment.first && nextNode->getNumber() == segment.second) ||
                (currentNode->getNumber() == segment.second && nextNode->getNumber() == segment.first))
            {
                isSegmentOnShortestPath = true;
                break;
            }

            currentNode = nextNode;
        }

        // Output segment with attribute
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
    nodes[startNodeNumber-1].setShortestPathToNode(0);

    //std::cout << "\n--------------------------------------Input nodes--------------------------------------" << std::endl;

    // Create a vector of Wrapper objects, one for each node
    std::vector<Wrapper> pq;
    for (unsigned int i=0; i<nodes.size();i++)
    {
        pq.push_back(Wrapper(i+1, &nodes[i]));
    }

    // Convert vector into a priority queue using make_heap
    std::make_heap(pq.begin(), pq.end(), compareByShortestPath);

    //Start timing
    std::clock_t start = std::clock();
    for(unsigned int i = 0; i<nodes.size(); i++)
    {
        std::pop_heap(pq.begin(), pq.end(), compareByShortestPath); //Pop current shortest node
        Wrapper topElement = pq.back(); //Get the wrapper
        pq.pop_back(); //Remove the last wrapper - and the node which we won't re-add as it is visited
        topElement.getNodePtr()->setVisited(); //Set node to visited
        std::vector<unsigned int> connectedNodes = topElement.getNodesConnectedNodes(); //Get Node's connected nodes

        for (unsigned int i = 0; i < connectedNodes.size(); ++i)
        {
            if(nodes[connectedNodes[i]-1].getVisited()==0) //If the connected node has not been visited, it will still be on the stack, so must be removed and re-added
            {
                nodes[connectedNodes[i]-1].setShortestPathToMin(); // Need to use this so the value is stored in temp
                std::push_heap(pq.begin(),pq.begin()+nodes[connectedNodes[i]-1].getWrapper()->getNumber(),compareByShortestPath); //Moves the node to the top
                std::pop_heap(pq.begin(), pq.end(), compareByShortestPath); //Moves the node to the last wrapper (back of the vector)
                Wrapper removed = pq.back();
                pq.pop_back(); //Removes the last wrapper which now contains the node

                topElement.getNodePtr()->updateConnectedNodeLength(i); //Update this connected node's length
                pq.push_back(removed); //Re add to heap
                std::push_heap(pq.begin(),pq.end(),compareByShortestPath); //Push it into the correct place

            }
            else //It has been visited, so simply update value if smaller but do not re-add
                {
                    nodes[connectedNodes[i]-1].updateConnectedNodeLengths();
                }
        }
    }
    std::clock_t end = std::clock();

    std::cout << "\n--------------------------------------Dijkstra's using priority queue--------------------------------------" << std::endl;

    std::cout << "Shortest path length: " << nodes[endNodeNumber - 1].getShortestPathToNode() << std::endl;
    std::cout << "Time: " << (double)(((end - start) / (double)CLOCKS_PER_SEC)) << "s" << '\n';

    // Find and print the path to end node. This debug code is included as it may be useful for future developers
    /*std::cout << "\nShortest path to Node "<< endNodeNumber <<":\n";
    Node* currentNode = &nodes[endNodeNumber - 1];
    while (currentNode != nullptr) {
    currentNode->printNode();
        currentNode = currentNode->getPrevNode();
    }*/

    outputFile(nodes,endNodeNumber);

    return 0;
}
