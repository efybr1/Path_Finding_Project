//---------------------------------------------------------------------------------------//
//------------------------------------Node.h---------------------------------------------//
//---------------------------------------------------------------------------------------//
// Author: Ben Richards

//------------------------------------Includes-------------------------------------------//
#include "node.h"

//-------------------------------Standard functions--------------------------------------//

// Constructor definition, initialises all variables, and sets the shortestPathToNode to be
// infinity (which in C++ is simply the max a double can store)
Node::Node(unsigned int num, double x_coord, double y_coord)
    : number(num), x(x_coord), y(y_coord), shortestPathToNode(std::numeric_limits<double>::max()), visited(false), prevNode(nullptr) {}

// Copy constructor assignment operator overload, simply creates an exact copy
Node& Node::operator=(const Node& existingNode) {
    if (this != &existingNode) { // If trying to = itself, don't.
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

//----------------------------------Sets and Gets----------------------------------------//
// Sets
void Node::setConnection(Node* connectedNode) {
    connectedNodes.push_back(connectedNode);
}
void Node::setShortestPathToNode(double length) {
    shortestPathToNode = length;
}
void Node::setVisited() {
    visited = true;
}
void Node::setWrapper(Wrapper* wrapper) {
    wrapperPtr = wrapper;
}
void Node::setTemp(){
    temp = shortestPathToNode;
}
// Gets
unsigned int Node::getNumber() const {
    return number;
}
double Node::getShortestPathToNode() const {
    return shortestPathToNode;
}
Node* Node::getPrevNode() const {
    return prevNode;
}
Wrapper* Node::getWrapper() const {
    return wrapperPtr;
}
std::vector<Node*> Node::getConnectedNodes() const {
    return connectedNodes;
}

//--------------------------------Utility functions--------------------------------------//
// Function to calculate the distance between this Node and another Node.
double Node::pythagoras(const Node* connectedNode) const {
    return std::sqrt(std::pow(connectedNode->x - x, 2) + std::pow(connectedNode->y - y, 2));
}
// Calculates the distances between this Node and all it's connected Nodes.
void Node::calculateEdges() {
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        double distance = pythagoras(connectedNodes[i]);
        connectedNodeLengths.push_back(distance);
    }
}
// Used in Dijkstra's algorithm to update all the Nodes connected to this one. It will
// Calculate the shortest path through this Node for each of the connected Nodes, and if
// that calculated distance is smaller than a Node's current shortest path, it will be
// overwritten.
void Node::updateConnectedNodeLengths() {
    for (size_t i = 0; i < connectedNodes.size(); ++i) {
        Node* connectedNode = connectedNodes[i];
        double newLength = shortestPathToNode + connectedNodeLengths[i]; // Calculate length to connected node through this node.

        // Only update if the new length is smaller than the newly calculated value - to maintain the shortest path
        if (newLength < connectedNode->shortestPathToNode) {
            connectedNode->shortestPathToNode = newLength;
            connectedNode->prevNode = this;
        }
    }
}
void Node::setAllConnectedNodesToMin(){
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        connectedNodes[i]->setTemp();
        connectedNodes[i]->setShortestPathToNode(-1);
    }
}

// Prints
void Node::printNode() const {
    std::cout << "Node number: " << number << ", x: " << x << ", y: " << y << ", shortest path length: " << shortestPathToNode << ", visited: " << visited << std::endl;
}
void Node::printConnectedNodes() const {
    if (!connectedNodes.empty()) {
        for (size_t i = 0; i < connectedNodes.size(); i++) {
            std::cout << "Connected Node: " << connectedNodes[i]->number << std::endl;
        }
    }
}
void Node::printDistancesToConnections() const {
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        double distance = connectedNodeLengths[i];
        std::cout << "Distance from node " << number << " to node " << connectedNodes[i]->number << ": " << distance << std::endl;
    }
}
