#include "node.h"

//Constructor
Node::Node(unsigned int num, double x_coord, double y_coord)
    : number(num), x(x_coord), y(y_coord), shortestPathToNode(std::numeric_limits<double>::max()), visited(false), prevNode(nullptr) {}

// Assignment operator overload
Node& Node::operator=(const Node& existingNode) {
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

//Sets/gets
void Node::setConnection(Node* connectedNode) {
    connectedNodes.push_back(connectedNode);
}

void Node::setShortestPathToNode(double length) {
    shortestPathToNode = length;
}

void Node::setVisited() {
    visited = true;
}

unsigned int Node::getNumber() {
    return number;
}

//Utility Functions
void Node::calculateEdges() {
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        double distance = pythagoras(connectedNodes[i]);
        connectedNodeLengths.push_back(distance);
    }
}

double Node::pythagoras(const Node* connectedNode) const {
    return std::sqrt(std::pow(connectedNode->x - x, 2) + std::pow(connectedNode->y - y, 2));
}

void Node::updateConnectedNodeLengths() {
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
void Node::printNode() {
    std::cout << "Node number: " << number << ", x:  " << x << ", y: " << y << ", shortest path length: " << shortestPathToNode << ", visited:  " << visited << std::endl;
}

void Node::printConnectedNodes() {
    if (connectedNodes.size() > 0) {
        for (size_t i = 0; i < connectedNodes.size(); i++) {
            std::cout << "Connected Node: " << connectedNodes[i]->number << std::endl;
        }
    }
}

void Node::printDistancesToConnections() {
    for (size_t i = 0; i < connectedNodes.size(); i++) {
        double distance = connectedNodeLengths[i];
        std::cout << "Distance from node " << number << " to node " << connectedNodes[i]->number << ": " << distance << std::endl;
    }
}

