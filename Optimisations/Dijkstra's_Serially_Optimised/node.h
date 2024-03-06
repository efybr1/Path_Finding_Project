//---------------------------------------------------------------------------------------//
//------------------------------------Node.h---------------------------------------------//
//---------------------------------------------------------------------------------------//
// Author: Ben Richards

#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

//------------------------------------Includes-------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

// Class to encapsulate the concept of a Node for creating a weighted map for Dijkstra's to
// solve. The class contains the unique identifier, coordinates, the Node's current shortest
// path from the designated start node, the location of the Node before it in this shortest
// path, and all the associated standard and utility functions.
// Note, this has been modified for use in a specialised priority queue where the Node is
// wrapped in a  custom 'Wrapper' object, so it also contains a pointer to this object

class Wrapper; // Forward declaration of Wrapper class

class Node {
private:
    unsigned int number; // Node index
    Wrapper* wrapperPtr; // Pointer to the Wrapper that will hold the node
    double x, y, shortestPathToNode; // x, y position, also shortest path distance
    bool visited; // bool to make sure node isn't visited twice, this should be taken out for optimisation
    Node* prevNode; // Previous node variable to be able to trace the shortest path back for output

    double temp;

    std::vector<Node*> connectedNodes; // Connected node pointer array
    std::vector<double> connectedNodeLengths; // Connected node lengths

    double pythagoras(const Node* connectedNode) const;

public:
    // Constructor
    Node(unsigned int num, double x_coord, double y_coord);

    // Copy constructor assignment operator overload
    Node& operator=(const Node& existingNode);

    // Sets
    void setConnection(Node* connectedNode);
    void setShortestPathToNode(double length);
    void setVisited();
    void setWrapper(Wrapper* wrapper);
    void setTemp();

    // Gets
    unsigned int getNumber() const;
    double getShortestPathToNode() const;
    Node* getPrevNode() const;
    Wrapper* getWrapper() const;
    std::vector<Node*> getConnectedNodes() const;

    // Utility Functions
    void calculateEdges();
    void updateConnectedNodeLengths();
    void setAllConnectedNodesToMin();

    // Prints
    void printNode() const;
    void printConnectedNodes() const;
    void printDistancesToConnections() const;
};

#endif // NODE_H_INCLUDED
