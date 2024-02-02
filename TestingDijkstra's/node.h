#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>   // Include for std::sort, std::unique
#include <math.h>


//TO ADD//

//Private variables (issue with comparator function) along with sets and gets if needed
//Destructor


class Node {
public:
    unsigned int number; //Node index
    double x, y, shortestPathToNode; //x, y position, also shortest path distance
    bool visited; //bool to make sure node isn't visited twice, this should be taken out for optimization
    Node* prevNode; //Previous node variable to be able to trace the shortest path back for output.

    std::vector<Node*> connectedNodes; //Connected node pointer array
    std::vector<double> connectedNodeLengths; //Connected node lengths

    //Constructor
    Node(unsigned int num, double x_coord, double y_coord);

    //Sets/gets
    void setConnection(Node* connectedNode);
    void setShortestPathToNode(double length);
    void setVisited();

    unsigned int getNumber();
    double getShortestPathToNode();

    // Assignment operator overload
    Node& operator=(const Node& existingNode);

    //Utility Functions
    void calculateEdges();
    double pythagoras(const Node* connectedNode) const;
    void updateConnectedNodeLengths();

    //Prints
    void printNode();
    void printConnectedNodes();
    void printDistancesToConnections();
};

#endif // NODE_H_INCLUDED
