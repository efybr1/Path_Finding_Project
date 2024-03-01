#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

//Class to encapsulate the concept of a node
class Node {
public:
    unsigned int number; //Node index
    double x, y; //x,y position
    std::vector<Node*> connectedNodes; //Connected node pointer array
    Node(unsigned int num, double x_coord, double y_coord):number(num),x(x_coord),y(y_coord){} //Basic constructor to initialise index and coordinates
    void setConnection(Node* connectedNode)
    {
        connectedNodes.push_back(connectedNode);
    }
    void print()
    {
        std::cout << "Number: " << number << ", x:  " << x << ", y: " << y << std::endl;
    }
    void printConnectedNodes()
    {
        if (connectedNodes.size()>0)
        {
            for(int i=0;i<connectedNodes.size();i++)
            {
                std::cout << "Connected Node: " <<connectedNodes[i]->number << std::endl;
            }
        }
    }
};

//Node needs to have a vector of the edges


int main() {

    Node a = Node(1,1,1);
    Node b = Node(2,2,2);
    Node c = Node (3,3,3);

    a.print();
    a.setConnection(&b);
    a.setConnection(&c);
    a.printConnectedNodes();

    //The question is how to efficiently parse the ele and node files for the connections.
    //Triangles list these connections but repeated?
    return 0;
}
