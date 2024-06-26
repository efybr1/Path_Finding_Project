//---------------------------------------------------------------------------------------//
//----------------------------------Wrapper.cpp------------------------------------------//
//---------------------------------------------------------------------------------------//
// Author: Ben Richards

//------------------------------------Includes-------------------------------------------//
#include "wrapper.h"

// Implementations of the Wrapper class, see Wrapper.h for details

//-------------------------------Standard functions--------------------------------------//

// Constructor definition, sets the Wrapper pointer in the associated Node to point back to
// this object, while also initialising the Wrapper's unique identifier and the pointer to
// it's associated Node object.
Wrapper::Wrapper(int num, Node* node) : number(num), nodePtr(node) {
    if (nodePtr)
    {
        nodePtr->setWrapper(this);
    }
}

//Standard copy constructor to create a copy of the object, making sure to maintain the wrapper-node link
Wrapper::Wrapper(const Wrapper& other) : number(other.number), nodePtr(other.nodePtr) {
    if (nodePtr)
        nodePtr->setWrapper(this); // Update the associated Node's wrapper pointer to point to this Wrapper
}

//This allows for the ovveride of the std::swap function to only swap the node pointers
//When Node A is set = to node B, all that happens is the node pointers are swapped and their wrapper-node links updated
Wrapper& Wrapper::operator=(const Wrapper& other) {
    if (this != &other)
        {
        // Swap the node pointers
        nodePtr = other.nodePtr;
        // Update the associated Node's wrapper pointers
        if (nodePtr){
            nodePtr->setWrapper(this);
            }
        }
    return *this;
}

//----------------------------------Sets and Gets----------------------------------------//
//Gets
unsigned int Wrapper::getNumber() const {
    return number;
}
Node* Wrapper::getNodePtr() const {
    return nodePtr;
}
//--------------------------------Utility functions--------------------------------------//

void Wrapper::updateNodeWrapperPtr() {
    nodePtr->setWrapper(this); // Update the associated Node's wrapper pointer to point to this Wrapper
}

// Print wrapper information, and that of it's associated Node if it has one
void Wrapper::printWrapper() const {
    std::cout << "Wrapper number: " << number;
    if (nodePtr != nullptr)
        {
        std::cout << ", Associated Node number: " << nodePtr->getNumber();
        }
    else
        {
        std::cout << ", No associated Node";
        }
    std::cout << std::endl;
}
//Return connected nodes
std::vector<unsigned int> Wrapper::getNodesConnectedNodes() const {

    std::vector<unsigned int> connectedNodeNumbers;
    if (nodePtr) {
        std::vector<Node*> connectedNodes = nodePtr->getConnectedNodes();
        for (const auto& connectedNode : connectedNodes) {
            connectedNodeNumbers.push_back(connectedNode->getNumber());
        }
    }
    return connectedNodeNumbers;
}
