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
        //std::cout << "Changed wrapper" << number << "'s node (" << nodePtr->getNumber()<< ") to point to wrapper " << nodePtr->getWrapper()->getNumber() << std::endl;
    }
}

Wrapper::Wrapper(const Wrapper& other) : number(other.number), nodePtr(other.nodePtr) {
    if (nodePtr)
        nodePtr->setWrapper(this); // Update the associated Node's wrapper pointer to point to this Wrapper
    //std::cout << "Copy changed wrapper" << number << "'s node (" << nodePtr->getNumber()<< ") to point to wrapper " << nodePtr->getWrapper()->getNumber() << std::endl;
}

Wrapper& Wrapper::operator=(const Wrapper& other) {
    if (this != &other)
        {
        // Swap the node pointers
        nodePtr = other.nodePtr;
        // Update the associated Node's wrapper pointers
        if (nodePtr){
            nodePtr->setWrapper(this);//std::cout << "1. Assignment changed wrapper" << number << "'s node (" << nodePtr->getNumber()<< ") to point to wrapper " << nodePtr->getWrapper()->getNumber() << std::endl;
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

//Swap function which just swaps the node pointers of two Wrapper objects
void Wrapper::swap(Wrapper& other) {
    std::swap(nodePtr, other.nodePtr);
}

void Wrapper::getNodesConnectedNodes(){

    std::cout<< "Success!" << std::endl;
}
