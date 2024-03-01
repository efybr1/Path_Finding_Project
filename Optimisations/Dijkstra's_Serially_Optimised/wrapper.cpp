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
        nodePtr->setWrapper(this);
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

