//---------------------------------------------------------------------------------------//
//-----------------------------------Wrapper.h-------------------------------------------//
//---------------------------------------------------------------------------------------//
// Author: Ben Richards

// Wrapper is a class to 'wrap' around a Node object. This is to add a layer of abstraction
// in between the Node objects and the vector that holds them. This then allows make_heap
// to simply move the node pointers around inside the Wrapper objects (1 per node), to
// allow for access into the heap once it is constructed

#ifndef WRAPPER_H_INCLUDED
#define WRAPPER_H_INCLUDED

//------------------------------------Includes-------------------------------------------//
#include "node.h"


class Wrapper {
private:
    unsigned int number; // Unique identifier
    Node* nodePtr; // Pointer to the Node the wrapper will contain

public:
    // Constructor
    Wrapper(int num, Node* node);

    // Gets
    unsigned int getNumber() const;
    Node* getNodePtr() const;

    // Utility functions
    void printWrapper() const;

    //Swapper function for wrappers. Called only from an specialisation of std::swap
    void swap(Wrapper& other);
};

#endif // WRAPPER_H_INCLUDED
