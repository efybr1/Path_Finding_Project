//---------------------------------------------------------------------------------------//
//File Manipulation Program
//Ben Richards

//---------------------------------Description-------------------------------------------//
//Program to add node and ele file together for viewing in "veeMeshViewer" software.
//Requires change to first line of the node file.
//And also the auto-generated comment in the node file to be deleted for some reason.
//Node and ele file then written to a .tri file and this is the output of the program.

//----------------------------------Includes---------------------------------------------//
#include <iostream>
#include <fstream>
#include <string>

int main() {
//----------------------------------Open files-------------------------------------------//

    // File paths as variables as they are used in multiple locations in the program
    const char* nodeFilePath = "C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/SquareThreeHoles.1.node";
    const char* eleFilePath = "C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/SquareThreeHoles.1.ele";
    const char* triFilePath = "C:/Users/richa/OneDrive - The University of Nottingham/Documents/A_Year 4 EEC/A_Project/Meshes/SquareThreeHoles/SquareThreeHoles.1.tri";

    //The ele and node file to be added to the tri file
    std::ifstream nodeFile(nodeFilePath);
    std::ifstream eleFile(eleFilePath);
    std::ofstream triFile(triFilePath);

    //If any didn't open, stop the program
    if (!nodeFile.is_open())
    {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }
    if (!eleFile.is_open())
    {
        std::cerr << "Error opening eleFile file." << std::endl;
        return 1;
    }
    if (!triFile.is_open())
    {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

//------------------------------Manipulate node file-------------------------------------//
    // To hold the first line
    std::string firstLine;

    //Read first line and remove the last 4 characters from the first line.
    //This works as when file opens the file pointer defaults to the first line.
    std::getline(nodeFile, firstLine);
    firstLine.erase(firstLine.length() - 4);

    // Add the number 2 at the end of the (modified) first line
    firstLine += "2";

//--------------------------------Output node file---------------------------------------//

    //Count number of lines in node file as we need to remove last line by not writing it to tri file.
    std::string line;
    int numberOfLines = 0;
    while (std::getline(nodeFile, line)) //Loops through all lines
    {
        numberOfLines++;
    }
    //Close and reopen the node file because otherwise bad things happen - think because file pointer
    nodeFile.close();
    nodeFile.open(nodeFilePath);

    // Write the modified first line to the output tri file
    triFile << firstLine << std::endl;
    // Write the rest of the lines to the tri file, ignoring first as its already there
    // And ignoring last line as we want that removed because bad things happen with viewer software if left in.
    int linecount = 0;
    while (std::getline(nodeFile, line)) //Loops through all lines
    {
        if(linecount!=0 && linecount<numberOfLines)
        {
            triFile << line << std::endl;
        }
        linecount++;
    }

//---------------------------------Output ele file---------------------------------------//

    // Append the entire contents of the ele file to the tri file
    while (std::getline(eleFile, line)) //Loops through all lines
    {
        triFile << line << std::endl;
    }

//-----------------------------------Close files-----------------------------------------//
    std::cout << "Node and ele successfully combined." << std::endl;

    nodeFile.close();
    eleFile.close();
    triFile.close();

    return 0;
}

