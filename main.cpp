/**
 * CLion 2021.2.3<br>
 * Licensed to Gilberto Rouxinol<br>
 * For educational use only.<br>
 * <p>
 * Polytechnic Institute of Viseu<br>
 * School of Technology and Management of Viseu<br>
 * <p>
 * C++ File created by Gilberto Rouxinol on 2022<br>
 * Copyright © 2022 Gilberto Rouxinol<br>
 * All rights reserved<br><br>
 * <p>
 *
 * @author Gilberto Rouxinol
 * @version 05.05.2022.2
 */

/**
 * Website with all numerical libraries in use here:
 * https://people.math.sc.edu/Burkardt/cpp_src/cpp_src.html
 * r8lib.cpp and r8lib.hpp
 */

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sstream>


#include "r8lib.hpp"

#if defined(_WIN32) || defined(__CYGWIN__) // Windows (x86 or x64)
#define ENVIRONMENT "WINDOWS"
#elif defined(__linux__) // Linux
#define ENVIRONMENT "LINUX"
#elif defined(__APPLE__) && defined(__MACH__) // Mac OS
#define ENVIRONMENT "MACOS"
#elif defined(unix) || defined(__unix__) || defined(__unix) // Unix like OS
#define ENVIRONMENT "UNIX"
#else
#error Unknown environment!
#define ENVIRONMENT ERROR
#endif

#define SMALL_VALUE 0.0001

using namespace::std;




string *splitString(string strIn, char delimiter){
    // Figure out the dimension of the strings array spl[ ] by counting the tokens in the string strIn
    int numberTokens = 0;
    string strOut;
    string aux = strIn;
    stringstream ssAux(aux); // Breaking input into token
    while(getline(ssAux, strOut, delimiter)) {
        ++numberTokens;
    }

    // Building the strings array spl[ ]
    string *spl;
    spl = new string[numberTokens];
    stringstream ss(strIn); // Breaking input into token
    int i = 0;
    while(getline(ss, strOut, delimiter)){
        spl[i] = strOut;
        ++i;
    }

    return spl;
}




double equationMEd(double a, double b, double c, double x){ // a = distLoad[i] | b = vectEd[2] | c = vectEd[3]
    return -0.5*a*pow(x,2) + b*x + c;
}




double equationVEd(double a, double b, double x){ // a = distLoad[i] | b = vectEd[2]
    return -a*x + b;
}




double equationNEd(double a, double x){ // a = vectEd[1]
    return a;
}




double zeroInVEd(double a, double b, double lh){ // a = distLoad[i] | b = vectEd[2] | lh = length | x must be: 0 < x < length
    if( (abs(a) <= SMALL_VALUE) || (b/a < 0) || (b/a > lh))
        return -1; //not possible
    else
        return b/a;
}




double *zeroInMEd(double a, double b, double c, double lh){ // a = distLoad[i] | b = vectEd[2] | c = vectEd[3]
    double *sol; sol = new double[2]; sol[0] = 0; sol[1] = 0;
    if (abs(a) < SMALL_VALUE && abs(b) > SMALL_VALUE){
        sol[0] = -c/b;
        sol[1] = -1;
        return sol;
    }else if (abs(a) < SMALL_VALUE && abs(b) < SMALL_VALUE){
        sol[0] = -1;
        sol[1] = -1;
        return sol;
    }else{
        sol[0] = ( -b + sqrt(pow(b,2)-4*(-a/2)*c) )/( 2*(-a/2) );
        sol[1] = ( -b - sqrt(pow(b,2)-4*(-a/2)*c) )/( 2*(-a/2) );
        if(sol[0] < 0){sol[0] = -1;}
        if(sol[1] < 0){sol[1] = -1;}
        if(sol[0] > lh){sol[0] = -1;}
        if(sol[1] > lh){sol[1] = -1;}
        return sol;
    }
}




double deltaDeforma(double lgth, double modE, double inertia, double m1, double m2, double x){
    return pow(lgth,2)*(m1*(2*x-3*pow(x,2)+pow(x,3)) + m2*(x-pow(x,3)))/(6*modE*inertia);
}




int main(int argc, const char * argv[]) {




//**********************************************************************************************************************
//      Read Input file save in the InputStructures folder
//**********************************************************************************************************************
    string nameIn;
    string pathFileIn = "../InputStructures/";

    cout << "Enter the name of the Input file saved in" << endl;
    cout << pathFileIn << endl;
    cout << "folder: " << endl;
    getline(cin, nameIn);
    pathFileIn = pathFileIn + nameIn;

    ifstream myReadFile(pathFileIn, ios::in);
    if(!myReadFile){
        cout << "File not found" << endl; exit(0);
    }




//**********************************************************************************************************************
//      Creating output file and saving it in the OutputStructures folder
//**********************************************************************************************************************
    string nameOut;
    string pathFileOutFolder = "../OutputStructures/";
    string pathFileOutGeneral;
    string pathFileOutDeformabilityPython;
    string pathFileOutDeformabilityLatexMP;
    string pathFileOutMomentLatexMP;
    string pathFileOutShearLatexMP;
    string pathFileOutNormalLatexMP;

    cout << "The Output files will be saved in the folder:" << endl;
    cout << pathFileOutFolder << endl;

    nameIn[0] = toupper(nameIn[0]);
    nameOut = "out" + nameIn;
    pathFileOutGeneral = pathFileOutFolder + nameOut;
    pathFileOutDeformabilityPython = pathFileOutFolder + "pythonOutDeformation" + nameIn;
    pathFileOutDeformabilityLatexMP = pathFileOutFolder + "latexMPOutDeformation" + nameIn;
    pathFileOutMomentLatexMP = pathFileOutFolder + "latexMPOutMoment" + nameIn;
    pathFileOutShearLatexMP = pathFileOutFolder + "latexMPOutShear" + nameIn;
    pathFileOutNormalLatexMP = pathFileOutFolder + "latexMPOutNormal" + nameIn;


    ofstream svCout(pathFileOutGeneral, ios::out);
    ofstream svPythonDefCout(pathFileOutDeformabilityPython, ios::out);
    ofstream svLatexDefCout(pathFileOutDeformabilityLatexMP, ios::out);
    ofstream svLatexMomCout(pathFileOutMomentLatexMP, ios::out);
    ofstream svLatexSheCout(pathFileOutShearLatexMP, ios::out);
    ofstream svLatexNorCout(pathFileOutNormalLatexMP, ios::out);




//**********************************************************************************************************************
//      Geometric characteristics, material properties, structural support and loads
//**********************************************************************************************************************
    string lineRead;
    getline (myReadFile, lineRead);
    getline (myReadFile, lineRead);
    int readMode = stoi(lineRead);
    svCout << "Mode of the data input:  1 automatic processing |  2 pre-processing" << endl;
    svCout << readMode << endl;




    //  Reading strutures
    int numberNodes;
    double **nodes;

    int numberElements;
    int **elements;

    double modE;
    double modG;

    int numberMaterialsAreas;
    double *materialsAreasValues;
    double *materialsArea;

    int numberMaterialsInertias;
    double *materialsInertiasValues;
    double *materialsInertia;

    int numberFixities;
    int **fixities;
    int *nodesFix;

    int numberDistLoad;
    double *distLoad;

    int numberNodalLoad;
    double **nodalLoad;




    //  Calculating strutures
    int **degreesOfFreedom;

    double **stiffnessMatrix;
    double *f0LoadVector;
    double *fLoadVector;

    int **incidence;

    double **localStiffnessMatrix;
    double **globalStiffnessMatrix;

    double *localF0LoadVector;
    double *globalF0LoadVector;

    double **transformation;
    double **transformationTranspose;


    double zoomFactorD;
    double zoomFactorM;
    double zoomFactorV;
    double zoomFactorN;


    //******************************************************************************************************************
    // 1 - Automatic Processing ****************************************************************************************
    //******************************************************************************************************************
    if ( readMode == 1 ){
        //**************************************************************************************************************
        //      SI units system and structure schema
        //**************************************************************************************************************
        //  m    - meter
        //  kPa  - kiloPascal
        //  kN   - kiloNewton
        //
        //  Matrix and vector numbering starts with 0 not with 1, e.g., node 1 is 0, node 2 is 1, ..., node 16 is 15
        //
        //
        //                   Numbering the nodes                                           Numbering the elements
        //
        //           13----------14----------15----------16                      -----19----------20-----------21-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                      15          16          17          18
        //           |           |           |           |                       |           |           |           |
        //           9----------10----------11----------12                       -----12-----------13----------14-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       8           9          10          11
        //           |           |           |           |                       |           |           |           |
        //           5-----------6-----------7-----------8                       ------5-----------6------------7-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       1           2           3           4
        //           |           |           |           |                       |           |           |           |
        //         --1--       --2--       --3--       --4--                   -----       -----       -----       -----
        //
        //
        //                   Load the nodes - F                                          Load the elements - F0
        //
        //                                                                          50 kN/m      50 kN/m    50 kN/m
        //   60 kN                                                                vvvvvvvvvvv vvvvvvvvvvv vvvvvvvvvvv
        //   ---->   -------------------------------------                       -------------------------------------
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       |  50 kN/m  |  50 kN/m  |  50 kN/m  |
        //   40 kN   |           |           |           |                       |vvvvvvvvvvv|vvvvvvvvvvv|vvvvvvvvvvv|
        //   ---->   -------------------------------------                       -------------------------------------
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       |  50 kN/m  |  50 kN/m  |  50 kN/m  |
        //   20 kN   |           |           |           |                       |vvvvvvvvvvv|vvvvvvvvvvv|vvvvvvvvvvv|
        //   ---->   -------------------------------------                       -------------------------------------
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       |           |           |           |
        //         -----       -----       -----       -----                   -----       -----       -----       -----
        //
        //
        //
        //
        //                      Numbering the DOF
        //
        //       25 26 27 ----28 29 30 ---- 31 32 33----34 35 36                 -----19----------20-----------21-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                      15          16          17          18
        //           |           |           |           |                       |           |           |           |
        //       13 14 15-----16 17 18 ----19 20 21----22 23 24                  -----12-----------13----------14-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       8           9          10          11
        //           |           |           |           |                       |           |           |           |
        //        1 2 3--------4 5 6-------7 8 9------10 11 12                   ------5-----------6------------7-----
        //           |           |           |           |                       |           |           |           |
        //           |           |           |           |                       1           2           3           4
        //           |           |           |           |                       |           |           |           |
        //         -----       -----       -----       -----                   -----       -----       -----       -----
        //


        //**************************************************************************************************************
        //      Geometry
        //**************************************************************************************************************
        int nSpan;
        int nStory;
        double heightStory;
        double widthFrame;

        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number of span
        getline (myReadFile, lineRead);
        nSpan = stoi(lineRead);
        svCout << nSpan << endl;

        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number of floors
        getline (myReadFile, lineRead);
        nStory = stoi(lineRead);
        svCout << nStory << endl;


        getline (myReadFile, lineRead);
        svCout << lineRead << lineRead << endl; // Height floor [m]
        getline (myReadFile, lineRead);
        heightStory = stod(lineRead);
        svCout << heightStory << endl;

        getline (myReadFile, lineRead);
        svCout << lineRead << lineRead << endl; // Width frame [m]
        getline (myReadFile, lineRead);
        widthFrame = stod(lineRead);
        svCout << widthFrame << endl;




        //**************************************************************************************************************
        //      Nodes coordinates
        //**************************************************************************************************************
        numberNodes = (nStory+1)*(nSpan+1);
        nodes = new double *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            nodes[i] = new double[2];
        }
        int k = 0;
        for (int i = 0; i < nStory+1; ++i) {
            for (int j = 0; j < nSpan+1; ++j) {
                nodes[k][0] = widthFrame*j;
                nodes[k][1] = heightStory*i;
                ++k;
            }
        }

        svCout << "Nodes Coordinates" << endl;
        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1 << setw(12) << nodes[i][0] << setw(12) << nodes[i][1] << endl;
        }




        //**************************************************************************************************************
        //      Elements nodes
        //**************************************************************************************************************
        numberElements = (nSpan+1)*(nStory) + (nSpan)*(nStory);
        elements = new int *[numberElements];
        for (int i = 0; i < numberElements; i++ ){
            elements[i] = new int[2];
        }
        int node = 1;
        int elem = 0;
        for (int j = 0; j < nStory; ++j) {
            for (int i = node; i < node + nSpan + 1; ++i) {
                elements[elem][0] = i;
                elements[elem][1] = elements[elem][0] + nSpan + 1;
                ++elem;
            }
            node = node + nSpan + 1;
            for (int i = node; i < node + nSpan; ++i) {
                elements[elem][0] = i;
                elements[elem][1] = elements[elem][0] + 1;
                ++elem;
            }
        }

        svCout << "Elements nodes" << endl;
        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << elements[i][0] << setw(12) << elements[i][1] << endl;
        }


        //**************************************************************************************************************
        //      Constants
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Elasticity modulus [kPa]
        getline (myReadFile, lineRead);
        modE = stod(lineRead);
        svCout << modE << endl;

        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Shear modulus [kPa]
        getline (myReadFile, lineRead);
        modG = stod(lineRead);
        svCout << modG << endl;




        //**************************************************************************************************************
        //      Elements Areas
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number areas
        getline (myReadFile, lineRead);
        numberMaterialsAreas = stoi(lineRead);
        svCout << numberMaterialsAreas << endl;

        materialsAreasValues = new double [numberMaterialsAreas];
        for (int i = 0; i < numberMaterialsAreas; i++ ){
            materialsAreasValues[i] = 0.0;
        }
        for (int i = 0; i < numberMaterialsAreas; ++i) {
            getline (myReadFile, lineRead);
            svCout << lineRead << endl; // Area i [m2]
            getline (myReadFile, lineRead);
            materialsAreasValues[i] = stod(lineRead);
            svCout << materialsAreasValues[i] << endl;
        }

        materialsArea = new double [numberElements];
        for (int i = 0; i < numberElements; i++ ){
            materialsArea[i] = 0.0;
        }
        int column = 0;
        int beam = nSpan + 1;
        for (int j = 0; j < nStory; ++j) {
            for (int i = column; i < beam; ++i) {
                materialsArea[i] = materialsAreasValues[0];
            }
            column = beam + nSpan;
            for (int i = beam; i < beam+nSpan; ++i) {
                materialsArea[i] = materialsAreasValues[1];
            }
            beam += 2*nSpan+1;
        }

        svCout << "Elements area" << endl;
        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << materialsArea[i] << endl;
        }




        //**************************************************************************************************************
        //      Elements Inertia
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number inertias
        getline (myReadFile, lineRead);
        numberMaterialsInertias = stoi(lineRead);
        svCout << numberMaterialsInertias << endl;

        materialsInertiasValues = new double [numberMaterialsInertias];
        for (int i = 0; i < numberMaterialsInertias; i++ ){
            materialsInertiasValues[i] = 0.0;
        }
        for (int i = 0; i < numberMaterialsInertias; ++i) {
            getline (myReadFile, lineRead);
            svCout <<  lineRead << endl; // Inertia i [m4]
            getline (myReadFile, lineRead);
            materialsInertiasValues[i] = stod(lineRead);
            svCout << materialsInertiasValues[i] << endl;
        }

        materialsInertia = new double [numberElements];
        for (int i = 0; i < numberElements; i++ ){
            materialsInertia[i] = 0.0;
        }
        column = 0;
        beam = nSpan + 1;
        for (int j = 0; j < nStory; ++j) {
            for (int i = column; i < beam; ++i) {
                materialsInertia[i] = materialsInertiasValues[0];
            }
            column = beam + nSpan;
            for (int i = beam; i < beam + nSpan; ++i) {
                materialsInertia[i] = materialsInertiasValues[1];
            }
            beam += 2 * nSpan + 1;
        }

        svCout << "Elements inertia" << endl;
        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << materialsInertia[i] << endl;
        }




        //**************************************************************************************************************
        //      Structural support or fixities
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number fixities
        getline (myReadFile, lineRead);
        numberFixities = stoi(lineRead);
        svCout << numberFixities << endl;

        fixities = new int *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            fixities[i] = new int[3];
        }
        for (int i = 0; i < numberNodes; i++ ){
            fixities[i][0] = 0;
            fixities[i][1] = 0;
            fixities[i][2] = 0;
        }

        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Node(s) fixitie(s)

        nodesFix = new int [numberFixities];
        for (int i = 0; i < numberFixities; i++ ){
            getline (myReadFile, lineRead);
            nodesFix[i] = stoi(lineRead);
        }

        for (int i = 0; i < numberFixities; ++i) {
            svCout << setw(5) << i+1 << setw(12) << nodesFix[i] << endl;
        }

        for (int i = 0; i < numberFixities; ++i) {
            fixities[i][0] = 1; fixities[i][1] = 1; fixities[i][2] = 1;
        }

        svCout << "Fixitie(s)" << endl;
        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1 << setw(12) << fixities[i][0] << setw(12) << fixities[i][1] <<setw(12) << fixities[i][2] <<endl;
        }




        //**************************************************************************************************************
        //      Distributed loads - To build the vector F0
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Number elements with distributed loads
        getline (myReadFile, lineRead);
        numberDistLoad = stoi(lineRead);
        svCout << numberDistLoad << endl;

        distLoad = new double [numberElements];
        for (int i = 0; i < numberElements; i++ ){
            distLoad[i] = 0.0;
        }

        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Element and load [kN/m]

        string *arrayStringsDL;
        for (int i = 0; i < numberDistLoad; ++i) {
            getline (myReadFile, lineRead);
            arrayStringsDL = splitString(lineRead, ' ');
            distLoad[stoi(arrayStringsDL[0])-1] = stod(arrayStringsDL[1]);
        }

        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << distLoad[i] <<endl;
        }




        //**************************************************************************************************************
        //      Nodal loads - To build the vector F
        //**************************************************************************************************************
        getline (myReadFile, lineRead);
        svCout <<lineRead << endl; // Number nodes with nodal loads
        getline (myReadFile, lineRead);
        numberNodalLoad = stoi(lineRead);
        svCout << numberNodalLoad << endl;

        nodalLoad = new double *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            nodalLoad[i] = new double[3];
        }
        for (int i = 0; i < numberNodes; i++ ){
            nodalLoad[i][0] = 0;
            nodalLoad[i][1] = 0;
            nodalLoad[i][2] = 0;
        }


        getline (myReadFile, lineRead);
        svCout << lineRead << endl; // Node and load [kN]
        string *arrayStringsNL;
        for (int i = 0; i < numberNodalLoad; ++i) {
            getline (myReadFile, lineRead);
            arrayStringsNL = splitString(lineRead, ' ');
            nodalLoad[stoi(arrayStringsNL[0])-1][0] = stod(arrayStringsNL[1]);
            nodalLoad[stoi(arrayStringsNL[0])-1][1] = stod(arrayStringsNL[2]);
            nodalLoad[stoi(arrayStringsNL[0])-1][2] = stod(arrayStringsNL[3]);
        }

        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 3; ++j) {
                svCout << setw(12) << nodalLoad[i][j];
            }
            svCout << endl;
        }

        getline (myReadFile, lineRead);// Zoom Factor Deformability
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorD = stod(lineRead);
        svCout << zoomFactorD << endl;

        getline (myReadFile, lineRead);// Zoom Factor Normal Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorN = stod(lineRead);
        svCout << zoomFactorN << endl;

        getline (myReadFile, lineRead);// Zoom Factor Transverse Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorV = stod(lineRead);
        svCout << zoomFactorV << endl;

        getline (myReadFile, lineRead);// Zoom Factor Bending Moment Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorM = stod(lineRead);
        svCout << zoomFactorM << endl;

        svCout << "Reading completed successfully" << endl;
    }




    //******************************************************************************************************************
    // 2 - Pre-processing **********************************************************************************************
    //******************************************************************************************************************
    if ( readMode == 2 ){
        //**************************************************************************************************************
        //      SI units system and structure schema
        //**************************************************************************************************************
        //  m    - metro
        //  kPa  - kiloPascal
        //  kN   - kiloNewton
        //

        //**************************************************************************************************************
        //      Geometry - Nodes coordinates
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Number of nodes
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        numberNodes = stoi(lineRead);
        svCout << numberNodes << endl;

        nodes = new double *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            nodes[i] = new double[2];
        }

        getline (myReadFile, lineRead); // List of nodes
        svCout << lineRead << endl;
        string *strNodes;
        strNodes = new string[3];
        for (int i = 0; i < numberNodes; ++i) {
            getline (myReadFile, lineRead);
            strNodes = splitString(lineRead, ' ');
            nodes[stoi(strNodes[0])-1][0] = stod(strNodes[1]);
            nodes[stoi(strNodes[0])-1][1] = stod(strNodes[2]);
        }

        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1 << setw(12) << nodes[i][0] << setw(12) << nodes[i][1] << endl;
        }




        //**************************************************************************************************************
        //      Geometry - Elements nodes
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Number of Elements
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        numberElements = stoi(lineRead);
        svCout << numberElements << endl;

        elements = new int *[numberElements];
        for (int i = 0; i < numberElements; i++ ){
            elements[i] = new int[2];
        }

        getline (myReadFile, lineRead); // List of elements
        svCout << lineRead << endl;
        string *strElem;
        strElem = new string[3];
        for (int j = 0; j < numberElements; ++j) {
            getline (myReadFile, lineRead);
            strElem = splitString(lineRead, ' ');
            elements[stoi(strElem[0])-1][0] = stoi(strElem[1]);
            elements[stoi(strElem[0])-1][1] = stoi(strElem[2]);
        }

        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << elements[i][0] << setw(12) << elements[i][1] << endl;
        }




        //**************************************************************************************************************
        //      Constants
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Elasticity modulus [kPa]
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        modE = stod(lineRead);
        svCout << modE << endl;

        getline (myReadFile, lineRead); // Shear modulus [kPa]
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        modG = stod(lineRead);
        svCout << modG << endl;




        //**************************************************************************************************************
        //      Elements Areas
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Elements areas [m2]
        svCout << lineRead << endl;

        string *strAre;
        strAre = new string[2];

        materialsArea = new double [numberElements];
        for (int j = 0; j < numberElements; ++j) {
            getline (myReadFile, lineRead);
            strAre = splitString(lineRead, ' ');
            materialsArea[stoi(strAre[0])-1] = stod(strAre[1]);
        }

        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << materialsArea[i] << endl;
        }




        //**************************************************************************************************************
        //      Elements Inertia
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Elements inertia [m4]
        svCout << lineRead << endl;

        string *strIne;
        strIne = new string[2];

        materialsInertia = new double [numberElements];
        for (int j = 0; j < numberElements; ++j) {
            getline (myReadFile, lineRead);
            strIne = splitString(lineRead, ' ');
            materialsInertia[stoi(strIne[0])-1] = stod(strIne[1]);
        }

        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << materialsInertia[i] << endl;
        }




        //**************************************************************************************************************
        //      Structural support or fixities
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Fixitie(s)
        svCout << lineRead << endl;

        string *strFix;
        strFix = new string[4];

        fixities = new int *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            fixities[i] = new int[3];
        }

        for (int j = 0; j < numberNodes; ++j) {
            getline (myReadFile, lineRead);
            strFix = splitString(lineRead, ' ');
            fixities[stoi(strFix[0])-1][0] = stoi(strFix[1]);
            fixities[stoi(strFix[0])-1][1] = stoi(strFix[2]);
            fixities[stoi(strFix[0])-1][2] = stoi(strFix[3]);
        }

        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1 << setw(12) << fixities[i][0] << setw(12) << fixities[i][1] <<setw(12) << fixities[i][2] <<endl;
        }



        //**************************************************************************************************************
        //      Distributed loads - To build the vector F0
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Distributed loads [kN/m]
        svCout << lineRead << endl;

        string *strDistLoad;
        strDistLoad = new string[2];

        distLoad = new double [numberElements];
        for (int j = 0; j < numberElements; ++j) {
            getline (myReadFile, lineRead);
            strDistLoad = splitString(lineRead, ' ');
            distLoad[stoi(strDistLoad[0])-1] = stod(strDistLoad[1]);
        }

        for (int i = 0; i < numberElements; ++i) {
            svCout << setw(5) << i+1 << setw(12) << distLoad[i] <<endl;
        }




        //**************************************************************************************************************
        //      Nodal loads - - To build the vector F
        //**************************************************************************************************************
        getline (myReadFile, lineRead); // Nodal loads [kN]
        svCout << lineRead << endl;

        string *strNodLoad;
        strNodLoad = new string[4];

        nodalLoad = new double *[numberNodes];
        for (int i = 0; i < numberNodes; i++ ){
            nodalLoad[i] = new double[3];
        }
        for (int j = 0; j < numberNodes; ++j) {
            getline (myReadFile, lineRead);
            strNodLoad = splitString(lineRead, ' ');
            nodalLoad[stoi(strNodLoad[0])-1][0] = stod(strNodLoad[1]);
            nodalLoad[stoi(strNodLoad[0])-1][1] = stod(strNodLoad[2]);
            nodalLoad[stoi(strNodLoad[0])-1][2] = stod(strNodLoad[3]);

        }

        for (int i = 0; i < numberNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 3; ++j) {
                svCout << setw(12) << nodalLoad[i][j];
            }
            svCout << endl;
        }

        getline (myReadFile, lineRead);// Zoom Factor Deformability
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorD = stod(lineRead);
        svCout << zoomFactorD << endl;

        getline (myReadFile, lineRead);// Zoom Factor Normal Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorN = stod(lineRead);
        svCout << zoomFactorN << endl;

        getline (myReadFile, lineRead);// Zoom Factor Transverse Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorV = stod(lineRead);
        svCout << zoomFactorV << endl;

        getline (myReadFile, lineRead);// Zoom Factor Bending Moment Effort
        svCout << lineRead << endl;
        getline (myReadFile, lineRead);
        zoomFactorM = stod(lineRead);
        svCout << zoomFactorM << endl;

        svCout << "Reading completed successfully" << endl;
    }

    myReadFile.close();




//**********************************************************************************************************************
//     SOLUTION EQUATION: [K]{d} + {F0} = {F}
//**********************************************************************************************************************




    //******************************************************************************************************************
    // List of DOF - 3 for each node
    //******************************************************************************************************************
    degreesOfFreedom = new int *[numberNodes];
    for (int i = 0; i < numberNodes; i++ ){
        degreesOfFreedom[i] = new int[3];
    }
    for (int i = 0; i < numberNodes; i++ ){
        degreesOfFreedom[i][0] = 0;
        degreesOfFreedom[i][1] = 0;
        degreesOfFreedom[i][2] = 0;
    }
    int dof = 1;
    for(int i = 0; i < numberNodes; ++i){
        if(fixities[i][0] == 0){
            degreesOfFreedom[i][0] = dof;
            ++dof;
        }
        if(fixities[i][1] == 0){
            degreesOfFreedom[i][1] = dof;
            ++dof;
        }
        if(fixities[i][2] == 0){
            degreesOfFreedom[i][2] = dof;
            ++dof;
        }
    }

    svCout << "Degrees Of Freedom - DOF" << endl;
    for (int i = 0; i < numberNodes; ++i) {
        svCout << setw(5) << i+1 << setw(12) << degreesOfFreedom[i][0] << setw(12) << degreesOfFreedom[i][1] << setw(12) << degreesOfFreedom[i][2] <<endl;
    }


    //******************************************************************************************************************
    // Matrix reset to zero - Stiffness Matrix - [k] (size equal to dofXdof)
    //******************************************************************************************************************
    dof = degreesOfFreedom[numberNodes-1][2];
    int maxInDegreesOfFreedom = 0;
    for (int i = 0; i < numberNodes; ++i){
        for(int j = 0; j < 3; ++j){
            if(degreesOfFreedom[i][j] > maxInDegreesOfFreedom){
                maxInDegreesOfFreedom = degreesOfFreedom[i][j];
            }
        }
    }
    dof = maxInDegreesOfFreedom;
    svCout << "DOF " << endl;
    svCout << setw(12) << dof << endl;

    stiffnessMatrix = new double *[dof];
    for (int i = 0; i < dof; i++ ){
        stiffnessMatrix[i] = new double[dof];
    }
    for (int i = 0; i < dof; i++ ){
        for (int j = 0; j < dof; j++ ){
            stiffnessMatrix[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Vector reset to zero - Load Vector - Force reactions of the elements with distributed load - {F0}
    // (size equal to dof)
    //******************************************************************************************************************
    f0LoadVector = new double [dof];
    for (int i = 0; i < dof; i++ ){
        f0LoadVector[i] = 0.0;
    }




    //******************************************************************************************************************
    // Vector reset to zero - Load Vector - Force nodal - {F} (size equal to dof)
    //******************************************************************************************************************
    fLoadVector = new double [dof];
    for (int i = 0; i < dof; i++ ){
        fLoadVector[i] = 0.0;
    }




    //******************************************************************************************************************
    // Incidence Matrix - Mapping the local degrees of freedom with the global degrees of freedom
    // (size equal to numberElementsXdof)
    //******************************************************************************************************************
    incidence = new int *[numberElements];
    for (int i = 0; i < numberElements; i++ ){
        incidence[i] = new int[dof];
    }
    for (int i = 0; i < numberElements; i++ ){
        for (int j = 0; j < dof; j++ ){
            incidence[i][j] = 0;
        }
    }

    int nodeI;
    int nodeJ;
    int dof1, dof2, dof3, dof4, dof5, dof6;
    for(int i = 0; i < numberElements; i++){

        nodeI = elements[i][0];
        dof1 = degreesOfFreedom[nodeI-1][0];
        dof2 = degreesOfFreedom[nodeI-1][1];
        dof3 = degreesOfFreedom[nodeI-1][2];

        nodeJ = elements[i][1];
        dof4 = degreesOfFreedom[nodeJ-1][0];
        dof5 = degreesOfFreedom[nodeJ-1][1];
        dof6 = degreesOfFreedom[nodeJ-1][2];

        if( dof1 != 0) incidence[i][dof1-1] = 1;
        if( dof2 != 0) incidence[i][dof2-1] = 2;
        if( dof3 != 0) incidence[i][dof3-1] = 3;
        if( dof4 != 0) incidence[i][dof4-1] = 4;
        if( dof5 != 0) incidence[i][dof5-1] = 5;
        if( dof6 != 0) incidence[i][dof6-1] = 6;
    }

    svCout << "Incidence Matrix " << endl;
    for (int i = 0; i < numberElements; ++i) {
        svCout << setw(5) << i+1;
        for (int j = 0; j < dof; ++j) {
            svCout << setw(12) << incidence[i][j];
        }
        svCout << endl;
    }




    //******************************************************************************************************************
    // Assembly Stiffness Matrix (add all [k_G_i] to [k])
    // and assembly Load Vector (add all {F0_G_i} to {F0})
    // taking into account the incidence matrix
    // Where:
    // [k] as said before, is the Stiffness Matrix of the structure
    // {F0} as said before, is the Load Vector of the structure
    // [k_G_i] is the Stiffness Global Matrix of the element i and is equal to:
    //         [K_G_i] = [t_i]^T[k_L_i][t_i]
    // [k_L_i] is the Stiffness Local Matrix of the element i
    // [t] is the transformation matrix of the element i
    // [t]^T is the transpose of the [t]
    // {F0_G_i} is the Local Load Vector fo the element i and is equal to:
    //         {F0_G_i} = [t_i]^T{F0_L_i}
    // {F0_L_i} is the Local Vector Load of the element i
    //******************************************************************************************************************
    double x0, y0, x1, y1, length, theta;
    double k1, k2, k3, k4, k5;
    double PI = 4.0*atan(1.0);




    //******************************************************************************************************************
    // Matrix reset to zero -  Local Stiffness Matrix
    //******************************************************************************************************************
    localStiffnessMatrix = new double *[6];
    for (int i = 0; i < 6; i++ ){
        localStiffnessMatrix[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            localStiffnessMatrix[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Matrix reset to zero -  Global Stiffness Matrix
    //******************************************************************************************************************
    globalStiffnessMatrix = new double *[6];
    for (int i = 0; i < 6; i++ ){
        globalStiffnessMatrix[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            globalStiffnessMatrix[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Vector reset to zero -  Local F0 Load Vector
    //******************************************************************************************************************
    localF0LoadVector = new double [6];
    for (int i = 0; i < 6; i++ ){
        localF0LoadVector[i] = 0.0;
    }




    //******************************************************************************************************************
    // Vector reset to zero -  Global F0 Load Vector
    //******************************************************************************************************************
    globalF0LoadVector = new double [6];
    for (int i = 0; i < 6; i++ ){
        globalF0LoadVector[i] = 0.0;
    }




    //******************************************************************************************************************
    // Matrix reset to zero -  Transformation
    //******************************************************************************************************************
    transformation = new double *[6];
    for (int i = 0; i < 6; i++ ){
        transformation[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            transformation[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Matrix reset to zero -  Transpose Transformation
    //******************************************************************************************************************
    transformationTranspose = new double *[6];
    for (int i = 0; i < 6; i++ ){
        transformationTranspose[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            transformationTranspose[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Matrix reset to zero -  Auxiliary Matrix - [aux] = [k_L_i][t_i]
    //******************************************************************************************************************
    double **aux;
    aux = new double *[6];
    for (int i = 0; i < 6; i++ ){
        aux[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            aux[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Stiffness Matrix Assembling and Load Matrix Assembling
    // Loop through all elements i and make the Stiffness Matrix and the Load Vector of the structure
    //******************************************************************************************************************
    int linekElement;
    int columnkElement;
    for(int i = 0; i < numberElements; i++){

        svCout << "ELEMENT " << setw(5) << i+1 << endl;

        x0 = nodes[elements[i][0]-1][0];
        y0 = nodes[elements[i][0]-1][1];
        x1 = nodes[elements[i][1]-1][0];
        y1 = nodes[elements[i][1]-1][1];

        length = sqrt(pow(x1-x0,2)+pow(y1-y0,2)); //cout << length << endl;
        theta = abs(x1-x0) <= 1/10000.0 ? PI/2: atan((y1-y0)/(x1-x0)); //cout << theta << endl;

        svCout << "Length " << endl;
        svCout << setw(12) << length << endl;
        svCout << "Theta [degrees]" << endl;
        svCout << setw(12) << theta*180.0/PI << endl;

        k1 = modE*materialsArea[i]/length;
        k2 = 12*modE*materialsInertia[i]/(pow(length,3));
        k3 =  6*modE*materialsInertia[i]/(pow(length,2));
        k4 =  4*modE*materialsInertia[i]/length;
        k5 = k4/2;

        localStiffnessMatrix[0][0] = k1;
        localStiffnessMatrix[0][1] = 0;
        localStiffnessMatrix[0][2] = 0;
        localStiffnessMatrix[0][3] = -k1;
        localStiffnessMatrix[0][4] = 0;
        localStiffnessMatrix[0][5] = 0;

        localStiffnessMatrix[1][0] = 0;
        localStiffnessMatrix[1][1] = k2;
        localStiffnessMatrix[1][2] = k3;
        localStiffnessMatrix[1][3] = 0;
        localStiffnessMatrix[1][4] = -k2;
        localStiffnessMatrix[1][5] = k3;

        localStiffnessMatrix[2][0] = 0;
        localStiffnessMatrix[2][1] = k3;
        localStiffnessMatrix[2][2] = k4;
        localStiffnessMatrix[2][3] = 0;
        localStiffnessMatrix[2][4] = -k3;
        localStiffnessMatrix[2][5] = k5;

        localStiffnessMatrix[3][0] = -k1;
        localStiffnessMatrix[3][1] = 0;
        localStiffnessMatrix[3][2] = 0;
        localStiffnessMatrix[3][3] = k1;
        localStiffnessMatrix[3][4] = 0;
        localStiffnessMatrix[3][5] = 0;

        localStiffnessMatrix[4][0] = 0;
        localStiffnessMatrix[4][1] = -k2;
        localStiffnessMatrix[4][2] = -k3;
        localStiffnessMatrix[4][3] = 0;
        localStiffnessMatrix[4][4] = k2;
        localStiffnessMatrix[4][5] = -k3;

        localStiffnessMatrix[5][0] = 0;
        localStiffnessMatrix[5][1] = k3;
        localStiffnessMatrix[5][2] = k5;
        localStiffnessMatrix[5][3] = 0;
        localStiffnessMatrix[5][4] = -k3;
        localStiffnessMatrix[5][5] = k4;

        svCout << "Local Stiffness Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << localStiffnessMatrix[i][j];
            }
            svCout << endl;
        }

        localF0LoadVector[0] = 0.0;
        localF0LoadVector[1] = distLoad[i]*length/2;
        localF0LoadVector[2] = distLoad[i]*pow(length,2)/12;
        localF0LoadVector[3] = 0.0;
        localF0LoadVector[4] = distLoad[i]*length/2;
        localF0LoadVector[5] = -distLoad[i]*pow(length,2)/12;

        svCout << "Local F0 Load" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(12) << localF0LoadVector[i] << endl;
        }

        for (int i = 0; i < 6; i++ ){
            for (int j = 0; j < 6; j++ ){
                transformation[i][j] = 0.0;
            }
        }
        transformation[0][0] = cos(theta);
        transformation[0][1] = sin(theta);
        transformation[1][0] = -sin(theta);
        transformation[1][1] = cos(theta);
        transformation[2][2] = 1;
        transformation[3][3] = cos(theta);
        transformation[3][4] = sin(theta);
        transformation[4][3] = -sin(theta);
        transformation[4][4] = cos(theta);
        transformation[5][5] = 1;

        svCout << "Transformation Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << transformation[i][j];
            }
            svCout << endl;
        }

        for (int i = 0; i < 6; i++ ){
            for (int j = 0; j < 6; j++ ){
                transformationTranspose[i][j] = 0.0;
            }
        }
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                transformationTranspose[i][j] = transformation[j][i];
            }
        }

        svCout << "Transpose Transformation Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << transformationTranspose[i][j];
            }
            svCout << endl;
        }

        //Set [globalStiffness] = [transformation]^T * [localStiffness] * [transformation]
        //                                        ij jk ij jk ij jk  ij jk ij jk ij jk  ij jk ij jk ij jk
        // | 11 12 13 | x | 11 12 13 | = [ik] = | 11x11+12x21+13x31  11x12+12x22+13x32  11x13+12x23+13x33 |
        // | 21 22 23 |   | 21 22 23 |          | 21x11+22x21+23x31  21x12+22x22+23x32  21x13+22x23+23x33 |
        // | 31 32 33 |   | 31 32 33 |          | 31x11+32x21+33x31  31x12+32x22+33x32  31x13+32x23+33x33 |
        //Auxiliary matrix -  [aux] = [localStiffness][transformation]
        for (int i = 0; i < 6; i++) {
            for (int k = 0; k < 6; k++) {
                aux[i][k] = 0;
                for (int j = 0; j < 6; j++) {
                    aux[i][k] += localStiffnessMatrix[i][j] * transformation[j][k];
                }
            }
        }
        for (int i = 0; i < 6; i++) {
            for (int k = 0; k < 6; k++) {
                globalStiffnessMatrix[i][k] = 0;
                for (int j = 0; j < 6; j++) {
                    globalStiffnessMatrix[i][k] = globalStiffnessMatrix[i][k] + transformationTranspose[i][j] * aux[j][k];
                }
            }
        }
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                if (abs(globalStiffnessMatrix[i][j]) <= 1.0/10000.0)
                    globalStiffnessMatrix[i][j] = 0;
            }
        }

        svCout << "Global Stiffness Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << globalStiffnessMatrix[i][j];
            }
            svCout << endl;
        }

        //Set globalF0LoadVector - globalF0LoadVector
        // | 11 12 13 | x | 1 | = | 11x1 + 12x2 + 13x3 |
        // | 21 22 23 |   | 2 |   | 21x1 + 22x2 + 23x3 |
        // | 31 32 33 |   | 3 |   | 31x1 + 32x2 + 33x3 |
        for (int i = 0; i < 6; i++) {
            globalF0LoadVector[i] = 0;
            for (int j = 0; j < 6; j++) {
                globalF0LoadVector[i] += transformationTranspose[i][j] * localF0LoadVector[j];
            }
        }

        svCout << "Global F0 Load" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(12) << globalF0LoadVector[i] << endl;
        }

        //Assembly Stiffness Matrix: j - line | k - column
        linekElement = 0;
        columnkElement = 0;
        for(int j = 0; j < dof; j++){
            if(incidence[i][j] != 0){
                linekElement = incidence[i][j]-1;
                for(int k = 0; k < dof; k++){
                    if(incidence[i][k] != 0){
                        columnkElement = incidence[i][k]-1;
                        stiffnessMatrix[j][k] = stiffnessMatrix[j][k] + globalStiffnessMatrix[linekElement][columnkElement];
                    }
                }
            }
        }

        //Assembly F0 Load
        linekElement = 0;
        for(int j = 0; j < dof; j++){
            if(incidence[i][j] != 0){
                linekElement = incidence[i][j]-1;
                f0LoadVector[j] = f0LoadVector[j] + globalF0LoadVector[linekElement];
            }
        }

    } // END for(int i = 0; i < numberElements; i++)

    //Map F forces to DOF
    for (int i = 0; i < numberNodes; ++i) {
        if(degreesOfFreedom[i][0] != 0){
            fLoadVector[degreesOfFreedom[i][0]-1] = nodalLoad[i][0];
        }
        if(degreesOfFreedom[i][1] != 0){
            fLoadVector[degreesOfFreedom[i][1]-1] = nodalLoad[i][1];
        }
        if(degreesOfFreedom[i][2] != 0){
            fLoadVector[degreesOfFreedom[i][2]-1] = nodalLoad[i][2];
        }
    }

    svCout << "Stiffness Matrix (A)" << endl;
    for (int i = 0; i < dof; ++i) {
        svCout << setw(5) << i+1;
        for (int j = 0; j < dof; ++j) {
            svCout << setw(12) << stiffnessMatrix[i][j];
        }
        svCout << endl;
    }

    svCout << "F0 Load Matrix" << endl;
    for (int i = 0; i < dof; ++i) {
        svCout << setw(5) << i+1 << setw(12) << f0LoadVector[i] << endl;
    }

    svCout << "F Load Matrix" << endl;
    for (int i = 0; i < dof; ++i) {
        svCout << setw(5) << i+1 << setw(12) << fLoadVector[i] << endl;
    }




    //******************************************************************************************************************
    // Solution system:  K.x = b   with b = F - F0
    //******************************************************************************************************************
    double *fMinusf0;
    fMinusf0 = new double [dof];
    for (int i = 0; i < dof; i++ ){
        fMinusf0[i] = 0.0;
        fMinusf0[i] = fLoadVector[i] - f0LoadVector[i];
    }
    svCout << "Vector B (B = F - F0) " << endl;
    for (int i = 0; i < dof; ++i) {
        svCout << setw(5) << i+1 << setw(12) << fMinusf0[i] << endl;
    }

    double *solutionSystem;
    solutionSystem = new double [dof];
    for (int i = 0; i < dof; i++ ){
        solutionSystem[i] = 0.0;
    }

    solutionSystem = r8rmat_fs_new ( dof, stiffnessMatrix, fMinusf0);

    svCout << "Solution System (AX = B)" << endl;
    for (int i = 0; i < dof; ++i) {
        svCout << setw(5) << i+1 << setw(16) << solutionSystem[i] << endl;
    }




    //******************************************************************************************************************
    // Auxiliary matrix to follow the material strength convention
    //******************************************************************************************************************
    double **strengthMaterialMatrix; //
    strengthMaterialMatrix = new double *[6];
    for (int i = 0; i < 6; i++ ){
        strengthMaterialMatrix[i] = new double[6];
    }
    for (int i = 0; i < 6; i++ ){
        for (int j = 0; j < 6; j++ ){
            strengthMaterialMatrix[i][j] = 0;
        }
    }
    strengthMaterialMatrix[0][0] = -1;
    strengthMaterialMatrix[1][1] = 1;
    strengthMaterialMatrix[2][2] = -1;
    strengthMaterialMatrix[3][3] = 1;
    strengthMaterialMatrix[4][4] = -1;
    strengthMaterialMatrix[5][5] = 1;

    svCout << "Strength Material Matrix" << endl;
    for (int i = 0; i < 6; ++i) {
        svCout << setw(5) << i+1;
        for (int j = 0; j < 6; ++j) {
            svCout << setw(12) << strengthMaterialMatrix[i][j];
        }
        svCout << endl;
    }


    //******************************************************************************************************************
    // Reset all matrix:
    // Back to top: From Global Displacement Matrix to Local Displacement Matrix
    //******************************************************************************************************************
    double *globalDisplMatrix;
    globalDisplMatrix = new double [6];
    for (int i = 0; i < 6; i++ ){
        globalDisplMatrix[i] = 0.0;
    }
    double *localDisplMatrix;
    localDisplMatrix = new double [6];
    for (int i = 0; i < 6; i++ ){
        localDisplMatrix[i] = 0.0;
    }
    double *vectorEd;
    vectorEd = new double [6];
    for (int i = 0; i < 6; i++ ){
        vectorEd[i] = 0.0;
    }
    double *vectorEdStrengthMaterial;
    vectorEdStrengthMaterial = new double [6];
    for (int i = 0; i < 6; i++ ){
        vectorEdStrengthMaterial[i] = 0.0;
    }
    double *zeros;
    zeros = new double[2];
    zeros[0] = 0; zeros[1] = 0;

    double **transf2DMatrix;
    int it2DM = 4;
    transf2DMatrix = new double *[it2DM];
    for (int i = 0; i < it2DM; i++ ){
        transf2DMatrix[i] = new double[it2DM];
    }
    for (int i = 0; i < it2DM; i++ ){
        for (int j = 0; j < it2DM; j++ ){
            transf2DMatrix[i][j] = 0.0;
        }
    }




    //******************************************************************************************************************
    // Deformability and Diagrams
    //******************************************************************************************************************
    int numberSnippet = 10; // Splitting the length of the element into several sections
    int nLines = 0;
    int nNodes = 0;
    double xDef = 0.0, yDef = 0.0; // Coordinates of a point of the final local deformed
    double deltaDef1 = 0.0; // The initial displacement at x = 0 (Displacement constant along the element)
    double deltaDef2 = 0.0; // The displacement due to the difference in extremity displacements - x = 0 and x = L
    double deltaDef3 = 0.0; // The displacement due to the deformation of the element only
    double deltaX = 0.0; // Increment of the x coordinate to calculate the next point (xDef, yDef)
    double xi = 0.0; // The xi parameter for the function deltaDeformationFunction() - Value between 0.0 and 1.0  (xi = xdef/length)

    svLatexDefCout << "path listMPDef[];" << endl;
    svLatexMomCout << "path listMPMom[];" << endl;
    svLatexSheCout << "path listMPShe[];" << endl;
    svLatexNorCout << "path listMPNor[];" << endl;
    svLatexMomCout << "path mpMom[];" << endl;
    svLatexSheCout << "path mpShe[];" << endl;
    svLatexNorCout << "path mpNor[];" << endl;
    int kMom = 0;
    int kShe = 0;
    int kNor = 0;
    for(int i = 0; i < numberElements; i++){

        svCout << "ELEMENT " << setw(5) << i+1 << endl;

        x0 = nodes[elements[i][0]-1][0];
        y0 = nodes[elements[i][0]-1][1];
        x1 = nodes[elements[i][1]-1][0];
        y1 = nodes[elements[i][1]-1][1];

        length = sqrt(pow(x1-x0,2)+pow(y1-y0,2)); //cout << length << endl;
        theta = abs(x1-x0) <= 1.0/10000.0 ? PI/2: atan((y1-y0)/(x1-x0)); //cout << theta << endl;

        svCout << "Length " << endl;
        svCout << setw(12) << length << endl;
        svCout << "Theta [degrees]" << endl;
        svCout << setw(12) << theta*180.0/PI << endl;

        k1 = modE*materialsArea[i]/length;
        k2 = 12*modE*materialsInertia[i]/(pow(length,3));
        k3 =  6*modE*materialsInertia[i]/(pow(length,2));
        k4 =  4*modE*materialsInertia[i]/length;
        k5 = k4/2;

        localStiffnessMatrix[0][0] = k1;
        localStiffnessMatrix[0][1] = 0;
        localStiffnessMatrix[0][2] = 0;
        localStiffnessMatrix[0][3] = -k1;
        localStiffnessMatrix[0][4] = 0;
        localStiffnessMatrix[0][5] = 0;

        localStiffnessMatrix[1][0] = 0;
        localStiffnessMatrix[1][1] = k2;
        localStiffnessMatrix[1][2] = k3;
        localStiffnessMatrix[1][3] = 0;
        localStiffnessMatrix[1][4] = -k2;
        localStiffnessMatrix[1][5] = k3;

        localStiffnessMatrix[2][0] = 0;
        localStiffnessMatrix[2][1] = k3;
        localStiffnessMatrix[2][2] = k4;
        localStiffnessMatrix[2][3] = 0;
        localStiffnessMatrix[2][4] = -k3;
        localStiffnessMatrix[2][5] = k5;

        localStiffnessMatrix[3][0] = -k1;
        localStiffnessMatrix[3][1] = 0;
        localStiffnessMatrix[3][2] = 0;
        localStiffnessMatrix[3][3] = k1;
        localStiffnessMatrix[3][4] = 0;
        localStiffnessMatrix[3][5] = 0;

        localStiffnessMatrix[4][0] = 0;
        localStiffnessMatrix[4][1] = -k2;
        localStiffnessMatrix[4][2] = -k3;
        localStiffnessMatrix[4][3] = 0;
        localStiffnessMatrix[4][4] = k2;
        localStiffnessMatrix[4][5] = -k3;

        localStiffnessMatrix[5][0] = 0;
        localStiffnessMatrix[5][1] = k3;
        localStiffnessMatrix[5][2] = k5;
        localStiffnessMatrix[5][3] = 0;
        localStiffnessMatrix[5][4] = -k3;
        localStiffnessMatrix[5][5] = k4;

        svCout << "Local Stiffness Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << localStiffnessMatrix[i][j];
            }
            svCout << endl;
        }

        localF0LoadVector[0] = 0.0;
        localF0LoadVector[1] = distLoad[i]*length/2;
        localF0LoadVector[2] = distLoad[i]*pow(length,2)/12;
        localF0LoadVector[3] = 0.0;
        localF0LoadVector[4] = distLoad[i]*length/2;
        localF0LoadVector[5] = -distLoad[i]*pow(length,2)/12;

        svCout << "Local Load" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(12) << localF0LoadVector[i] << endl;
        }

        for (int i = 0; i < 6; i++ ){
            for (int j = 0; j < 6; j++ ){
                transformation[i][j] = 0.0;
            }
        }
        transformation[0][0] = cos(theta);
        transformation[0][1] = sin(theta);
        transformation[1][0] = -sin(theta);
        transformation[1][1] = cos(theta);
        transformation[2][2] = 1;
        transformation[3][3] = cos(theta);
        transformation[3][4] = sin(theta);
        transformation[4][3] = -sin(theta);
        transformation[4][4] = cos(theta);
        transformation[5][5] = 1;

        svCout << "Transformation" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 6; ++j) {
                svCout << setw(12) << transformation[i][j];
            }
            svCout << endl;
        }

        for (int i = 0; i < 6; i++ ){
            globalDisplMatrix[i] = 0.0;
        }
        int iGlobElem;
        for (int k = 0; k < dof; k++){
            if(incidence[i][k] != 0){
                iGlobElem = incidence[i][k]-1;
                globalDisplMatrix[iGlobElem] = solutionSystem[k];
            }
        }

        svCout << "Global Displacement Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(16) << globalDisplMatrix[i] << endl;
        }

        for (int i = 0; i < 6; i++) {
            localDisplMatrix[i] = 0;
            for (int j = 0; j < 6; j++) {
                localDisplMatrix[i] += transformation[i][j] * globalDisplMatrix[j];
            }
        }

        svCout << "Local Displacement Matrix" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(16) << localDisplMatrix[i] << endl;
        }

        for (int i = 0; i < 6; i++) {
            vectorEd[i] = localF0LoadVector[i];
            for (int j = 0; j < 6; j++) {
                vectorEd[i] += localStiffnessMatrix[i][j] * localDisplMatrix[j];
            }
        }
        svCout << "Vector Ed" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(12) << vectorEd[i] << endl;
        }

        for (int i = 0; i < 6; i++) {
            vectorEdStrengthMaterial[i] = 0;
            for (int j = 0; j < 6; j++) {
                vectorEdStrengthMaterial[i] += strengthMaterialMatrix[i][j] * vectorEd[j];
            }
        }

        svCout << "Vector Ed - Strength Materials" << endl;
        for (int i = 0; i < 6; ++i) {
            svCout << setw(5) << i+1 << setw(12) << vectorEdStrengthMaterial[i] << endl;
        }

        zeros = zeroInMEd(distLoad[i], vectorEdStrengthMaterial[1], vectorEdStrengthMaterial[2], length);
        svCout << "MEd equal to zero at" << endl;
        if(zeros[0] != -1){
            svCout << setw(12) << zeros[0] <<endl;
        }else{
            svCout << "         ---" << endl;
        }
        if(zeros[1] != -1){
            svCout << setw(12) << zeros[1] <<endl;
        }else{
            svCout << "         ---" << endl;
        }

        double xCancelsVEd = zeroInVEd(distLoad[i], vectorEdStrengthMaterial[1], length);
        double MEdMax;
        if (xCancelsVEd != -1){
            MEdMax = equationMEd(distLoad[i], vectorEdStrengthMaterial[1], vectorEdStrengthMaterial[2], xCancelsVEd);
            svCout << "VEd equal to zero at" << endl;
            svCout << setw(12) << xCancelsVEd << endl;
            svCout << "Maximum MEd - in the element" << endl;
            svCout << setw(12) << MEdMax << endl;
        } else {
            svCout << "VEd equal to zero at" << endl;
            svCout << "         ---" << endl;
            svCout << "Maximum MEd - at one end of the element" << endl;
            MEdMax = abs(vectorEdStrengthMaterial[2])>=abs(vectorEdStrengthMaterial[5]) ? vectorEdStrengthMaterial[2]:vectorEdStrengthMaterial[5];
            svCout << setw(12) << MEdMax << endl;
        }




        //**************************************************************************************************************
        //          Deformability
        //**************************************************************************************************************
        /*
         * double zoomFactor = 1.0; // Displacement amplification
         * int numberSnippet = 1; // Splitting the length of the element into several sections
         * int nLines = 0;
         * int nNodes = 0;
         * double xDef = 0.0, yDef = 0.0; // Coordinates of a point of the final local deformed
         * double deltaDef1 = 0.0; // The initial displacement at x = 0 (Displacement constant along the element)
         * double deltaDef2 = 0.0; // The displacement due to the difference in extremity displacements - x = 0 and x = L
         * double deltaDef3 = 0.0; // The displacement due to the deformation of the element only
         * double deltaX = 0.0; // Increment of the x coordinate to calculate the next point (xDef, yDef)
         * double xi = 0.0; // The xi parameter for the function deltaDeformationFunction()
         */
        transf2DMatrix[0][0] = cos(theta);  transf2DMatrix[0][1] = sin(theta);
        transf2DMatrix[1][0] = -sin(theta); transf2DMatrix[1][1] = cos(theta);
        transf2DMatrix[2][2] = cos(theta);  transf2DMatrix[2][3] = sin(theta);
        transf2DMatrix[3][2] = -sin(theta); transf2DMatrix[3][3] = cos(theta);

        // 1    2    3     4     5
        // *----*----*-----*-----*
        nLines = numberSnippet; // Element i subdivided into nLines lines
        nNodes = nLines + 1; // Element i subdivided into nLines lines has nNodes nodes

        svCout << "Deformability" << endl;
        svCout << "Lenght" << endl;
        svCout << setw(5) << length << endl;
        svCout << "Snippet - Number" << endl;
        svCout << setw(5) << numberSnippet << endl;
        svCout << "Nodes in the division - Number" << endl;
        svCout << setw(5) << nNodes << endl;

        // Coordinate local of point 1 (x0, y0) and point 2 (x1, y1) of section m of element i deformed
        // [nLineElementDeformability] = | x0 y0 |
        //                               |  ...  |
        double **nLineElementDeformability; // In local coordinates
        nLineElementDeformability = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            nLineElementDeformability[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                nLineElementDeformability[i][j] = 0.0;
            }
        }

        // Coordinate global of point 1 (x0, y0) and point 2 (x1, y1) of section m of element i deformed
        double **nLineElementDeformabilityGlobal; // In global coordinates
        nLineElementDeformabilityGlobal = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            nLineElementDeformabilityGlobal[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                nLineElementDeformabilityGlobal[i][j] = 0.0;
            }
        }

        svCout << "Zoom factor - Deformability" << endl;
        svCout << setw(5) << zoomFactorD << endl;
        deltaX = (double) length/nLines;
        svCout << "Delta x" << endl;
        svCout << setw(5) << deltaX << endl;

        xi = 0.0;

        xDef = 0;

        deltaDef1 = localDisplMatrix[1];
        deltaDef2 = (localDisplMatrix[4] - localDisplMatrix[1])*xi;
        deltaDef3 = deltaDeforma(length, modE, materialsInertia[i],vectorEdStrengthMaterial[2],vectorEdStrengthMaterial[5],xi);
        yDef = zoomFactorD*(deltaDef1 + deltaDef2 + deltaDef3);

        // Local
        nLineElementDeformability[0][0] = xDef + zoomFactorD*localDisplMatrix[0];
        nLineElementDeformability[0][1] = yDef;

        // Global
        nLineElementDeformabilityGlobal[0][0] = x0 + transf2DMatrix[0][0] * nLineElementDeformability[0][0] + transf2DMatrix[1][0] * nLineElementDeformability[0][1];
        nLineElementDeformabilityGlobal[0][1] = y0 + transf2DMatrix[0][1] * nLineElementDeformability[0][0] + transf2DMatrix[1][1] * nLineElementDeformability[0][1];

        for (int m = 1; m < nNodes; m++){

            xDef += deltaX;
            xi = xDef/length; // greek letter csi
            deltaDef1 = localDisplMatrix[1];
            deltaDef2 = (localDisplMatrix[4] - localDisplMatrix[1])*xi;
            deltaDef3 = deltaDeforma(length, modE, materialsInertia[i],vectorEdStrengthMaterial[2],vectorEdStrengthMaterial[5],xi);
            yDef = zoomFactorD*(deltaDef1 + deltaDef2 + deltaDef3);

            // Local
            nLineElementDeformability[m][0] = xDef+ zoomFactorD*localDisplMatrix[0];
            nLineElementDeformability[m][1] = yDef;

            //Global
            nLineElementDeformabilityGlobal[m][0] = x0 + transf2DMatrix[0][0] * nLineElementDeformability[m][0] + transf2DMatrix[1][0] * nLineElementDeformability[m][1];
            nLineElementDeformabilityGlobal[m][1] = y0 + transf2DMatrix[0][1] * nLineElementDeformability[m][0] + transf2DMatrix[1][1] * nLineElementDeformability[m][1];

        }

        svCout << "Local - Line Element Deformability" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << nLineElementDeformability[i][j];
            }
            svCout << endl;
        }
        svCout << "Global - Line Element Deformability" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << nLineElementDeformabilityGlobal[i][j];
            }
            svCout << endl;
        }


        string listx = "defElem" + to_string(i+1) + "x = [";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listx += to_string(nLineElementDeformabilityGlobal[i][0]) + ", ";
            }else{
                listx += to_string(nLineElementDeformabilityGlobal[i][0]) + "]";
            }
        }
        svPythonDefCout << listx << endl;

        // For use in python
        string listy = "defElem" + to_string(i+1) + "y = [";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listy += to_string(nLineElementDeformabilityGlobal[i][1]) + ", ";
            }else{
                listy += to_string(nLineElementDeformabilityGlobal[i][1]) + "]";
            }
        }
        svPythonDefCout << listy << endl;



        // For use in MetaPost in LaTex
        string listMPDef = "listMPDef[" + to_string(i) + "] = ";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listMPDef += "(" + to_string(nLineElementDeformabilityGlobal[i][0]) + "*scale+ offsetX, " + to_string(nLineElementDeformabilityGlobal[i][1]) + "*scale+ offsetY)--";
            }else{
                listMPDef += "(" + to_string(nLineElementDeformabilityGlobal[i][0]) + "*scale+ offsetX, " + to_string(nLineElementDeformabilityGlobal[i][1]) + "*scale+ offsetY);";
            }
        }
        svLatexDefCout << listMPDef << endl;


        //**************************************************************************************************************
        //          Effort diagrams - M
        //**************************************************************************************************************
        nLines = numberSnippet;
        nNodes = nLines + 1;
        deltaX = (double) length/nLines;

        svCout << "Effort diagram - M" << endl;
        svCout << "Lenght" << endl;
        svCout << setw(5) << length << endl;
        svCout << "Snippet - Number" << endl;
        svCout << setw(5) << numberSnippet << endl;
        svCout << "Nodes in the division - Number" << endl;
        svCout << setw(5) << nNodes << endl;
        svCout << "Zoom factor - M" << endl;
        svCout << setw(5) << zoomFactorM << endl;
        svCout << "Delta x" << endl;
        svCout << setw(5) << deltaX << endl;

        double **localDiagM; // In local coordinates
        localDiagM = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            localDiagM[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                localDiagM[i][j] = 0.0;
            }
        }
        double **globalDiagM; // In global coordinates
        globalDiagM = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalDiagM[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalDiagM[i][j] = 0.0;
            }
        }
        double **globalAxisDiagM; // In global coordinates
        globalAxisDiagM = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalAxisDiagM[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalAxisDiagM[i][j] = 0.0;
            }
        }

        // Initialization of the diagrams calculation
        xDef = 0.0;
        for (int m = 0; m < nNodes; m++){
            // Local
            localDiagM[m][0] = xDef;
            // (-1.0) -  Because of the convention of moments. The direction of the ordinate axis is switched.
            localDiagM[m][1] = zoomFactorM*(-1.0)*equationMEd(distLoad[i], vectorEdStrengthMaterial[1], vectorEdStrengthMaterial[2], xDef);
            xDef += deltaX;

            // Global
            globalDiagM[m][0] = x0 + transf2DMatrix[0][0] * localDiagM[m][0] + transf2DMatrix[1][0] * localDiagM[m][1];
            globalDiagM[m][1] = y0 + transf2DMatrix[0][1] * localDiagM[m][0] + transf2DMatrix[1][1] * localDiagM[m][1];

            // Global Axis
            globalAxisDiagM[m][0] = x0 + transf2DMatrix[0][0] * localDiagM[m][0];
            globalAxisDiagM[m][1] = y0 + transf2DMatrix[0][1] * localDiagM[m][0];
        }

        svCout << "Local - Diagram M" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << localDiagM[i][j];
            }
            svCout << endl;
        }

        svCout << "Global - Diagram M" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalDiagM[i][j];
            }
            svCout << endl;
        }

        svCout << "Global Axis - Diagram M" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalAxisDiagM[i][j];
            }
            svCout << endl;
        }

        // For use in MetaPost in LaTex
        string listMPMom = "listMPMom[" + to_string(i) + "] = ";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listMPMom += "(" + to_string(globalDiagM[i][0]) + "*scale+ offsetX, " + to_string(globalDiagM[i][1]) + "*scale+ offsetY)--";
            }else{
                listMPMom += "(" + to_string(globalDiagM[i][0]) + "*scale+ offsetX, " + to_string(globalDiagM[i][1]) + "*scale+ offsetY);";
            }
        }
        svLatexMomCout << listMPMom << endl;
        svLatexMomCout << "draw listMPMom[" + to_string(i) + "] withcolor red;" << endl;

        // For use in MetaPost in LaTex
        for (int i = 0; i < nNodes; ++i) {
            svLatexMomCout << "mpMom[" + to_string(kMom) + "] = (" + to_string(globalAxisDiagM[i][0]) + "*scale+ offsetX, " + to_string(globalAxisDiagM[i][1]) + "*scale+ offsetY)--(" + to_string(globalDiagM[i][0]) + "*scale+ offsetX, " + to_string(globalDiagM[i][1]) + "*scale+ offsetY);" << endl;
            svLatexMomCout << "draw mpMom["+ to_string(kMom)+"] withcolor red;" << endl;
            ++kMom;
        }


        //**************************************************************************************************************
        //          Effort diagrams - V
        //**************************************************************************************************************
        nLines = numberSnippet ;
        nNodes = nLines + 1;
        deltaX = (double) length/nLines;

        svCout << "Effort diagram - V" << endl;
        svCout << "Lenght" << endl;
        svCout << setw(5) << length << endl;
        svCout << "Snippet - Number" << endl;
        svCout << setw(5) << numberSnippet << endl;
        svCout << "Nodes in the division - Number" << endl;
        svCout << setw(5) << nNodes << endl;
        svCout << "Zoom factor - V" << endl;
        svCout << setw(5) << zoomFactorV << endl;
        svCout << "Delta x" << endl;
        svCout << setw(5) << deltaX << endl;

        double **localDiagV; // In local coordinates
        localDiagV = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            localDiagV[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                localDiagV[i][j] = 0.0;
            }
        }
        double **globalDiagV; // In global coordinates
        globalDiagV = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalDiagV[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalDiagV[i][j] = 0.0;
            }
        }
        double **globalAxisDiagV; // In global coordinates
        globalAxisDiagV = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalAxisDiagV[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalAxisDiagV[i][j] = 0.0;
            }
        }

        // Initialization of the diagrams calculation
        xDef = 0.0;
        for (int m = 0; m < nNodes; m++){
            // Local
            localDiagV[m][0] = xDef;
            localDiagV[m][1] = zoomFactorV*equationVEd(distLoad[i], vectorEdStrengthMaterial[1], xDef);
            xDef += deltaX;

            // Global
            globalDiagV[m][0] = x0 + transf2DMatrix[0][0] * localDiagV[m][0] + transf2DMatrix[1][0] * localDiagV[m][1];
            globalDiagV[m][1] = y0 + transf2DMatrix[0][1] * localDiagV[m][0] + transf2DMatrix[1][1] * localDiagV[m][1];

            // Global Axis
            globalAxisDiagV[m][0] = x0 + transf2DMatrix[0][0] * localDiagV[m][0];
            globalAxisDiagV[m][1] = y0 + transf2DMatrix[0][1] * localDiagV[m][0];
        }

        svCout << "Local - Diagram V" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << localDiagV[i][j];
            }
            svCout << endl;
        }

        svCout << "Global - Diagram V" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalDiagV[i][j];
            }
            svCout << endl;
        }

        svCout << "Global Axis - Diagram V" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalAxisDiagV[i][j];
            }
            svCout << endl;
        }



        // For use in MetaPost in LaTex
        string listMPShe = "listMPShe[" + to_string(i) + "] = ";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listMPShe += "(" + to_string(globalDiagV[i][0]) + "*scale+ offsetX, " + to_string(globalDiagV[i][1]) + "*scale+ offsetY)--";
            }else{
                listMPShe += "(" + to_string(globalDiagV[i][0]) + "*scale+ offsetX, " + to_string(globalDiagV[i][1]) + "*scale+ offsetY);";
            }
        }
        svLatexSheCout << listMPShe << endl;
        svLatexSheCout << "draw listMPShe[" + to_string(i) + "] withcolor green;" << endl;


        // For use in MetaPost in LaTex
        for (int i = 0; i < nNodes; ++i) {
            svLatexSheCout << "mpShe[" + to_string(kShe) + "] = (" + to_string(globalAxisDiagV[i][0]) + "*scale+ offsetX, " + to_string(globalAxisDiagV[i][1]) + "*scale+ offsetY)--(" + to_string(globalDiagV[i][0]) + "*scale+ offsetX, " + to_string(globalDiagV[i][1]) + "*scale+ offsetY);" << endl;
            svLatexSheCout << "draw mpShe["+ to_string(kShe)+"] withcolor green;" << endl;
            ++kShe;
        }

        //**************************************************************************************************************
        //          Effort diagrams - N
        //**************************************************************************************************************
        nLines = numberSnippet ;
        nNodes = nLines + 1;
        deltaX = (double) length/nLines;

        svCout << "Effort diagram - N" << endl;
        svCout << "Lenght" << endl;
        svCout << setw(5) << length << endl;
        svCout << "Snippet - Number" << endl;
        svCout << setw(5) << numberSnippet << endl;
        svCout << "Nodes in the division - Number" << endl;
        svCout << setw(5) << nNodes << endl;
        svCout << "Zoom factor - N" << endl;
        svCout << setw(5) << zoomFactorN << endl;
        svCout << "Delta x" << endl;
        svCout << setw(5) << deltaX << endl;

        double **localDiagN; // In local coordinates
        localDiagN = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            localDiagN[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                localDiagN[i][j] = 0.0;
            }
        }
        double **globalDiagN; // In global coordinates
        globalDiagN = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalDiagN[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalDiagN[i][j] = 0.0;
            }
        }
        double **globalAxisDiagN; // In global coordinates
        globalAxisDiagN = new double *[nNodes];
        for (int i = 0; i < nNodes; i++ ){
            globalAxisDiagN[i] = new double[2];
        }
        for (int i = 0; i < nNodes; i++ ){
            for (int j = 0; j < 2; j++ ){
                globalAxisDiagN[i][j] = 0.0;
            }
        }

        // Initialization of the diagrams calculation
        xDef = 0.0;
        for (int m = 0; m < nNodes; m++){
            // Local
            localDiagN[m][0] = xDef;
            localDiagN[m][1] = zoomFactorN*equationNEd(vectorEdStrengthMaterial[0], xDef);
            xDef += deltaX;

            // Global
            globalDiagN[m][0] = x0 + transf2DMatrix[0][0] * localDiagN[m][0] + transf2DMatrix[1][0] * localDiagN[m][1];
            globalDiagN[m][1] = y0 + transf2DMatrix[0][1] * localDiagN[m][0] + transf2DMatrix[1][1] * localDiagN[m][1];

            // Global Axis
            globalAxisDiagN[m][0] = x0 + transf2DMatrix[0][0] * localDiagN[m][0];
            globalAxisDiagN[m][1] = y0 + transf2DMatrix[0][1] * localDiagN[m][0];
        }

        svCout << "Local - Diagram N" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << localDiagN[i][j];
            }
            svCout << endl;
        }

        svCout << "Global - Diagram N" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalDiagN[i][j];
            }
            svCout << endl;
        }

        svCout << "Global Axis - Diagram N" << endl;
        for (int i = 0; i < nNodes; ++i) {
            svCout << setw(5) << i+1;
            for (int j = 0; j < 2; ++j) {
                svCout << setw(16) << globalAxisDiagN[i][j];
            }
            svCout << endl;
        }

        // For use in MetaPost in LaTex
        string listMPNor = "listMPNor[" + to_string(i) + "] = ";
        for (int i = 0; i < nNodes; ++i) {
            if (i != nNodes-1){
                listMPNor += "(" + to_string(globalDiagN[i][0]) + "*scale+ offsetX, " + to_string(globalDiagN[i][1]) + "*scale+ offsetY)--";
            }else{
                listMPNor += "(" + to_string(globalDiagN[i][0]) + "*scale+ offsetX, " + to_string(globalDiagN[i][1]) + "*scale+ offsetY);";
            }
        }
        svLatexNorCout << listMPNor << endl;
        svLatexNorCout << "draw listMPNor[" + to_string(i) + "] withcolor blue;" << endl;


        // For use in MetaPost in LaTex
        for (int i = 0; i < nNodes; ++i) {
            svLatexNorCout << "mpNor[" + to_string(kNor) + "] = (" + to_string(globalAxisDiagN[i][0]) + "*scale+ offsetX, " + to_string(globalAxisDiagN[i][1]) + "*scale+ offsetY)--(" + to_string(globalDiagN[i][0]) + "*scale+ offsetX, " + to_string(globalDiagN[i][1]) + "*scale+ offsetY);" << endl;
            svLatexNorCout << "draw mpNor["+ to_string(kNor)+"] withcolor blue;" << endl;
            ++kNor;
        }
    }


    //Python
    for (int i = 0; i < numberElements; ++i) {
        svPythonDefCout << "pltpy.plot(defElem" + to_string(i + 1) + "x, defElem" + to_string(i + 1) + "y,'y')" << endl;
    }

    //Latex
    for (int i = 0; i < numberElements; ++i) {
        svLatexDefCout << "draw listMPDef[" + to_string(i) + "] withcolor blue;" << endl;
    }



    svCout.close();
    svPythonDefCout.close();
    svLatexDefCout.close();
    svLatexMomCout.close();
    svLatexSheCout.close();
    svLatexNorCout.close();

    return 0;
}