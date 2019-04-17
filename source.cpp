#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class TheMatrix { //defining matrix class
public:
    int numrow, numcol; //row and column numbers of the matrix class
    void allocate(); //function to allocate memory dynamically
    double** currmatx; //double to allocate memory dynamically
    TheMatrix(int numrow, int numcol); //constructor with number of rows and columns parameters   
    TheMatrix(int numrow); //another constructor with one parameter (to be able to create identity matrix, we are going to get there)
    TheMatrix transpose(); //transpose function
};

void TheMatrix::allocate() { //allocates space for TheMatrix class
    currmatx = new double*[numrow]; 
    for (int i = 0; i < numrow; i++)
        currmatx[i] = new double[numcol];
};

TheMatrix::TheMatrix(int row, int col) : numrow(row), numcol(col) { //the constructor of matrix class with row and column parameters
    allocate(); //only allocates space with given row and column numbers
}

TheMatrix::TheMatrix(int row) : numrow(row) { //the constructor of matrix class with one parameter
    allocate(); //allocates space with given row and column numbers
    //makes every element equal to zero. going to be used to create identity matrix later
    for (int i = 0; i < numrow; i++) { 
        for (int j = 0; j < numrow; j++) {
            currmatx[i][j] = 0; 
        }
    }
}

TheMatrix transpose(TheMatrix mat) { //takes the transpose of a given matrix
    TheMatrix transposed(mat.numcol, mat.numrow);
    for (int i = 0; i < mat.numrow; ++i) {
        for (int j = 0; j < mat.numcol; ++j) {
            transposed.currmatx[j][i] = mat.currmatx[i][j];
        }
    }
    return transposed; //returns the transposed matrix
}

double norminf(TheMatrix mat) { //finds a matrix's norm at infinity
	//the norm of a matrix in infinity is equal to maximum of its absolute row sums
    double norminf = 0;
    double rowsum = 0;
    for (int i = 0; i < mat.numrow; ++i) {
        rowsum = 0;
        for (int j = 0; j < mat.numcol; ++j) {
            rowsum = rowsum + fabs(mat.currmatx[i][j]); 
        }
        if (rowsum > norminf) { //makes the norminf equal to the biggest rowsum at the end of the iteration
            norminf = rowsum;
        }
    }
    return norminf; //returns the norm at infinity
}

TheMatrix operator*(const TheMatrix& mat1, const TheMatrix& mat2) { //the multiplication operator for two matrices
    TheMatrix temp(mat1.numrow); //creates the product matrix which is temp, and has elements all zero
    //takes the product of two matrices with dimensions nxm and mxk, giving a result matrix of nxk
    for (int i = 0; i < temp.numrow; ++i) {
        for (int j = 0; j < temp.numcol; ++j) {
            for (int k = 0; k < mat1.numcol; ++k) {
                temp.currmatx[i][j] = temp.currmatx[i][j] + (mat1.currmatx[i][k] * mat2.currmatx[k][j]);
            }
        }
    }
    return temp; //returns the product matrix
}

TheMatrix operator*(const TheMatrix& mat, double factor) { //the multiplication operator for a matrix and a scalar
    for (int i = 0; i < mat.numrow; ++i) {
        for (int j = 0; j < mat.numcol; ++j) {
            mat.currmatx[i][j] = mat.currmatx[i][j] * factor; //takes the product of each element with the factor
        }
    }
    return mat; //returns the product matrix
}

TheMatrix operator/(const TheMatrix& mat, double divisor) { //the division operator of a matrix by a scalar
    for (int i = 0; i < mat.numrow; ++i) {
        for (int j = 0; j < mat.numcol; ++j) {
            mat.currmatx[i][j] = mat.currmatx[i][j] / divisor; //divides each element of the matrix by the divisor
        }
    }
    return mat; //returns the result matrix
}

TheMatrix operator+(const TheMatrix& mat1, const TheMatrix& mat2) { //the addition operator of two matrices
    for (int i = 0; i < mat1.numrow; ++i) {
        for (int j = 0; j < mat1.numcol; ++j) {
            mat1.currmatx[i][j] = mat1.currmatx[i][j] + mat2.currmatx[i][j]; //adds each element of the second matrix to the first one's
        }
    }
    return mat1; //returns the result matrix
}

TheMatrix operator-(const TheMatrix& mat1, const TheMatrix& mat2) { //the substraction operator of two matrices
    for (int i = 0; i < mat1.numrow; ++i) {
        for (int j = 0; j < mat1.numcol; ++j) {
            mat1.currmatx[i][j] = mat1.currmatx[i][j] - mat2.currmatx[i][j]; //substracts each element of the second matrix from the first one's
        }
    }
    return mat1; //returns the result matrix
}

int main(int argc, char** argv) { //main function
    //these three variables are uninitialized because they will be determined by the command line arguments
    char* namein; //keeps the name of the txt file of input matrix
    char* nameout; //keeps the name of the txt file of output
    double tolerance; //keeps the tolerance value
    if (argc != 4) { //this program requires 3 inputs, which corresponds to argc = 4. it quits otherwise
        cout << "Exactly 3 inputs are required.";
        return 0;
    }
    namein = argv[1]; //the name of the input file
    tolerance = atof(argv[2]); //the value of tolerance
    nameout = argv[3]; //the name of the output file
    ifstream readfile(namein); //will read input file
    if (readfile.is_open() == 0) { //prints out an error message and quits if there is a problem with opening the file
        cout << "Input file could not be opened.";
        return 0;
    }
    string line; //will be used to find the dimensions of the input matrix
    int n = 0;
    while (getline(readfile, line)) {
        n++; //n becomes the dimensions of the input matrix, since its row and column numbers are equal
    }
    readfile.clear(); 
    readfile.seekg(0, ios::beg); //after finding the dimension, sets the cursor to the beginning.
    TheMatrix matA(n, n); //creates the input matrix
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) {
            readfile >> matA.currmatx[i][j];
        }
    }
    TheMatrix idenmat(n); //n by n identity matrix, is going to be used during the power iteration method
    for (int i = 0; i < n; i++) {
        idenmat.currmatx[i][i] = 1; //sets the elements at the diagonal equal to 1
    }
    readfile.close(); //closes the read file 
    TheMatrix eigenvector(n, 1); //creates a vector that is going to be the eigenvector
    for (int i = 0; i < eigenvector.numrow; ++i) {
        for (int j = 0; j < eigenvector.numcol; ++j) {
            eigenvector.currmatx[i][j] = 1; //sets the elements of the eigenvector equal to 1, in order not to have a problem during division
        }
    }
    double eigenvalue = 0; //creates a double that is going to be the biggest eigenvalue
    double tempvalue = 0; //is going to be used in the do while below

    //power iteration algorithm
    do {
        eigenvector = matA * eigenvector; //multiplies the first placeholder vector with the input matrix
        tempvalue = eigenvalue; //sets its value to another variable
        eigenvalue = 0; //sets the eigenvalue equal to zero at each iteration, otherwise errors occur
        eigenvalue = norminf(eigenvector); //eigenvalue becomes the norm of eigenvector
        eigenvector = eigenvector / eigenvalue; //eigenvector becomes a unit vector
    } while (fabs(eigenvalue - tempvalue) > tolerance); //while the difference between the last value of eigenvalue and the previous one is greater than the tolerance
    TheMatrix matforsign = matA * eigenvector; //matforsign matrix is going to be used to determine the sign of eigenvalue
    //if the signs of the first element of the eigenvector and the first element of the matrix which is the multiplication above are different
    //inverts the sign of the eigenvalue
    if ((matforsign.currmatx[0][0] * eigenvector.currmatx[0][0]) < 0) 
        eigenvalue = eigenvalue * (-1);
    ofstream writefile(nameout); //now that the first eigenvalue and eigenvector are found, writes them into to output file
    if (writefile.is_open() == 0) { //prints out an error message and quits if a problem is detected during opening the file
        cout << "Output file could not be opened.";
        return 0;
    }
    writefile<<"Eigenvalue#1:"<<eigenvalue<<"\n"; //writes eigenvalue into the output file
    for(int i=0;i<n;i++){//writes the eigenvector into the output file
        writefile<<eigenvector.currmatx[i][0]<<"\n";       
    }

    //hotelling's deflation algorithm   
    double alpha = 0;
    for (int i = 0; i < eigenvector.numrow; ++i) { //the variable alpha becomes the addition of square of each element of the eigenvector
        alpha = alpha + eigenvector.currmatx[i][0] * eigenvector.currmatx[i][0];
    } 
    alpha = sqrt(alpha); //then its square is taken, and now is the norm of eigenvector
    TheMatrix uniteigvec(n, 1); //creates a unit vector with the same dimensions as eigenvector
    uniteigvec = eigenvector; 
    uniteigvec = uniteigvec / alpha; //now uniteigvec is the unit vector of eigenvector
    TheMatrix deflation(n, n); //creates a deflation matrix
    deflation = uniteigvec * transpose(uniteigvec); //deflation matrix becomes the product of the unit eigenvector and its transpose
    deflation = deflation * eigenvalue; //then it is multiplied with eigenvalue
    deflation = matA - deflation; //gives the final deflation matrix
    double eigenvalue2, tempvalue2; //the second eigenvalue
    TheMatrix eigenvector2(n, 1); //the second eigenvector
    for (int i = 0; i < eigenvector2.numrow; ++i) { //same a before, setsthe elements of the eigenvector equal to 1, in order not to have a problem during division
        for (int j = 0; j < eigenvector2.numcol; ++j) {
            eigenvector2.currmatx[i][j] = 1;
        }
    }
    do { //same process as the finding the first eigenvector
        eigenvector2 = deflation * eigenvector2;
        tempvalue2 = eigenvalue2;
        eigenvalue2 = norminf(eigenvector2);
        eigenvector2 = eigenvector2 / eigenvalue2;
    } while (fabs(eigenvalue2 - tempvalue2) > tolerance); //now the second eigenvalue and eigenvectors are found
    writefile<<"Eigenvalue#2:"<<eigenvalue2<<"\n"; //prints out the second eigenvalue
    return 0; //terminates
}
