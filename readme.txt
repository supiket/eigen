This program takes an input matrix and finds its dominant eigenvalue, eigenvector and its second dominant eigenvalue.
Three command line parameters are expected, respectively:
1) The name of the file that the matrix will be read from,
2) The value of tolerance for the eigenvalue,
3) The name of the file that the output will be written to.
Normalized power iteration algorithm and Hotelling's deflation are implemented during the computations.
Normalized power iteration algorithm finds a matrix's dominant eigenvector and eigenvalue.
Hotteling's deflation is a method of deflation that finds a matrix's desired eigenvector and eigenvalue.
It is used to find second dominant eigenvalue in this project.
The output has the form (don't mind these brackets: [])
[OUTPUT]
Eigenvalue #1: [dominant eigenvalue]
[dominant eigenvector's first element
second element
third element
.
.
.]
Eigenvalue #2: [second dominant eigenvalue]
[/OUTPUT]
The input matrix is assumed to be a square matrix, and should be so.
The output is written to a txt file with the name given as a parameter.
