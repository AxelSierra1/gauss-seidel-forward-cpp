#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

int significant_digits = 5; // precision control (significant digits)

// Function Prototypes (Translated)
double** Read_Matrix(const char* FileName, int* rows, int* cols);
double* Read_Vector(const char* FileName, int* size);
double* Subtract_Vector(double* vec1, double* vec2, int* size);
double Vector_Norm(double* vec, int* size); // Euclidean (L2) norm

void Print_Vector(double* vec, int* size);

double* gauss_seidel_forward(double** A, double* X, double* b, int n, double tolerance){
    double* X_prev = new double[n]; // Stores the previous iteration's vector
    int iteration = 0;
    double current_error; // Norm of the difference between iterations

    do{
        iteration++;

        // Copy current X to X_prev before calculating the new X
        for(int i = 0; i < n; i++){
            X_prev[i] = X[i];
            // Note: Original code unnecessarily reset X[i] to 0 here. 
            // The calculation below correctly overwrites it. Keeping original logic:
            X[i] = 0; 
        }

        // Gauss-Seidel Iteration Formula: X[i] = (1/A[i][i]) * (b[i] - sum(A[i][j]*X[j] for j<i) - sum(A[i][j]*X_prev[j] for j>i))
        // This implementation updates X[i] immediately and uses the updated values for subsequent calculations in the same iteration.
        for(int i = 0; i < n; i++){
            double sigma = 0.0; 
            // Use updated values for j < i
            for(int j = 0; j < i; j++){
                 sigma += A[i][j] * X[j];
            }
            // Use previous iteration's values for j > i
             for(int j = i + 1; j < n; j++){
                 sigma += A[i][j] * X_prev[j];
            }
            // Calculate the new X[i]
            X[i] = (b[i] - sigma) / A[i][i];
        }

        // Calculate the norm of the difference between the current and previous iteration
        double* difference_vector = Subtract_Vector(X_prev, X, &n);
        current_error = Vector_Norm(difference_vector, &n);
        delete[] difference_vector; // Clean up the temporary difference vector

        // PRINT iterations details ==============================================================
        cout << "Previous X (X_prev): " << endl;
        Print_Vector(X_prev, &n);
        cout << endl;

        cout << "Current X:" << endl;
        Print_Vector(X, &n);
        cout << endl;

        cout << fixed << setprecision(significant_digits + 1) << "|| X_prev - X ||: " << current_error << endl << endl;
        // =======================================================================================

        if(iteration == 10000){ // Maximum allowed iterations limit
            cout << "Maximum iterations reached (10000), convergence not achieved within limit.\n";
            break;
        }
    }
    while(current_error > tolerance); // Continue if the error is greater than the tolerance

    if (current_error <= tolerance) {
        cout << "Tolerance reached after (" << iteration << ") iterations:\n\n";
    }

    delete[] X_prev; // Clean up the storage for the previous vector
    return X; // Return the calculated solution vector
}



int main(){
    // Calculate tolerance based on the desired number of significant digits
    double tolerance = (double)1 / (pow(10, significant_digits));  

    int *n_ptr = new int; // Pointer to store the dimension read from files
    
    // Read matrix A from file
    // Note: Assumes matrix is square, using n_ptr for both rows and columns
    double** A = Read_Matrix("data/matrix.txt", n_ptr, n_ptr); 
    
    // Read vector b from file
    double* b = Read_Vector("data/vector.txt", n_ptr);

    // Read initial guess vector x0 from file
    double* x0 = Read_Vector("data/initial_vec.txt", n_ptr);

    // Check if reads were successful (basic check: n > 0)
    if (*n_ptr <= 0) {
        cerr << "Error reading matrix/vector dimensions from file." << endl;
        // Perform necessary cleanup before exiting
        delete n_ptr;
        // Potentially A, b, x0 might be null or partially allocated, handle carefully
        // Basic cleanup assuming they might have been allocated if n was read > 0 initially:
        if (A) { /* Add deletion logic if needed */ }
        if (b) delete[] b;
        if (x0) delete[] x0;
        return 1; // Indicate error
    }
    
    int n_val = *n_ptr; // Get the actual dimension value

    // Solve the system using Gauss-Seidel
    double* solution = gauss_seidel_forward(A, x0, b, n_val, tolerance);

    cout << "Approximate solution to the system: \n";
    Print_Vector(solution, &n_val);
    cout << endl;


    // Memory Cleanup
    // Delete matrix A
    for(int i = 0; i < n_val; i++){
        delete[] A[i];
    }
    delete[] A;

    // Delete dimension pointer
    delete n_ptr; 

    // Delete vectors
    delete[] b;
    // x0 is now pointing to the solution returned by gauss_seidel_forward,
    // so deleting solution means deleting the memory originally allocated for x0.
    // No need to delete x0 separately IF gauss_seidel_forward returns the modified input x0.
    // If gauss_seidel_forward allocates new memory for the result, you'd delete both x0 and solution.
    // Based on the implementation, it modifies and returns the input X (which was x0).
    // So, deleting solution (or x0) is correct here. Let's stick to deleting `solution`.
    delete[] solution; 

    return 0; // Indicate successful execution
}

double** Read_Matrix(const char* FileName, int* rows, int* cols){
    ifstream fin(FileName);
    if (!fin.is_open()) {
        cerr << "Error: Could not open matrix file: " << FileName << endl;
        *rows = 0;
        *cols = 0;
        return nullptr; // Indicate failure
    }

    fin >> *rows;
    fin >> *cols;

    if (*rows <= 0 || *cols <= 0) {
         cerr << "Error: Invalid matrix dimensions read from file: " << *rows << "x" << *cols << endl;
         fin.close();
         return nullptr; // Indicate failure
    }

    double** A = new double* [*rows];

    for(int i = 0; i < *rows; i++){
        A[i] = new double[*cols];
        for(int j = 0; j < *cols; j++){
            if (!(fin >> A[i][j])) { // Check if read was successful
                 cerr << "Error: Failed to read value for A[" << i << "][" << j << "] from file " << FileName << endl;
                 // Cleanup partially allocated memory
                 for(int k = 0; k < i; ++k) {
                     delete[] A[k];
                 }
                 delete[] A;
                 *rows = 0;
                 *cols = 0;
                 fin.close();
                 return nullptr; // Indicate failure
            }
        }
    }
    fin.close();
    return A;
};

double* Read_Vector(const char* FileName, int* size){
    ifstream fin(FileName);
     if (!fin.is_open()) {
        cerr << "Error: Could not open vector file: " << FileName << endl;
        *size = 0;
        return nullptr; // Indicate failure
    }

    fin >> *size;

     if (*size <= 0) {
         cerr << "Error: Invalid vector size read from file: " << *size << endl;
         fin.close();
         return nullptr; // Indicate failure
    }

    double* v = new double [*size];

    for(int i = 0; i < *size; i++){
         if (!(fin >> v[i])) { // Check if read was successful
             cerr << "Error: Failed to read value for v[" << i << "] from file " << FileName << endl;
             // Cleanup allocated memory
             delete[] v;
             *size = 0;
             fin.close();
             return nullptr; // Indicate failure
         }
    }
    fin.close();

    return v;
};

double* Subtract_Vector(double* vec1, double* vec2, int* size){
    double* result = new double [*size];
    for(int i = 0; i < *size; i++){
        result[i] = vec1[i] - vec2[i]; 
    }
    return result;
}

double Vector_Norm(double* vec, int* size){
    double norm_sq = 0.0; // Use double for intermediate sum
    for(int i = 0; i < *size; i++){
        norm_sq += (vec[i] * vec[i]);
    }
    return sqrt(norm_sq);
}

void Print_Vector(double* vec, int* size){
    cout << "[";
    for(int i = 0; i < *size; i++){
        // Set precision for each element
        cout << fixed << setprecision(significant_digits) << vec[i];
        if(i < *size - 1){
            cout << ",\n "; // Print newline and space for alignment (optional)
        }
    }
     cout << "]" << endl; // Close bracket and add final newline
};