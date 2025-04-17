# Gauss-Seidel Linear System Solver
This C++ program implements the **Gauss-Seidel iterative method** to find an approximate solution for a system of linear equations of the form **Ax = b**.
The program reads the coefficient matrix (A), the constant vector (b), and an initial guess vector (x0) from separate text files. It then iteratively refines the solution until a specified tolerance is met or a maximum number of iterations is reached.

## Features
*   Solves square systems of linear equations (Ax = b) using the Gauss-Seidel method.
*   Reads input matrix A, vector b, and initial guess x0 from specified text files (`data/matrix.txt`, `data/vector.txt`, `data/initial_vec.txt`).
*   Allows configuration of numerical precision (significant digits) for calculations and output.
*   Calculates tolerance based on the specified number of significant digits.
*   Outputs iteration details, including the previous vector (X_prev), current vector (X), and the norm of their difference (`|| X_prev - X ||`).
*   Includes a safeguard against infinite loops by limiting the maximum number of iterations.
*   Basic file I/O error checking.

## Input File Format

The program expects input files to be located in a subdirectory named `data/` relative to the executable.

1.  **`data/matrix.txt` (Coefficient Matrix A)**
    *   Line 1: Integer `N` (number of rows)
    *   Line 2: Integer `N` (number of columns - must be equal to rows for a square system)
    *   Lines 3 to N+2: Each line contains `N` space-separated floating-point numbers representing a row of the matrix A.

    *Example (`data/matrix.txt` for a 3x3 system):*
    ```
    3
    3
    10.0 -1.0 2.0
    -1.0 11.0 -1.0
    2.0 -1.0 10.0
    ```

2.  **`data/vector.txt` (Constant Vector b)**
    *   Line 1: Integer `N` (size of the vector)
    *   Lines 2 to N+1: Each line contains one floating-point number representing an element of the vector b.

    *Example (`data/vector.txt` for N=3):*
    ```
    3
    6.0
    25.0
    -11.0
    ```

3.  **`data/initial_vec.txt` (Initial Guess Vector x0)**
    *   Line 1: Integer `N` (size of the vector)
    *   Lines 2 to N+1: Each line contains one floating-point number representing an element of the initial guess vector x0. Often, this is a zero vector.

    *Example (`data/initial_vec.txt` for N=3):*
    ```
    3
    0.0
    0.0
    0.0
    ```

**Important:** The dimension `N` must be consistent across all three files.

*   ## Structure
```
gauss-seidel-forward-cpp/
├── data/
│   ├── initial_vec.txt
│   ├── matrix.txt
│   └── vector.txt
├── main.cpp
└── README.md
```

## How to compile and run
Create the build directory if it doesn't exist
```
mkdir -p build
```
Compile
```
g++ -o build/gauss_seidel main.cpp -Iinclude
```
Run
```
./build/gauss_seidel
```
## Configuration
- Precision: You can change the number of significant digits used for calculations and output by modifying the significant_digits global variable at the top of the main.cpp file. The tolerance for convergence is automatically calculated based on this value (tolerance = 1 / 10^significant_digits).
- Maximum Iterations: The maximum number of iterations is currently hardcoded to 10000 within the gauss_seidel_forward function.

## Algorithm Notes
The Gauss-Seidel method is an iterative technique. It works best for matrices that are strictly diagonally dominant or symmetric positive definite, as convergence is guaranteed in these cases. For other matrices, the method might converge slowly, diverge, or oscillate.
