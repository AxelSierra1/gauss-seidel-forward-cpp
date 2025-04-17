# Gauss-Seidel Linear System Solver
This C++ program implements the **Gauss-Seidel iterative method** to find an approximate solution for a system of linear equations of the form **Ax = b**.
The program reads the coefficient matrix (A), the constant vector (b), and an initial guess vector (x0) from separate text files. It then iteratively refines the solution until a specified tolerance is met or a maximum number of iterations is reached.

## Features
*   Solves square systems of linear equations (Ax = b) using the Gauss-Seidel method.
*   Reads input matrix A, vector b, and initial guess x0 from specified text files (`data/matriz.txt`, `data/vector.txt`, `data/vec_inicial.txt`).
*   Allows configuration of numerical precision (significant digits) for calculations and output.
*   Calculates tolerance based on the specified number of significant digits.
*   Outputs iteration details, including the previous vector (X_prev), current vector (X), and the norm of their difference (`|| X_prev - X ||`).
*   Includes a safeguard against infinite loops by limiting the maximum number of iterations.
*   Basic file I/O error checking.

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
