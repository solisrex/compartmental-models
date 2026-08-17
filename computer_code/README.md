# Running the Python scripts

To run these Python scripts, you'll need to have the following packages installed:

1. `sympy` (this package handles all the symbolic operations),
2. `lark` (this one parses user inputs),
3. `nicegui` (this one handles the user interface).

# What the scripts do and why they are useful

The `interactive_solver.py` script allows a user to assign values to variables. It can also update the system of equations after making such assignments. This is useful when searching for the solutions to a system of equations.

The `matrix_manipulator.py` script allows a user to interact with matrices. Some of the classes of Darboux polynomials we found were expressed in terms of the kernel of a matrix. This script helps find such solutions.

The `integral_finder.py` script allows a user to select and combine different Darboux polynomials together. Our method for finding conserved quantities is based on the existence of Darboux polynomials with linearly dependent cofactors. The integral finder script helps with this process.
