# Finite-Difference-Temperature-Distribution-of-a-Plate
A temperature distribution of a square plate with different temperatures on each side is calculated with Finite Difference method at the equilibrium state.

The length of the square, the grid length and the temperatures of the sides are determined at the beginning. By default, these values would be 1 meter for the leght with 0.1 meter grids. The temperatures will be 50 for the right, 25 for the left, 100 for the up and 0 for the down side. These values can be changed manually. 

By using the parameters, the number of nodes are calculated and two matrices are defined for the solution. Step by step, the corners, the sides and the interior is calculated with finite difference method. The temperature difference problem is solved with partial pivoting and backpropogation. The table is formed and plotted at the end. 
