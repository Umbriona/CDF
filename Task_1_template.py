# MTF072 Computational Fluid Dynamics
# Task 1: diffusion equation
# Template prepared by:
# Gonzalo Montero Villar
# Department of Mechanics and Maritime Sciences
# Division of Fluid Dynamics
# villar@chalmers.se
# November 2019
# Packages needed
import numpy as np
import matplotlib.pyplot as plt
def main():
    #================= Inputs =====================
    # Geometric inputs
    nI =  # number of nodes X direction.
    nJ =  # number of nodes Y direction.
    grid_type = 'equidistant' # this sets equidistant mesh sizing or non-
equidistant
    xL =  # length of the domain in X direction
    yL =  # length of the domain in Y direction
    # Solver inputs
    nIterations  =  # maximum number of iterations
    resTolerance =  # convergence criteria for residuals each variable
    # ================ Code =======================
    # For all the matrices the first input makes reference to the x coordinate
    # and the second input to the y coordinate, [i+1] is east and [j+1] north
    # Allocate all needed variables
    cI = nI + 1                    # number of cells in the X direction. Cells 
                                   # added in the boundaries
    cJ = nJ + 1                    # number of cells in the Y direction. Cells 
                                   # added in the boundaries
    coeffsT = np.zeros((cI,cJ,5))  # coefficients for temperature
                                   # E, W, N, S and P
    S_U     = np.zeros((cI,cJ))    # source term for temperature
    S_P     = np.zeros((cI,cJ))    # source term for temperature
    T       = np.zeros((cI,cJ))    # temperature matrix
    k       = np.zeros((cI,cJ))    # coefficient of conductivity
    q       = np.zeros((cI,cJ,2))  # heat flux, first x and then y component
    residuals = [] # List containing the value of the residual for each iteration
    # Generate mesh and compute geometric variables
    # Allocate all variables matrices
    xCoords_C = np.zeros((cI,cJ)) # X coords of the cells
    yCoords_C = np.zeros((cI,cJ)) # Y coords of the cells
    xCoords_N = np.zeros((nI,nJ)) # X coords of the nodes
    yCoords_N = np.zeros((nI,nJ)) # Y coords of the nodes
    dxe_C     = np.zeros((cI,cJ)) # X distance to east cell
    dxw_C     = np.zeros((cI,cJ)) # X distance to west cell
    dyn_C     = np.zeros((cI,cJ)) # Y distance to north cell
    dys_C     = np.zeros((cI,cJ)) # Y distance to south cell
    dx_C      = np.zeros((cI,cJ)) # X size of the cell
    dy_C      = np.zeros((cI,cJ)) # Y size of the cell
    if grid_type == 'equidistant':
        # Cell size
        dx = xL/(nI - 1)
        dy = yL/(nJ - 1)
        # Fill the coordinates
        for i in range(nI):
            for j in range(nJ):
                # For the nodes
                xCoords_N[i,j] = i*dx
                yCoords_N[i,j] = j*dy
                # For the cells
                if i > 0:
                    xCoords_C[i,j] = 0.5*(xCoords_N[i,j] + xCoords_N[i-1,j])
                if i == (nI-1) and j>0:
                    yCoords_C[i+1,j] = 0.5*(yCoords_N[i,j] + yCoords_N[i,j-1])
                if j >0:
                    yCoords_C[i,j] = 0.5*(yCoords_N[i,j] + yCoords_N[i,j-1])
                if j == (nJ-1) and i>0:
                    xCoords_C[i,j+1] = 0.5*(xCoords_N[i,j] + xCoords_N[i-1,j])
                # Fill dx_C and dy_C
                if i>0:
                    dx_C[i,j] = xCoords_N[i,j] - xCoords_N[i-1,j]
                if j>0:
                    dy_C[i,j] = yCoords_N[i,j] - yCoords_N[i,j-1]
    elif grid_type == 'non-equidistant':
        rx = 1.15
        ry = 1.15
        # Fill the necessary code to generate a non equidistant grid and
        # fill the needed matrixes for the geometrical quantities
    xCoords_C[-1,:] = xL
    yCoords_C[:,-1] = yL
    # Fill dxe, dxw, dyn and dys
    for i in range(1,cI - 1):
        for j in range(1,cJ - 1):
            dxe_C[i,j] = 
            dxw_C[i,j] = 
            dyn_C[i,j] = 
            dys_C[i,j] = 
    # Initialize variable matrices and boundary conditions
    # Looping
    for iter in range(nIterations):
        # Update conductivity coefficient matrix, k
        # Update source term matrix according to your case
        # Compute coefficients (taking into account boundary conditions)
        for i in range(1,cI-1):
            for j in range(1,cJ-1):    
                coeffsT[i,j,0] = 
                coeffsT[i,j,1] = 
                coeffsT[i,j,2] = 
                coeffsT[i,j,3] = 
                coeffsT[:,:,4] =        
        # Solve for T Gauss-Seidel
        # Copy T to boundaries where homegeneous Neumann needs to be applied
        # Compute residuals (taking into account normalization)
        r = 0
        residuals.append(r)
        print('iteration: %d\nresT = %.5e\n\n'  % (iter, residuals[-1]))
        #  Check convergence
        if resTolerance>residuals[-1]:
            break
    # Compute heat fluxes
    for i in range(1,cI-1):
        for j in range(1,cJ-1):
            q[i,j,0] = 
            q[i,j,1] = 
    # Plotting section (these are some examples, more plots might be needed)
    # Plot results
    plt.figure()
    # Plot mesh
    plt.subplot(2,2,1)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Computational mesh')
    plt.axis('equal')
    # Plot temperature contour
    plt.subplot(2,2,2)
    plt.title('Temperature [ºC]')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    # Plot residual convergence
    plt.subplot(2,2,3)
    plt.title('Residual convergence')
    plt.xlabel('iterations')
    plt.ylabel('residuals [-]')
    plt.title('Residual')
    # Plot heat fluxes
    plt.subplot(2,2,4)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Heat flux')
    plt.axis('equal')
    plt.show()
if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
