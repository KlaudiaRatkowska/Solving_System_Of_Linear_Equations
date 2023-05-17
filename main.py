import methods
from matrix_vector import *
import matplotlib.pyplot as plt

if __name__ == '__main__':

    size = 939
    matrix = Matrix(size)
    findCourage = True

    #taskA
    A = matrix.generate_matrix(11, -1, -1)
    b = matrix.generate_vector()

    #taskB
    res_jacobi = methods.jacobi(A, b, findCourage)
    res_gauss_seidl = methods.gauss_seidl(A, b, findCourage)
    plt.plot(res_jacobi, color="pink", label="Jacobi residuum")
    plt.plot(res_gauss_seidl, color="green", label="Gauss-Seidl residuum")
    plt.xlabel("iterations")
    plt.ylabel("Residuum norm")
    plt.legend()
    plt.yscale('log')
    plt.title("Comparison of residuum in Jacobi and Gauss-Seidl methods")
    plt.savefig("taskB.png")
    plt.show()
    plt.close()


    #taskC
    C = matrix.generate_matrix(3, -1, -1)
    findCourage = True
    jacobi_courage = methods.jacobi(C, b, findCourage)
    gauss_courage = methods.gauss_seidl(C, b, findCourage)
    plt.plot(jacobi_courage, color="blue", label="Jacobi")
    plt.plot(gauss_courage, color="purple", label="Gauss-Seidl")
    plt.xlabel("Iterations")
    plt.ylabel("Residuum norm")
    plt.legend()
    plt.yscale('log')
    plt.title("Convergence of methods in task C")
    plt.savefig("taskC.png")
    plt.show()
    plt.close()

    #taskD
    methods.lu_decomposition(C, b)

    #taskE
    N = [100, 500, 1000, 2000, 3000]
    jacobi = []
    gauss_seidl = []
    lu = []
    i = 0
    new_matrix = []
    findCourage = False
    for size in N:
        print("Matrix size: ", size)
        new_matrix.append(Matrix(size))
        matrixA = new_matrix[i].generate_matrix(11, -1, -1)
        vectorb = new_matrix[i].generate_vector()

        jacobi.append(methods.jacobi(matrixA, vectorb, findCourage))
        gauss_seidl.append(methods.gauss_seidl(matrixA, vectorb, findCourage))
        lu.append(methods.lu_decomposition(matrixA, vectorb))
        i += 1


    plt.plot(N, jacobi, color="blue", label="Jacobi")
    plt.plot(N, gauss_seidl, color="orange", label="Gauss-Seidl")
    plt.plot(N, lu, color="green", label="LU Decomposition")
    plt.legend()
    plt.xlabel("matrix size")
    plt.ylabel("Time [s]")
    plt.title("Methods execution time")
    plt.savefig("methods.png")
    plt.show()
    plt.close()