import time
import math
from matrix_vector import *

residuum = pow(10, -9)


def norm(res):
    sum = 0
    for e in res:
        sum += pow(e, 2)
    return math.sqrt(sum)


def jacobi(A, b, findCourage):
    print("Jacobi method")
    iter = 0
    size = len(A)
    matrix = Matrix(size)
    tempA = matrix.temp_matrix(A)
    tempb = matrix.temp_vector(b)
    vec_x = matrix.vec_zero(len(tempA))
    tempx = matrix.temp_vector(vec_x)
    courage = []

    start = time.time()
    while True:
        for i in range(len(tempA)):
            value = tempb[i]
            for j in range(len(tempA)):
                if i != j:
                    value -= tempA[i][j] * vec_x[j]

            value /= tempA[i][i]
            tempx[i] = value

        vec_x = matrix.temp_vector(tempx)
        dot = matrix.dot_product(tempA, vec_x)
        res = matrix.subtract(tempb, dot)
        res_norm = norm(res)
        courage.append(res_norm)
        if iter >= 50 and findCourage :
            return courage

        if res_norm < residuum:
            break

        iter += 1

    finish = time.time() - start
    print("Number of iterations: ", iter)
    print("Time: ", finish)
    print("Residuum: ", norm(res))
    print(" ")
    if findCourage:
        return courage
    return finish


def gauss_seidl(A, b, findCourage):
    print("Gauss-Seidl method")
    iter = 0
    size = len(A)
    matrix = Matrix(size)
    tempA = matrix.temp_matrix(A)
    tempb = matrix.temp_vector(b)
    vec_x = matrix.vec_zero(len(tempA))
    courage = []

    start = time.time()
    while True:
        for i in range(len(tempA)):
            value = tempb[i]
            for j in range(len(tempA)):
                if i != j:
                    value -= tempA[i][j] * vec_x[j]

            value /= tempA[i][i]
            vec_x[i] = value

        dot = matrix.dot_product(tempA, vec_x)
        res = matrix.subtract(tempb, dot)
        res_norm = norm(res)
        courage.append(res_norm)
        if iter >= 50 and findCourage :
            return courage

        if res_norm < residuum:
            break

        iter += 1

    finish = time.time() - start
    print("Number of iterations: ", iter)
    print("Time: ", finish)
    print("Residuum: ", norm(res))
    print(" ")
    if findCourage:
        return courage
    return finish


def lu_decomposition(C, b):
    print("LU decomposition method")
    size = len(C)
    matrix = Matrix(size)
    tempA = matrix.temp_matrix(C)
    tempb = matrix.temp_vector(b)
    tempx = matrix.vector_one(len(C))
    tempy = matrix.vec_zero(len(C))
    matrixL = matrix.matrixL(matrix.vector_one(len(C)))
    matrixU = matrix.matrix_zero(len(C), len(C))

    start = time.time()
    for i in range(len(C)):
        for j in range(i + 1):
            matrixU[j][i] += tempA[i][j]
            for k in range(j):
                matrixU[j][i] -= matrixL[j][k] * matrixU[k][i]

        for j in range(i + 1, len(C)):
            for k in range(i):
                matrixL[j][i] -= matrixL[j][k] * matrixU[k][i]

            matrixL[j][i] += tempA[j][i]
            matrixL[j][i] /= matrixU[i][i]

    for i in range(len(C)):
        value = tempb[i]
        for j in range(i):
            if i != j:
                value -= matrixL[i][j] * tempy[j]
        tempy[i] = value / matrixL[i][i]

    for i in range(len(C) - 1, -1, -1):
        value = tempy[i]
        for j in range(i + 1, len(C)):
            if i != j:
                value -= matrixU[i][j] * tempx[j]
        tempx[i] = value / matrixU[i][i]

    dot = matrix.dot_product(tempA, tempx)
    res = matrix.subtract(dot, tempb)
    finish = time.time() - start
    print("Time: ", finish)
    print("Residuum: ", norm(res))
    print(" ")
    return finish
