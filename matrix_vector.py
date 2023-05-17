from math import sin

class Matrix:
    def __init__(self, N):
        self.N = N

    def generate_matrix(self, a1, a2, a3):
        A = []
        for i in range(self.N):
            r = []
            for j in range(self.N):
                if i-2 == j or i+2 == j:
                    r.append(a3)
                elif i-1 == j or i+1 == j:
                    r.append(a2)
                elif i == j:
                    r.append(a1)
                else:
                    r.append(int(0))
            A.append(r)
        return A

    def generate_vector(self):
        vectorB = []
        for i in range(self.N):
            vectorB.append(sin(i*9))
        return vectorB

    def temp_matrix(self, A):
        M = []
        for i in A:
            r = []
            for e in i:
                r.append(e)
            M.append(r)
        return M

    def temp_vector(self, b):
        v = []
        for i in range(len(b)):
            v.append(b[i])
        return v

    def vec_zero(self, M):
        vec = M*[0]
        return vec

    def vector_one(self, vec):
        return vec*[1]

    def matrix_zero(self, a, b):
        matrix = []
        for i in range(a):
            r = a*[0]
            matrix.append(r)
        return matrix

    def dot_product(self, A, b):
        tempA = A
        tempb = b
        vec = self.vec_zero(len(A))

        for i in range(len(A)):
            for j in range(len(A[0])):
                vec[i] += tempA[i][j]*tempb[j]
        return vec

    def subtract(self, b, dot):
        vec = self.temp_vector(dot)
        for i in range(len(b)):
            vec[i] -= b[i]
        return vec

    def matrixL(self, x):
        matrix = self.matrix_zero(len(x), len(x))
        for i in range (len(x)):
            matrix[i][i] = x[i]

        return matrix