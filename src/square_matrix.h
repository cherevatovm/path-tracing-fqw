#ifndef SQ_MAT_H
#define SQ_MAT_H
#include <math.h> 

template <int N>
class SquareMatrix {
public:
    double contents[N][N];

    SquareMatrix() {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                contents[i][j] = (i == j) ? 1 : 0;
    }

    SquareMatrix(std::initializer_list<std::initializer_list<double>> mat) {
        size_t i = 0, j = 0;
        for (const auto& row : mat) {
            j = 0;
            for (double val : row) {
                contents[i][j] = val;
                ++j;
            }
            ++i;
        }
    }

    SquareMatrix(const double mat[N][N]) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                contents[i][j] = mat[i][j];
    }

    SquareMatrix operator+(const SquareMatrix& m) const {
        SquareMatrix res = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                res.contents[i][j] += m.contents[i][j];
        return res;
    }

    SquareMatrix operator*(const SquareMatrix& m) const {
        SquareMatrix res = zero_matr();
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k)
                    res.contents[i][j] += contents[i][k] * m.contents[k][j];
        return res;
    }

    SquareMatrix operator*(double s) const {
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.contents[i][j] *= s;
        return r;
    }

    SquareMatrix operator/(double s) const {
        SquareMatrix r = *this;
        if (abs(s) > 0.00001) {
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    r.contents[i][j] /= s;
        }
        return r;
    }

    bool operator==(const SquareMatrix<N>& m2) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (contents[i][j] != m2.contents[i][j])
                    return false;
        return true;
    }

    bool operator<(const SquareMatrix<N>& m2) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                if (contents[i][j] < m2.contents[i][j])
                    return true;
                if (contents[i][j] > m2.contents[i][j])
                    return false;
            }
        return false;
    }

    SquareMatrix transpose() {
        SquareMatrix transposed;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                transposed.contents[j][i] = contents[i][j];
            }
        }
        return transposed;
    }

    static SquareMatrix zero_matr() {
        SquareMatrix m;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m.contents[i][j] = 0;
        return m;
    }

    double determinant3x3() {
        //if (N == 3)
        return
            contents[0][0] * (contents[1][1] * contents[2][2] - contents[1][2] * contents[2][1]) -
            contents[0][1] * (contents[1][0] * contents[2][2] - contents[1][2] * contents[2][0]) +
            contents[0][2] * (contents[1][0] * contents[2][1] - contents[1][1] * contents[2][0]);
        //else
    }

    SquareMatrix<3> inverse3x3() {
        double inv_det = 1.0 / determinant3x3();
        SquareMatrix<3> inv_matrix;

        inv_matrix.contents[0][0] = (contents[1][1] * contents[2][2] - contents[1][2] * contents[2][1]) * inv_det;
        inv_matrix.contents[0][1] = (contents[0][2] * contents[2][1] - contents[0][1] * contents[2][2]) * inv_det;
        inv_matrix.contents[0][2] = (contents[0][1] * contents[1][2] - contents[0][2] * contents[1][1]) * inv_det;

        inv_matrix.contents[1][0] = (contents[1][2] * contents[2][0] - contents[1][0] * contents[2][2]) * inv_det;
        inv_matrix.contents[1][1] = (contents[0][0] * contents[2][2] - contents[0][2] * contents[2][0]) * inv_det;
        inv_matrix.contents[1][2] = (contents[0][2] * contents[1][0] - contents[0][0] * contents[1][2]) * inv_det;

        inv_matrix.contents[2][0] = (contents[1][0] * contents[2][1] - contents[1][1] * contents[2][0]) * inv_det;
        inv_matrix.contents[2][1] = (contents[0][1] * contents[2][0] - contents[0][0] * contents[2][1]) * inv_det;
        inv_matrix.contents[2][2] = (contents[0][0] * contents[1][1] - contents[0][1] * contents[1][0]) * inv_det;

        return inv_matrix;
    }

    void print() {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j)
                std::cout << contents[i][j] << " ";
            std::cout << std::endl;
        }
    }
};

#endif
