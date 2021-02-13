#include <iostream>
#include "Matrix.h"



int main() {

    const int M = MT_ROWS;
    const int N = MT_COLS;

    auto **pMatrix1 = new Matrix *[MT_ROWS];
    for (int i = 0; i < MT_COLS; ++i) {
        pMatrix1[i] = new Matrix[N];
        for (int j = 0; j < N; ++j) {
            Matrix matrix(M, N);
            matrix.fill(2);
            pMatrix1[i][j]=matrix;
        }
    }

    const int K = MT_ROWS;
    auto **pMatrix2 = new Matrix *[N];
    for (int i = 0; i < N; ++i) {
        pMatrix2[i] = new Matrix[K];
        for (int j = 0; j < K; ++j) {
            Matrix matrix(N, K);
            matrix.fill(2);

            pMatrix2[i][j]=matrix;
        }
    }

    auto **result = new Matrix *[M];
    for (int i = 0; i < M; ++i) {
        result[i] = new Matrix[K];
        for (int j = 0; j < K; ++j) {
            Matrix matrix(M, K);
            matrix.fill(0.);
            result[i][j]=matrix;
        }
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < K; ++k) {
                Matrix m  = Matrix::multiply(pMatrix1[i][j], pMatrix2[j][k]);
                result[i][k].add(m);
                if (k==K-1) {
                    result[i][k].print();
                }
            }
        }
    }



}
