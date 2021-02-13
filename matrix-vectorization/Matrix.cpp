
#include <cstdlib>
#include <random>
#include <iostream>
#include <iomanip>
#include <immintrin.h>
#include "Matrix.h"

Matrix::Matrix(const Matrix& matrix) : rows(matrix.rows), cols(matrix.cols) {
    data = new float* [rows];
    
    for (int i = 0; i < rows; ++i) {
        data[i] = new float[cols];
        for (int j = 0; j < cols; ++j) {
            data[i][j] = matrix.getData()[i][j];
        }
    }

}

Matrix::~Matrix() {
    Matrix::clear();
}

Matrix Matrix::multiply(Matrix& m1, Matrix& m2)
{
    auto data = m1.data;
    auto data_ = m2.getData();
    auto cols = m1.cols;
    auto cols_ = m2.getCols();
    auto rows = m1.rows;
    auto resultData = allocate(rows, cols_);

    for (int i = 0; i < MT_M; ++i) {
        auto resultRow = resultData[i];
        auto row = data[i];
        for (int j = 0; j < MT_N; ++j) {
            auto val = row[j];
            auto col = data_[j];
#pragma loop(no_vector) 
            for (int k = 0; k < cols_; ++k) {
                resultRow[k] += val * col[k];
            }
        }
    }

    Matrix result(rows, cols_);
    result.fill(resultData);

    return result;
}

Matrix Matrix::multiplyOptimized(Matrix& m1, Matrix& m2)
{
    if (m1.cols != m2.rows)
        throw std::runtime_error("Can't multiply matrix");

    auto data = m1.getData();
    auto data_ = m2.getData();

    auto resultData = allocate(MT_M, MT_K);

    for (int i = 0; i < MT_M; ++i) {
        auto resultRow = resultData[i];
        auto row = data[i];
        for (int j = 0; j < MT_N; ++j) {
            auto val = row[j];
            __m256 row1 = _mm256_set_ps(val, val, val, val, val, val, val, val);
            for (int k = 0; k < MT_K; k+=8) {
                auto col = data_[j] + k;
                //std::cout << val << "_____" <<  j  << "______" << + k << "____" << col[0] <<  std::endl;
                //std::cout << "____" << resultRow[k] << std::endl;
                __m256 row2 = _mm256_load_ps(col);
                __m256 res_row = _mm256_load_ps(resultRow + k);
                __m256 fmadd = _mm256_fmadd_ps(row1, row2, res_row);
              
                _mm256_store_ps(resultRow + k, fmadd);
                //std::cout << resultRow[k] << std::endl;
                //float* ptr = (float*)&row2;
                //std::cout << ptr[0] << std::endl;
                //printf("%f %f %f %f %f %f %f %f\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
            }
        }
    }

    Matrix result(MT_M, MT_K);

    result.fill(resultData);
    return result;
}

float **Matrix::allocate(int rows, int cols) {
    auto data_ = new float *[rows];
    for (int i = 0; i < rows; ++i) {
        data_[i] = new float[cols];
        for (int j = 0; j < cols; ++j)
            data_[i][j] = 0;
    }
    return data_;
}

void Matrix::fill() {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] = getRandomFloat(MIN_VALUE, MAX_VALUE);
        }
    }
}

float Matrix::getRandomFloat(float min, float max) {
    
    float a = min + (float) (rand()) / ((float) (RAND_MAX / (max - min)));
    return a;
}

void Matrix::fill(float value) {
    for (int i = 0; i < rows ; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] = value;
        }
    }
}

void Matrix::clear() {
    for (int i = 0; i < rows; ++i) {
        delete[] data[i];
    }
    delete data;
}

const void Matrix::print() {
   ;
    for (int i = 0; i < rows; i++) {
        std::cout << "[";
        for (int j = 0; j < cols; j++) {
            std::cout << std::fixed << std::setprecision(2) << std::setw(5)  << data[i][j] << ' ';
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void Matrix::fill(float **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            data[i][j] = matrix[i][j];
        }
    }
}

Matrix Matrix::operator*(Matrix &matrix) {
    if (cols != matrix.getRows())
        throw std::runtime_error("Can't multiply matrix");

    auto data_ = matrix.getData();

    auto resultData = allocate(MT_M, MT_K);

    for (int i = 0; i < MT_M; ++i) {
        auto resultRow = resultData[i];
        auto row = data[i];
        for(int j = 0; j < MT_N; ++j){
            auto val = row[j];
            auto col = data_[j];
#pragma loop(no_vector) 
            for(int k = 0; k < MT_K; ++k){
                resultRow[k] += val * col[k];
            }
        }
    }

    Matrix result(MT_M, MT_K);
    
    result.fill(resultData);
    return result;

}

Matrix Matrix::operator+(Matrix &matrix1) {
    auto resultData = allocate(this->rows, this->cols);
    auto data1 = matrix1.getData();
    auto data2=data;
    for (int i=0; i<MT_M; ++i) {
        auto row1=data1[i];
        auto row2=data2[i];
#pragma loop(no_vector) 
        for (int j=0; j< MT_K; ++j) {
            resultData[i][j]=row1[j]+row2[j];
        }
    }

    Matrix result(rows, cols);
    result.fill(resultData);

    return result;
}

Matrix& Matrix::operator=(const Matrix &matrix) {
    this->cols=matrix.cols;
    this->rows=matrix.rows;
    this->data=allocate(matrix.rows, matrix.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            this->data[i][j] = matrix.getData()[i][j];
        }
    }
    return *this;
}

void Matrix::add(Matrix& m)
{
    auto data1 = m.getData();
    for (int i = 0; i < MT_M; ++i) {
        auto row1 = data1[i];
#pragma loop(no_vector) 
        for (int j = 0; j < MT_K; ++j) {
            data[i][j] += row1[j];
        }
    }
}

