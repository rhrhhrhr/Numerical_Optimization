//
// Created by 16058 on 2024/5/19.
//

#ifndef NUMERICAL_OPTIMIZATION_MATRIX_HPP
#define NUMERICAL_OPTIMIZATION_MATRIX_HPP

#include <cmath>
#include <cstring>
#include <iostream>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Declaration
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
namespace Matrix {
    template<typename T>
    class Matrix;

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix);

    // 防止编译器把 U 误认为Matrix<T>
    template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
    inline auto operator*(U lhs, const Matrix<T> &rhs) {
        return rhs * lhs;
    }

    template<typename T>
    T norm1(const Matrix<T> &matrix);

    template<typename T>
    auto norm2(const Matrix<T> &matrix);

    template<typename T>
    Matrix<T> identity(size_t n);

    template<typename T>
    class Matrix {
        static_assert(std::is_arithmetic_v<T>, "Matrix data can only be numeric type values.");

        // 便于与不同数据类型的矩阵计算
        template<typename U>
        friend
        class Matrix;

        friend std::ostream &operator
        <<<>(
        std::ostream &out,
        const Matrix<T> &matrix
        );

        friend T norm1<>(const Matrix<T> &matrix);

        friend auto norm2<>(const Matrix<T> &matrix);

        friend Matrix<T> identity<>(size_t n);

    private:
        size_t n_row, n_column, size, capacity;
        T *data;

        T vectorNorm1() const;

        T matrixNorm1() const;

        auto vectorNorm2() const;

        auto matrixNorm2() const;

        void rowExchange(size_t i, size_t j);

        template<typename U>
        void rowScale(size_t i, U scalar);

        template<typename U>
        void rowCombine(size_t i, size_t j, U scalar);

        size_t findMaxAbsInColumnBelowRow(size_t column_index, size_t row_index);

    public:
        Matrix();

        Matrix(size_t n_row, size_t n_column);

        Matrix(size_t n_row, size_t n_column, T value);

        Matrix(const std::initializer_list<T> &list);

        Matrix(const std::initializer_list<std::initializer_list<T>> &list);

        Matrix(const Matrix<T> &other);

        ~Matrix();

        [[nodiscard]] inline bool isVector() const {
            return n_row == 1 or n_column == 1;
        }

        [[nodiscard]] inline size_t getNRow() const {
            return n_row;
        }

        [[nodiscard]] inline size_t getNColumn() const {
            return n_column;
        }

        [[nodiscard]] inline size_t getSize() const {
            return size;
        }

        inline T *operator[](size_t i) {
            return data + i * n_column;
        }

        inline T const *operator[](size_t i) const {
            return data + i * n_column;
        }

        Matrix<T> transpose() const;

        auto inverse() const;

        template<typename U>
        explicit operator Matrix<U>() const;

        Matrix<T> &operator=(const Matrix<T> &other);

        template<typename U>
        auto operator+(const Matrix<U> &other) const;

        template<typename U>
        auto operator-(const Matrix<U> &other) const;

        Matrix<T> operator-() const;

        template<typename U>
        auto operator*(U other) const;

        template<typename U>
        auto operator*(const Matrix<U> &other) const;

        template<typename U>
        inline auto operator/(U other) const {
            return this->operator*(1 / other);
        }
    };
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Implementation
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
namespace Matrix {
    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix) {
        out << "{{";

        for (size_t i = 0; i < matrix.n_row; ++i) {
            if (i > 0) {
                out << " {";
            }

            for (size_t j = 0; j < matrix.n_column; ++j) {
                out << matrix.data[i * matrix.n_column + j];

                if (j != matrix.n_column - 1) {
                    out << ", ";
                }
            }

            if (i != matrix.n_row - 1) {
                out << "}\n";
            }
        }

        out << "}}";

        return out;
    }

    template<typename T>
    T norm1(const Matrix<T> &matrix) {
        return matrix.isVector() ? matrix.vectorNorm1() : matrix.matrixNorm1();
    }

    template<typename T>
    auto norm2(const Matrix<T> &matrix) {
        return matrix.isVector() ? matrix.vectorNorm2() : matrix.matrixNorm2();
    }

    template<typename T>
    Matrix<T> identity(size_t n) {
        Matrix<T> result(n, n);

        for (size_t i = 0; i < n; ++i) {
            result.data[i * result.n_column + i] = 1;
        }

        return result;
    }

    template<typename T>
    T Matrix<T>::vectorNorm1() const {
        T result = 0;

        for (size_t i = 0; i < size; ++i) {
            result += std::abs(data[i]);
        }

        return result;
    }

    template<typename T>
    T Matrix<T>::matrixNorm1() const {
        T result = 0, column_sum;

        for (size_t i = 0; i < n_column; ++i) {
            column_sum = 0;

            for (size_t j = 0; j < n_row; ++j) {
                column_sum += std::abs(data[j * n_column + i]);
            }

            result = (result < column_sum) ? column_sum : result;
        }

        return result;
    }

    template<typename T>
    auto Matrix<T>::vectorNorm2() const {
        auto result = std::pow(data[0], 2);

        for (int i = 1; i < size; ++i) {
            result += std::pow(data[i], 2);
        }

        result = std::sqrt(result);

        return result;
    }

    template<typename T>
    auto Matrix<T>::matrixNorm2() const {
        using U = std::conditional_t<std::is_integral_v<T>, double, T>;

        auto self_transpose_product = (n_row < n_column) ? (*this) * transpose() : transpose() * (*this);
        size_t min_dimension = self_transpose_product.n_row;

        Matrix<U> old_estimated_eigen_vector(min_dimension, 1, 1);
        Matrix<U> estimated_eigen_vector = self_transpose_product * old_estimated_eigen_vector;
        U old_result = estimated_eigen_vector.data[0] / old_estimated_eigen_vector.data[0];
        U result;
        U error;

        do {
            old_estimated_eigen_vector = estimated_eigen_vector;
            estimated_eigen_vector = self_transpose_product * estimated_eigen_vector;

            result = estimated_eigen_vector.data[0] / old_estimated_eigen_vector.data[0];
            error = result - old_result;
            old_result = result;

            estimated_eigen_vector = estimated_eigen_vector / estimated_eigen_vector.vectorNorm2();
        } while (error > 1e-5);

        return std::sqrt(result);
    }

    template<typename T>
    void Matrix<T>::rowExchange(size_t i, size_t j) {
        T temp;
        T *row_i = data + i * n_column, *row_j = data + j * n_column;

        for (size_t k = 0; k < n_column; ++k) {
            temp = row_i[k];
            row_i[k] = row_j[k];
            row_j[k] = temp;
        }
    }

    template<typename T>
    template<typename U>
    void Matrix<T>::rowScale(size_t i, U scalar) {
        T *row_i = data + i * n_column;

        for (int j = 0; j < n_column; ++j) {
            row_i[j] *= scalar;
        }
    }

    template<typename T>
    template<typename U>
    void Matrix<T>::rowCombine(size_t i, size_t j, U scalar) {
        T *row_i = data + i * n_column, *row_j = data + j * n_column;

        for (int k = 0; k < n_column; ++k) {
            row_i[k] += row_j[k] * scalar;
        }
    }

    template<typename T>
    size_t Matrix<T>::findMaxAbsInColumnBelowRow(size_t column_index, size_t row_index) {
        size_t max_abs_value_row_index = row_index;
        T max_abs_value = data[max_abs_value_row_index * n_column + column_index], abs_value_j;

        for (size_t j = row_index + 1; j < n_row; ++j) {
            abs_value_j = std::abs(data[j * n_column + column_index]);
            max_abs_value_row_index = (abs_value_j > max_abs_value) ? j : max_abs_value_row_index;
        }

        return max_abs_value_row_index;
    }

    template<typename T>
    Matrix<T>::Matrix() {
        n_row = 0;
        n_column = 0;
        size = 0;
        capacity = 0;
        data = nullptr;
    }

    template<typename T>
    Matrix<T>::Matrix(size_t n_row, size_t n_column) : n_row(n_row), n_column(n_column) {
        size = n_row * n_column;
        capacity = size;

        if (n_row != 0 and n_column != 0) {
            data = new T[n_row * n_column];
            memset(data, 0, sizeof(T) * n_row * n_column);
        } else {
            data = nullptr;
        }
    }

    template<typename T>
    Matrix<T>::Matrix(size_t n_row, size_t n_column, T value) : n_row(n_row), n_column(n_column) {
        size = n_row * n_column;
        capacity = size;
        data = new T[capacity];

        for (size_t i = 0; i < capacity; ++i) {
            data[i] = value;
        }
    }

    template<typename T>
    Matrix<T>::Matrix(const std::initializer_list<T> &list) {
        n_row = list.size();
        n_column = 1;
        size = n_row;
        capacity = n_row;

        data = new T[capacity];
        size_t i = 0;
        for (auto item: list) {
            data[i] = item;
            ++i;
        }
    }

    template<typename T>
    Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>> &list) {
        n_row = list.size();
        n_column = list.begin()->size();
        size = n_row * n_column;
        capacity = size;

        data = new T[capacity];
        size_t i = 0;
        for (auto sub_list: list) {
            if (sub_list.size() != n_column) {
                throw std::length_error("The input rows do not match!");
            }

            for (auto item: sub_list) {
                data[i] = item;
                ++i;
            }
        }
    }

    template<typename T>
    Matrix<T>::Matrix(const Matrix<T> &other) {
        n_row = other.n_row;
        n_column = other.n_column;
        size = n_row * n_column;
        capacity = size;

        data = new T[capacity];
        for (int i = 0; i < size; ++i) {
            data[i] = other.data[i];
        }
    }

    template<typename T>
    Matrix<T>::~Matrix() {
        delete[] data;
        data = nullptr;
    }

    template<typename T>
    Matrix<T> Matrix<T>::transpose() const {
        Matrix<T> result(n_column, n_row);

        for (size_t i = 0; i < n_row; ++i) {
            for (int j = 0; j < n_column; ++j) {
                result.data[j * n_row + i] = data[i * n_column + j];
            }
        }

        return result;
    }

    template<typename T>
    auto Matrix<T>::inverse() const {
        if (n_row != n_column) {
            throw std::length_error("Only square matrix has inverse!");
        }

        using U = std::conditional_t<std::is_integral_v<T>, double, T>;

        auto copy = static_cast<Matrix<U>>(*this);
        Matrix<U> result = identity < U > (n_row);

        size_t max_abs_index;
        U max_abs_value;
        U k;

        for (int i = 0; i < n_column; ++i) {
            max_abs_index = copy.findMaxAbsInColumnBelowRow(i, i);
            max_abs_value = std::abs(copy.data[max_abs_index * n_column + i]);

            if (max_abs_value == 0) {
                throw std::domain_error("The matrix is singular!");
            }

            if (max_abs_index != i) {
                copy.rowExchange(max_abs_index, i);
                result.rowExchange(max_abs_index, i);
            }

            k = 1. / max_abs_value;

            copy.rowScale(i, k);
            result.rowScale(i, k);

            for (int j = 0; j < n_row; ++j) {
                if (j != i) {
                    k = -copy.data[j * n_column + i];

                    copy.rowCombine(j, i, k);
                    result.rowCombine(j, i, k);
                }
            }
        }

        return result;
    }

    template<typename T>
    template<typename U>
    Matrix<T>::operator Matrix<U>() const {
        Matrix<U> result;

        result.n_row = n_row;
        result.n_column = n_column;
        result.size = size;
        result.capacity = size;

        result.data = new U[size];
        for (int i = 0; i < size; ++i) {
            result.data[i] = static_cast<U>(data[i]);
        }

        return result;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other) {
        if (this != &other) {
            n_row = other.n_row;
            n_column = other.n_column;

            if (capacity < other.size) {
                capacity = other.size;
                delete[] data;
                data = new T[other.size];
            }

            memcpy(data, other.data, sizeof(T) * other.size);
        }

        return *this;
    }

    template<typename T>
    template<typename U>
    auto Matrix<T>::operator+(const Matrix<U> &other) const {
        if (n_row != other.n_row or n_column != other.n_column) {
            throw std::length_error("");
        }

        Matrix<decltype(T() + U())> result(n_row, n_column);

        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] + other.data[i];
        }

        return result;
    }

    template<typename T>
    template<typename U>
    auto Matrix<T>::operator-(const Matrix<U> &other) const {
        if (n_row != other.n_row or n_column != other.n_column) {
            throw std::length_error("");
        }

        Matrix<decltype(T() - U())> result(n_row, n_column);

        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] - other.data[i];
        }

        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator-() const {
        Matrix<T> result(n_row, n_column);

        for (int i = 0; i < size; ++i) {
            result.data[i] = -data[i];
        }

        return result;
    }

    template<typename T>
    template<typename U>
    auto Matrix<T>::operator*(U other) const {
        Matrix<decltype(T() * U())> result(n_row, n_column);

        for (size_t i = 0; i < size; ++i) {
            result.data[i] = data[i] * other;
        }

        return result;
    }

    template<typename T>
    template<typename U>
    auto Matrix<T>::operator*(const Matrix<U> &other) const {
        if (n_column != other.n_row) {
            throw std::length_error("");
        }

        Matrix<decltype(T() * U())> result(n_row, other.n_column);
        T temp;

        for (size_t i = 0; i < n_row; ++i) {
            for (size_t k = 0; k < n_column; ++k) {
                temp = data[i * n_column + k];
                for (size_t j = 0; j < other.n_column; ++j) {
                    result.data[i * other.n_column + j] += temp * other.data[k * other.n_column + j];
                }
            }
        }

        return result;
    }
}

#endif //NUMERICAL_OPTIMIZATION_MATRIX_HPP
