//
// Created by 16058 on 2024/6/2.
//

#ifndef NUMERICAL_OPTIMIZATION_EXPRESSION_HPP
#define NUMERICAL_OPTIMIZATION_EXPRESSION_HPP

#include <iostream>

namespace Model {
    template<typename T>
    class Expression;

    template<typename T>
    struct is_expression : std::false_type {};

    template<typename T>
    struct is_expression<Expression<T>> : std::true_type {};

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Expression<T> &expression);

    template<typename T>
    constexpr bool is_expression_v = is_expression<T>::value;

    template<typename T>
    class Expression {
        friend std::ostream &operator
        <<<>(
        std::ostream &out,
        const Expression<T> &expression
        );
    private:
        T value;
    public:
        Expression() : value(0) {}

        explicit Expression(T value) : value(value) {};

        Expression(const Expression<T> &other) = default;

        ~Expression() = default;

        Expression<T> &operator=(const Expression<T> &other) = default;

        template<typename U>
        Expression<T> &operator=(U other);

        inline T getValue() const {
            return value;
        }
    };
}

namespace Model {
    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Expression<T> &expression) {
        out << expression.value;
        return out;
    }

    template<typename T>
    template<typename U>
    Expression<T> &Expression<T>::operator=(U other) {
        this->value = static_cast<T>(other);
        return *this;
    }
}

#endif //NUMERICAL_OPTIMIZATION_EXPRESSION_HPP
