#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <exception>
#include <cmath>
#include <GL/glut.h>


namespace Phi {
    namespace Math {
        constexpr std::size_t newtonSqrtIterations = 4;
        namespace Constant {
            constexpr double grav = 9.80665; // m/s^2
            constexpr double pi = 3.141593; // unitless
            constexpr double stdTemp = 20; // deg C
            constexpr double stdKelvin = 293; // K
            constexpr double stdPressure = 101325; // Pa
            constexpr double Kelvin = 273.15; // K @ 0degC
            constexpr double stdDens = 1.2041; // Kg/m^3
            constexpr double lightspeed = 299792458; // m/s
        }

        double sqrt(const double x) {
            if(!x) return 0;

            double n = 0;
            while(n * n < x) n += 0.1;
            for(size_t i = 0; i < newtonSqrtIterations; ++i) n -= ((n * n) - x) / (2 * n);
            return n;
        }

        class Polynomial {
        private:
            int64_t degree;
            std::vector<double> coefficients;
        public:
            Polynomial(const int64_t degree, const std::vector<double>& coefficients) : degree(degree), coefficients(coefficients) {}

            std::vector<double> solve() {
                std::vector<double> roots{};
                if(this->degree == 1) roots.push_back(-this->coefficients[1] / this->coefficients[0]);
                else if(degree == 2) {
                    roots.push_back((-this->coefficients[1] + Math::sqrt(this->coefficients[1] * this->coefficients[1] - 4 * this->coefficients[0] * this->coefficients[2])) / (2 * this->coefficients[0]) );
                    roots.push_back((-this->coefficients[1] - Math::sqrt(this->coefficients[1] * this->coefficients[1] - 4 * this->coefficients[0] * this->coefficients[2])) / (2 * this->coefficients[0]) );
                }
                return roots;
            }
        };
    };

    template<typename T, std::size_t n>
    class Vector {
    public:
        std::array<T, n> components;
        
        double distance(const Vector<T, n> other) {
            std::array<T, n> resultantLinearDistances = std::array<T, n>();
            for(std::size_t i = 0; i < components.size(); ++i) {
                resultantLinearDistances[i] = other.components[i] - this->components[i];
            }
            double runningTotal = 0;
            std::for_each(this->components.begin(), this->components.end(), [&](const T& i) {
                runningTotal += i * i;
            });
            return Math::sqrt(runningTotal);
        }

        void dump() {
            std::for_each(this->components.begin(), this->components.end(), [&](const T& i) {
                std::cout << i << ' ';
            });
            std::cout << '\n';
        }
    };

    using Vector2D = Vector<double, 2>;

    template<typename T, std::size_t rows, std::size_t columns>
    class Matrix {
    public:
        std::array<double, rows * columns> elements;

        Matrix() {
            elements.fill(0.0);
        }

        Matrix(const unsigned int seed, const T high, const T low) {
            srand(seed);
            std::for_each(this->elements.begin(), this->elements.end(), [&](const T& i) {
                i = T(rand() / RAND_MAX) * (high - low) + low;
            });
        }

        Matrix<T, rows, columns>& operator*(const Matrix<T, rows, columns>&& other) noexcept {
            Matrix<T, rows, columns> resultant = Matrix<T, rows, columns>();
            for(std::size_t i = 0; i < rows; ++i) {
                for(std::size_t j = 0; j < columns; ++j) {
                    const std::size_t index = j * rows + i;
                    resultant[index] = this->elements[index] * other.elements[index];
                }
            }
            return resultant;
        }

        Matrix<T, rows, columns>& operator*(const double&& other) noexcept {
            Matrix<T, rows, columns> resultant = Matrix<T, rows, columns>();
            for(std::size_t i = 0; i < rows; ++i) {
                for(std::size_t j = 0; j < columns; ++j) {
                    const std::size_t index = j * rows + i;
                    resultant[index] = this->elements[index] * other;
                }
            }
            return resultant;
        }

        // Accepts EQs in the form 2x2 and Ans in the form 2x1
        // Returns a 2x1 matrix with x and y solutions
        static Matrix<T, 2, 1> solve(const Matrix<T, 2, 2>& equations, const Matrix<T, 2, 1>& results) {
            Matrix<T, 2, 2> equations_copy = equations;
            // Swap topleft->bottomright negate bottomleft->topright
            equations_copy.elements = { { equations.elements[1][1], -equations.elements[0][1]}, {-equations.elements[1][0], equations.elements[0][0]} };
            equations_copy *= 1 / (equations.elements[1][1] * equations.elements[0][0] - equations.elements[0][1] * equations.elements[1][0]);
            return equations_copy * results;
        }

        void dump() {
            for(std::size_t i = 0; i < this->elements.size(); ++i) {
                std::cout << this->elements[i] << ' ';
                if(!(i % rows)) std::cout << '\n';
            }
        }
    };

    using Matrix2D = Matrix<double, 2, 2>;

    class Shape {
    public:
        enum class ShapeType {
            Circle,
            Rectangle,
            Triangle,
            Semicircle,
            Square,
            Polygon
        } type;

        std::vector<Vector2D> points;

        bool closed = true;

        double width; // m
        union {
            double length; // m
            double radius; // m
        };

        Shape(const double radius) : type(ShapeType::Circle), radius(radius) {
            findAreaPerim();
        }
        Shape(const ShapeType type, const double length, const double width = 0.0) : type(type), length(length), width(width == 0.0 ? length : width) {
            findAreaPerim();
        }
        Shape(const ShapeType type) : type(type) {
            findAreaPerim();
        }


    private:
        double area; // m^2
        double perimeter; // m

        void findAreaPerim() {
            perimeter = 0;
            area = 0;
            switch(this->type) {
                case ShapeType::Polygon: {
                    for(auto i = this->points.begin(); i != this->points.end() - 1; ++i) {
                        const Vector2D next_copy = *(i + 1);
                        perimeter += i->distance(next_copy);
                        area += abs((next_copy.components[1] - i->components[1]) / 2) * (next_copy.components[0] - i->components[0]);
                    }
                    if(this->closed) {
                        area += ((this->points.end()->components[1] - this->points.begin()->components[1]) / 2)
                                * (this->points.end()->components[0] - this->points.begin()->components[0]);
                        perimeter += this->points.begin()->distance(*this->points.end());
                    }
                    break;
                }
                case ShapeType::Circle: {
                    perimeter = 2 * Math::Constant::pi * radius;
                    area = Math::Constant::pi * radius * radius;
                    break;
                }
                case ShapeType::Square:
                case ShapeType::Rectangle: {
                    perimeter = 2 * (length + width);
                    area = length * width;
                    break;
                }
                case ShapeType::Triangle: {
                    std::cout << "case unsupported!  (weirdo)\n";
                    break;
                }
                case ShapeType::Semicircle: {
                    perimeter = 2 * radius + Math::Constant::pi * radius;
                    area = 0.5 * Math::Constant::pi * radius * radius;
                    break;
                }
            }
        }
    };

    namespace PrecomputedDragCoefficient {
        static constexpr double circle = 0.47;
        static constexpr double semicircle = 0.42;
        static constexpr double triangle = 0.5;
        static constexpr double cube = 1.05;
    }

    class Body {
    public:
        Shape shape;

        Vector2D position;
        double rotation = 0; // Radians rotation

        double mass; // kg
        double staticFriction = 0.0; // 0 < u < 1
        double kineticFriction = 0; // 0 < u < 1

        bool fixed = true; // Whether the body can move
        double angularDamp = 0.1; // Constant friction against angular motion  rad/s

        bool rigid = true; // Cannot deform

        Body(const Shape& shape, double mass, double staticFriction = 0.0, double kineticFriction = 0.0, double rotation = 0.0, double radius = 0.0, double length = 0.0, double width = 0.0) : mass(mass), staticFriction(staticFriction), kineticFriction(kineticFriction), rotation(rotation), shape(shape.type == Shape::ShapeType::Circle ? radius : shape.type == Shape::ShapeType::Square ? length : shape.type == Shape::ShapeType::Rectangle ? length, width : shape), position(){

        }
    };
};


void refresh() {

}
void reshape(int width, int height) {

}

int main(int argc, char** argv) {
    std::cout.precision(20);

    Phi::Body kekw = Phi::Body(Phi::Shape(3.0), 1.0);

    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(500, 400);
    glutCreateWindow("Calculator");
    glClearColor(255, 255, 255, 255);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 500, 0.0, 500);
    glutDisplayFunc(refresh);
    glutReshapeFunc(reshape);
    glutMainLoop();

	return 0;
}
