#include <iostream>
#include <GL/glut.h>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;

class body {
public:
    double mass; // kg
    enum class shapes {
        CIRCLE,
        RECTANGLE,
        TRIANGLE,
        SEMICIRCLE,
        SQUARE,
        POLYGON,
    } shape;
    
    struct form {
        vector<vector<double>> Points;
        bool close; // connect the first and last points
    } Form;
    double frictionS; // static friction 0 < u < 1
    double frictionK; // kinetic friction 0 < u < 1
    double perimeter; // m
    double area; // m^2
    double radius; // m
    double length; // rect//m
    double width; // rect//m
    bool rigid = true; // cannot deform
    bool Static; // cannot move
    double X; // x transform
    double Y; // y transform
    double rotation; // radians rotation
    double angularDamp; // constant friction against angular motion


private:
    void findAreaPerim(body X) {
        
    }
    vector<double> dragLookup = { 0.47, 0.42, 0.5, 1.05 }; // circle, semiC, triangle, cube
};
class Math {
public:
    double sqrt(double x) {
        if (x == 0) {
            return 0;
        }
        double n = 0;
        while (n * n < x) {
            n += 0.1;
        }
        n = n - ((n * n) - x) / (2 * n);
        n = n - ((n * n) - x) / (2 * n);
        n = n - ((n * n) - x) / (2 * n);
        n = n - ((n * n) - x) / (2 * n);
        return n;
    }
};
class Poly {
public:
    Math mathENG;
    vector<double> Solve(Poly nomial) {
        vector<double> Answer;
        if (degree == 1) {
            Answer.push_back(-nomial.poly[1] / nomial.poly[0]);
        }
        else if (degree == 2) {
            Answer.push_back( (-nomial.poly[1] + mathENG.sqrt( nomial.poly[1] * nomial.poly[1] - 4 * nomial.poly[0] * nomial.poly[2])) / (2*nomial.poly[0]) );
            Answer.push_back( (-nomial.poly[1] - mathENG.sqrt( nomial.poly[1] * nomial.poly[1] - 4 * nomial.poly[0] * nomial.poly[2])) / (2 * nomial.poly[0]) );
        }
        return Answer;
    }
    Poly(int _degree, vector<double> _poly) {
        degree = _degree;
        poly = _poly;
    }
    int degree;
    vector<double> poly;
};
class Matrix {
// Matrix object, n*m dimensional, with random init, clearing, mult, add, sub   - double only
public:
    int rows;
    int cols;
    vector<vector<double>> matrix; // main matrix object, two dimensional vector double

    // class initialiser, takes rows and cols and fills with 0's
    // row is the number of horizontals
    Matrix(int _rows, int _cols) {
        rows = _rows;
        cols = _cols;
        vector<double> tempvect;
        for (int i = 0; i < cols; i++) {
            tempvect.push_back(0);
        }

        for (int i = 0; i < rows; i++) {
            matrix.push_back(tempvect);
        }
    }

    // prints out the matrix with rows and cols organised as such
    void Print() {
        for (int i = 0; i < rows; i++) {
            cout << "\nrow: ";
            for (int j = 0; j < cols; j++) {
                cout << matrix[i][j] << " ";
            }
        }
    }

    // randomly fills the matrix with numbers from low to high with 6 significant figures, int or dec
    void randomfill(int low, int high) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = (double(rand() % 1000000) / 1000000) * (high - low) + low;
            }
        }
    }

    // multiplies matrix Mult1 by Mult2 and returns ans matrix without modifying originals
    Matrix mult(Matrix Mult1, Matrix Mult2) {
        Matrix answer = Matrix(Mult1.rows, Mult2.cols);
        for (int i = 0; i < Mult1.rows; i++) { // FIXXXXXXXXXXXXXX
            for (int j = 0; j < Mult2.cols; j++) {
                for (int k = 0; k < Mult2.rows && k < Mult1.cols; k++) {
                    answer.matrix[i][j] += Mult1.matrix[i][k] * Mult2.matrix[k][j];
                }
            }
        }
        return answer;
    }
    Matrix mult(Matrix Mult1, double Mult2) {
        Matrix answer = Mult1;
        for (int i = 0; i < Mult1.rows; i++) {
            for (int j = 0; j < Mult1.cols; j++) {
                answer.matrix[i][j] *= Mult2;
            }
        }
        return answer;
    }

    // accepts EQs in the form 2x2 and Ans in the form 2x1
    // returns a 2x1 matrix with x and y solutions
    Matrix Solve(Matrix EQs, Matrix Ans) {
        Matrix solved = Matrix(2, 1);
        if (EQs.rows == 2 && EQs.cols == 2 && Ans.rows == 2 && Ans.cols == 1) {
            Matrix TempEQ = EQs;
            TempEQ.matrix = { {EQs.matrix[1][1], -EQs.matrix[0][1]}, {-EQs.matrix[1][0], EQs.matrix[0][0]} };
            TempEQ = TempEQ.mult(TempEQ, 1 / (EQs.matrix[1][1] * EQs.matrix[0][0] - EQs.matrix[0][1] * EQs.matrix[1][0]));
            solved = TempEQ.mult(TempEQ, Ans);
        }
        return solved;
    }

private:

};
class Constant {
    public:
        double grav = 9.80665; // m/s^2
        double pi = 3.141593; // unitless
        double stdTemp = 20; // deg C
        double stdKelvin = 293; // K
        double stdPressure = 101325; // Pa
        double Kelvin = 273.15; // K @ 0degC
        double stdDens = 1.2041; // Kg/m^3
        double lightspeed = 299792458; // m/s


    private:

};
class Graphics {
    public:

};
void DisplayUpdated() {

}
void reshape(int width, int height) {

}
void VectPrint(vector<double> vect) {
    for (int i = 0; i < vect.size(); i++) {
        cout << vect[i] << "\n";
    }
}


int main(int argc, char** argv) {
    cout.precision(20);

    body kekw = body();
    kekw.shape = body::shapes::CIRCLE;
    switch (kekw.shape) {
        case body::shapes::CIRCLE: {
            cout << "kekw";
            break;
        };
        
    }
    Constant Const;
    Graphics graphic;
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(500, 400);
    glutCreateWindow("Calculator");
    glClearColor(255, 255, 255, 255);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 500, 0.0, 500);
    glutDisplayFunc(DisplayUpdated);
    glutReshapeFunc(reshape);
    glutMainLoop();
	return 0;
}