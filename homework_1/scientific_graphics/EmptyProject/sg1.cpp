#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

#include "glut.h"



// Precomputed constants:
double circle = atan(1) * 8;
double halfCircle = atan(1) * 4;
double tau = circle; // 2 * PI = TAU
double pi = halfCircle; // TAU / 2 = PI

// How often should the drawing algorithm sample the function.
double step = 0.05;

int defaultW = 1000, defaultH = 1000;
unsigned char prevKey;



// Task 1

double func1 (double x)
{
    if (x == 0)
        return 1.0;

    // distance from x to the nearest integer:
    double floorx = floor(x); // the largest integer less than or equal to x
    double ceilx = ceil(x); // the smallest integer greater than or equal to x
    double d = (x - floorx < ceilx - x) ? x - floorx : ceilx - x; // shortest between the 2 distances (the distance to the closest integer)

    return d / x;
}



void Display3 ()
{
    double xMax = 100.0;
    double yMax = 1.0;  // f(x) is at most 1.

    glColor3f(0, 0, 0);

    glBegin(GL_LINE_STRIP);

    for (double x = 0; x <= xMax; x += step)
    {
        double y = func1(x); // get the y coordinate of the point

        // normalize: x mapped to [0,1], y to [0,1] so it fits into openGL's system
        double nx = x / xMax;
        double ny = y / yMax;

        glVertex2d(nx, ny); // plot the point
    }

    glEnd(); // end drawing
}



// Task 3

void plot (double (*fx)(double, double, double),
    double (*fy)(double, double, double),
    double a, double b,
    double intervalStart, double intervalEnd,
    double plotStep, double scaleX, double scaleY,
    GLint primitive = GL_LINE_STRIP)
{
    glBegin(primitive);

    for (double t = intervalStart; t <= intervalEnd; t += plotStep)
    {
        double x = fx(a, b, t) / scaleX;
        double y = fy(a, b, t) / scaleY;

        glVertex2d(x, y);
    }

    glEnd();
}



// Task 2 - Circle Concoid

double circleConcoidX (double a, double b, double t)
{
    // calculate x
    return 2 * (a * cos(t) + b) * cos(t);
}



double circleConcoidY (double a, double b, double t)
{
    // calculate y
    return 2 * (a * cos(t) + b) * sin(t);
}



void Display4 ()
{
    glColor3f(0, 0, 0);

    plot(circleConcoidX, circleConcoidY, 0.3, 0.2, -pi, pi, step, 1.0, 1.0);
}



// Task 2 - Cicloid

double cicloidX (double a, double b, double t)
{
    // calculate x
    return a * t - b * sin(t);
}



double cicloidY (double a, double b, double t)
{
    // calculate y
    return a - b * cos(t);
}



void Display5 ()
{
    glColor3f(0, 0, 0);

    plot(cicloidX, cicloidY, 0.1, 0.2, -3*pi, 8 * pi, step, 1.0, 0.5);
}



// Task 2 - Epicicloid

double epicicloidX (double a, double b, double t)
{
    // calculate x
    return (a + b) * cos((b / a) * t) - b * cos(t + (b / a) * t);
}



double epicicloidY (double a, double b, double t)
{
    // calculate y
    return (a + b) * sin((b / a) * t) - b * sin(t + (b / a) * t);
}



void Display6 ()
{
    glColor3f(0, 0, 0);

    plot(epicicloidX, epicicloidY, 0.1, 0.3, 0, 2 * pi, step, 0.8, 0.8);
}



// Task 2 - Hipocicloid

double hipocicloidX (double a, double b, double t)
{
    // calculate x
    return (a - b) * cos((b / a) * t) - b * cos(t - (b / a) * t);
}



double hipocicloidY (double a, double b, double t)
{
    // calculate y
    return (a - b) * sin((b / a) * t) - b * sin(t - (b / a) * t);
}



void Display7 ()
{
    glColor3f(0, 0, 0);

    plot(hipocicloidX, hipocicloidY, 0.1, 0.3, 0, 2 * pi, step, 0.8, 0.8);
}



// Task 4 - Logarithmic Spiral

double polarLogSpiralX (double a, double dummy, double t)
{
    double r = a * exp(1 + t);

    return r * cos(t); // get the Cartesian coordinate
}



double polarLogSpiralY (double a, double dummy, double t)
{
    double r = a * exp(1 + t);

    return r * sin(t); // get the Cartesian coordinate
}



void Display8 ()
{
    glColor3f(0, 0, 0);

    plot(polarLogSpiralX, polarLogSpiralY, 0.02, 0, 0, 4, 0.05, 3.0, 3.0);
}



// Task 4 - Sine Polar Flower

double polarFlowerX (double a, double dummy, double t)
{
    double r = sin(a * t);

    return r * cos(t); // get the Cartesian coordinate
}



double polarFlowerY (double a, double dummy, double t)
{
    double r = sin(a * t);

    return r * sin(t); // get the Cartesian coordinate
}



void Display9 ()
{
    glColor3f(0, 0, 0);

    plot(polarFlowerX, polarFlowerY, 10, 0, 0, 2 * pi, 0.01, 1, 1);
}



// Task 5

double trisectrixX (double a, double dummy, double t)
{
    return a / (4 * cos(t) * cos(t) - 3);
}



double trisectrixY (double a, double dummy, double t)
{
    return (a * tan(t)) / (4 * cos(t) * cos(t) - 3);
}



void Display10 ()
{
    std::vector<std::pair<double, double>> points;
    double a = 0.2;
    double tmin = -pi / 2 + 0.01;
    double tmax = pi / 2 - 0.01;
    double singular1 = -pi / 6;
    double singular2 = pi / 6;
    double epsilon = 0.05;

    // Sample: if t is too close to either singularity, skip it.
    for (double t = tmin; t <= tmax; t += step)
    {
        if (fabs(t - singular1) < epsilon || fabs(t - singular2) < epsilon)
            continue;

        double x = trisectrixX(a, 0, t);
        double y = trisectrixY(a, 0, t);
        points.push_back({ x, y });
    }

    // Compute bounding box for scaling:
    double xmin = points[0].first, xmax = points[0].first;
    double ymin = points[0].second, ymax = points[0].second;
    for (auto& pt : points)
    {
        if (pt.first < xmin) xmin = pt.first;
        if (pt.first > xmax) xmax = pt.first;
        if (pt.second < ymin) ymin = pt.second;
        if (pt.second > ymax) ymax = pt.second;
    }

    // Use the maximum absolute value for scaling.
    double scaleX = fabs(xmax) > fabs(xmin) ? fabs(xmax) : fabs(xmin);
    double scaleY = fabs(ymax) > fabs(ymin) ? fabs(ymax) : fabs(ymin);
    if (scaleX == 0) scaleX = 1;
    if (scaleY == 0) scaleY = 1;

    // First, draw the edge (in red) using a line strip.
    glColor3f(1, 0, 0);
    glBegin(GL_LINE_STRIP);
    for (auto& pt : points)
    {
        glVertex2d(pt.first / scaleX, pt.second / scaleY);
    }
    glEnd();

    // Next, fill the interior using a triangle fan.
    // Compute the centroid.
    double cx = 0, cy = 0;
    for (auto& pt : points)
    {
        cx += pt.first;
        cy += pt.second;
    }
    cx /= points.size();
    cy /= points.size();

    glColor3f(0.8, 0.8, 1);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(cx / scaleX, cy / scaleY);
    for (auto& pt : points) {
        glVertex2d(pt.first / scaleX, pt.second / scaleY);
    }
    // Close the fan by connecting back to the first point.
    glVertex2d(points.front().first / scaleX, points.front().second / scaleY);
    glEnd();
}



void Display1 ()
{
    double xmax, ymax, xmin, ymin;
    double a = 1, b = 2;
    double localStep = 0.05;
    xmax = a - b - 1;
    xmin = a + b + 1;
    ymax = ymin = 0;
    for (double t = -pi / 2 + localStep; t < pi / 2; t += localStep) {
        double x1 = a + b * cos(t);
        double x2 = a - b * cos(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;
        xmax = (xmax < x2) ? x2 : xmax;
        xmin = (xmin > x2) ? x2 : xmin;
        double y1 = a * tan(t) + b * sin(t);
        double y2 = a * tan(t) - b * sin(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
        ymax = (ymax < y2) ? y2 : ymax;
        ymin = (ymin > y2) ? y2 : ymin;
    }
    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);
    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);
    for (double t = -pi / 2 + localStep; t < pi / 2; t += localStep) {
        double x1 = (a + b * cos(t)) / xmax;
        double y1 = (a * tan(t) + b * sin(t)) / ymax;
        glVertex2d(x1, y1);
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    for (double t = -pi / 2 + localStep; t < pi / 2; t += localStep) {
        double x2 = (a - b * cos(t)) / xmax;
        double y2 = (a * tan(t) - b * sin(t)) / ymax;
        glVertex2d(x2, y2);
    }
    glEnd();
}



void Display2 ()
{
    double xmax = 8 * pi;
    double ymax = exp(1.1);

    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);
    for (double x = 0; x < xmax; x += step) {
        double x1 = x / xmax;
        double y1 = (fabs(sin(x)) * exp(-sin(x))) / ymax;
        glVertex2d(x1, y1);
    }
    glEnd();
}



//----------------------------------------------------
// GLUT callbacks and main():
void init(void) {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glLineWidth(2);
    glPointSize(1);
    glEnable(GL_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_NICEST, GL_POINT_SMOOTH_HINT);
    glHint(GL_NICEST, GL_LINE_SMOOTH_HINT);
    glHint(GL_NICEST, GL_POLYGON_SMOOTH_HINT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
}



void Display (void)
{
    std::cout << "Call Display" << std::endl;
    glClear(GL_COLOR_BUFFER_BIT);

    // Use the key to decide which plot to display:
    switch (prevKey)
    {
        case '1':
            Display1();
            break;
        case '2':
            Display2();
            break;
        case '3':
            Display3(); // Task 1: function f(x)
            break;
        case '4':
            Display4(); // Circle Concoid
            break;
        case '5':
            Display5(); // Cicloid
            break;
        case '6':
            Display6(); // Epicicloid
            break;
        case '7':
            Display7(); // Hipocicloid
            break;
        case '8':
            Display8(); // Logarithmic spiral
            break;
        case '9':
            Display9(); // Sine polar flower
            break;
        case '0':
            Display10(); // Longchamps' Trisectrix
            break;
        default:
            break;
    }

    glFlush();
}



void Reshape (int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}



void KeyboardFunc (unsigned char key, int x, int y)
{
    prevKey = key;
    if (key == 27) // escape key
        exit(0);
    glutPostRedisplay();
}



void MouseFunc (int button, int state, int x, int y)
{
    std::cout << "Mouse button "
        << ((button == GLUT_LEFT_BUTTON) ? "left" : ((button == GLUT_RIGHT_BUTTON) ? "right" : "middle"))
        << " " << ((state == GLUT_DOWN) ? "pressed" : "released")
        << " at coordinates: " << x << " x " << y << std::endl;
}



int main (int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(defaultW, defaultH);
    glutInitWindowPosition(-1, -1);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutCreateWindow(argv[0]);
    init();

    glutReshapeFunc(Reshape);
    glutKeyboardFunc(KeyboardFunc);
    glutMouseFunc(MouseFunc);
    glutDisplayFunc(Display);

    glutMainLoop();

    return 0;
}