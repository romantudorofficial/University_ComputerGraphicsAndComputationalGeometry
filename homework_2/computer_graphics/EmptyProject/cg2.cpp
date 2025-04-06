#include "glut.h"
#include <cmath>

// Window size
int windowWidth = 800;
int windowHeight = 800;

// Grid properties
int gridSize = 20;
float cellSize = 0.05f;  // Size of each grid cell in normalized coordinates

// Function to draw a disc at a given integer grid position
void drawDisc(int x, int y, float radius = 0.02f) {
    glBegin(GL_POLYGON);
    for (int i = 0; i < 360; i++) {
        float theta = 2.0f * 3.1415926f * float(i) / 360.0f;
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex2f(x + dx, y + dy);
    }
    glEnd();
}

// Function to draw the raster grid with discs at intersections
void drawGrid() {
    glColor3f(0.8f, 0.8f, 0.8f);  // Light gray for the grid lines

    // Draw grid lines
    for (int i = -gridSize; i <= gridSize; i++) {
        glBegin(GL_LINES);
        glVertex2f(i * cellSize, -gridSize * cellSize);
        glVertex2f(i * cellSize, gridSize * cellSize);
        glEnd();

        glBegin(GL_LINES);
        glVertex2f(-gridSize * cellSize, i * cellSize);
        glVertex2f(gridSize * cellSize, i * cellSize);
        glEnd();
    }

    // Draw discs at the intersections of grid lines
    glColor3f(0.0f, 0.0f, 0.0f);  // Black for the discs
    for (int i = -gridSize; i <= gridSize; i++) {
        for (int j = -gridSize; j <= gridSize; j++) {
            drawDisc(i, j);
        }
    }
}

// Function to draw a geometric line segment with a thickness greater than 1 disc
void drawLine(float x1, float y1, float x2, float y2, float thickness = 0.04f) {
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue for the line

    // Draw a thick line as a collection of discs along the path
    int segments = 20;
    for (int i = 0; i < segments; i++) {
        float t = (float)i / (segments - 1);
        float x = (1 - t) * x1 + t * x2;
        float y = (1 - t) * y1 + t * y2;
        drawDisc(int(x), int(y), thickness);
    }
}

// Function to draw a circle made from straight lines (octagon approximation)
void drawCircle(float x, float y, float radius = 0.1f) {
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue for the circle

    int numSides = 8;  // Number of segments for the octagon approximation
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numSides; i++) {
        float theta = 2.0f * 3.1415926f * float(i) / float(numSides);
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex2f(x + dx, y + dy);
    }
    glEnd();
}

// Function to adjust the orthographic projection on window resize
void reshape(int w, int h) {
    // Prevent division by zero
    if (h == 0) h = 1;

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Set the aspect ratio to ensure square grid cells, regardless of window size
    if (w > h) {
        glOrtho(-gridSize * cellSize * float(w) / float(h), gridSize * cellSize * float(w) / float(h), -gridSize * cellSize, gridSize * cellSize, -1.0, 1.0);
    }
    else {
        glOrtho(-gridSize * cellSize, gridSize * cellSize, -gridSize * cellSize * float(h) / float(w), gridSize * cellSize * float(h) / float(w), -1.0, 1.0);
    }

    glMatrixMode(GL_MODELVIEW);
}

// Display function
void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw grid and other objects
    drawGrid();
    drawLine(-2.0f, -1.0f, 2.0f, 1.0f);  // Example line
    drawCircle(1.5f, 1.5f);  // Example circle

    glutSwapBuffers();
}

// Initialization
void init() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);  // White background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-gridSize * cellSize, gridSize * cellSize, -gridSize * cellSize, gridSize * cellSize);
}

// Main function
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Raster Grid and Geometric Primitives");

    init();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);  // Register reshape function

    glutMainLoop();

    return 0;
}
