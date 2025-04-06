#include "glut.h"
#include <cmath>
#include <algorithm>

// -----------------------------------------------------
// Global Parameters (Single coordinate system: grid intersections at integers)
// -----------------------------------------------------

int windowWidth = 800;
int windowHeight = 800;

int gridSize = 20;  // Grid from -gridSize to +gridSize (in integer units)
float margin = 1.0f;  // Margin (in the same units)

// Disc radius (in world units) for Bresenham "pixel" discs; increased for larger dots.
float bresenhamDiscRadius = 0.3f;

// -----------------------------------------------------
// Disc Drawing Function
// -----------------------------------------------------

// Draw a filled disc centered at (x, y) with the given radius using GL_POLYGON.
void drawDisc(int x, int y, float radius) {
    glBegin(GL_POLYGON);
    for (int i = 0; i < 360; i++) {
        float theta = 2.0f * 3.1415926f * float(i) / 360.0f;
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex2f(x + dx, y + dy);
    }
    glEnd();
}

// -----------------------------------------------------
// Grid Drawing
// -----------------------------------------------------

// Draw grid lines (light gray) at every integer coordinate.
// No discs are drawn at grid intersections.
void drawGrid() {
    glColor3f(0.8f, 0.8f, 0.8f);
    // Draw vertical grid lines.
    for (int i = -gridSize; i <= gridSize; i++) {
        glBegin(GL_LINES);
        glVertex2f(i, -gridSize);
        glVertex2f(i, gridSize);
        glEnd();
    }
    // Draw horizontal grid lines.
    for (int j = -gridSize; j <= gridSize; j++) {
        glBegin(GL_LINES);
        glVertex2f(-gridSize, j);
        glVertex2f(gridSize, j);
        glEnd();
    }
}

// -----------------------------------------------------
// Blue Circle Drawing (as an Octagon)
// -----------------------------------------------------

// Draw a blue octagon (an approximation of a circle) centered at (cx,cy)
// with the given radius. Vertices are computed in integer coordinates.
void drawBlueCircleOutline(int cx, int cy, int radius) {
    glColor3f(0.0f, 0.0f, 1.0f);
    const int numSides = 8;
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numSides; i++) {
        float theta = 2.0f * 3.1415926f * float(i) / numSides;
        int vx = cx + static_cast<int>(round(radius * cosf(theta)));
        int vy = cy + static_cast<int>(round(radius * sinf(theta)));
        glVertex2f(vx, vy);
    }
    glEnd();
}

// -----------------------------------------------------
// Bresenham’s Mid‑Point Line Rasterisation Algorithm
// -----------------------------------------------------

// The following four functions implement Bresenham's algorithm
// for different octants. All loop computations use only additions/subtractions.

// Octant 1: slope between 0 and 1 (increasing x and y).
void bresenhamOctant1(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0; // 0 <= dy <= dx
    int d = 2 * dy - dx;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        drawDisc(x, y, bresenhamDiscRadius);
        if (d > 0) {
            y++;
            d = d + 2 * (dy - dx);
        }
        else {
            d = d + 2 * dy;
        }
    }
}

// Octant 2: slope greater than 1 (increasing x and y; iterate over y).
void bresenhamOctant2(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0; // dy > dx
    int d = 2 * dx - dy;
    int x = x0;
    for (int y = y0; y <= y1; y++) {
        drawDisc(x, y, bresenhamDiscRadius);
        if (d > 0) {
            x++;
            d = d + 2 * (dx - dy);
        }
        else {
            d = d + 2 * dx;
        }
    }
}

// Octant 8: slope between 0 and -1 (increasing x, decreasing y).
void bresenhamOctant8(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0; // dy is negative; use |dy| = -dy
    int d = 2 * (-dy) - dx;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        drawDisc(x, y, bresenhamDiscRadius);
        if (d > 0) {
            y--;
            d = d + 2 * ((-dy) - dx);
        }
        else {
            d = d + 2 * (-dy);
        }
    }
}

// Octant 7: slope less than -1 (increasing x, decreasing y; iterate over y).
void bresenhamOctant7(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0; // dy is negative and |dy| > dx
    int d = 2 * dx - (-dy);
    int x = x0;
    for (int y = y0; y >= y1; y--) {
        drawDisc(x, y, bresenhamDiscRadius);
        if (d > 0) {
            x++;
            d = d + 2 * (dx - (-dy));
        }
        else {
            d = d + 2 * dx;
        }
    }
}

// Wrapper that selects the proper sub‑algorithm based on the slope.
// It also ensures drawing from left-to-right by swapping endpoints if necessary.
void drawBresenhamLine(int x0, int y0, int x1, int y1) {
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    if (dy >= 0) {
        if (dx >= dy) {
            bresenhamOctant1(x0, y0, x1, y1);
        }
        else {
            bresenhamOctant2(x0, y0, x1, y1);
        }
    }
    else { // dy < 0
        if (dx >= -dy) {
            bresenhamOctant8(x0, y0, x1, y1);
        }
        else {
            bresenhamOctant7(x0, y0, x1, y1);
        }
    }
}

// -----------------------------------------------------
// Combined Drawing: Blue Circle with Bresenham Rasterisation
// -----------------------------------------------------

// This function first draws the blue circle outline (an octagon) and then,
// for each edge of the octagon, applies Bresenham's algorithm to "light up"
// the grid intersections (by drawing black discs) that approximate that edge.
void drawCircleWithRasterization(int cx, int cy, int radius) {
    // Draw the blue outline.
    drawBlueCircleOutline(cx, cy, radius);

    // Compute the octagon's vertices.
    const int numSides = 8;
    int verticesX[numSides], verticesY[numSides];
    for (int i = 0; i < numSides; i++) {
        float theta = 2.0f * 3.1415926f * i / numSides;
        verticesX[i] = cx + static_cast<int>(round(radius * cosf(theta)));
        verticesY[i] = cy + static_cast<int>(round(radius * sinf(theta)));
    }
    // For each edge, run Bresenham's algorithm.
    for (int i = 0; i < numSides; i++) {
        int next = (i + 1) % numSides;
        drawBresenhamLine(verticesX[i], verticesY[i], verticesX[next], verticesY[next]);
    }
}

// -----------------------------------------------------
// Window Reshape Callback
// -----------------------------------------------------

// Set up an orthographic projection so that grid intersections at integer coordinates are shown.
void reshape(int w, int h) {
    if (h == 0)
        h = 1;
    windowWidth = w;
    windowHeight = h;
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float left = -gridSize - margin;
    float right = gridSize + margin;
    float bottom = -gridSize - margin;
    float top = gridSize + margin;

    float aspect = float(w) / float(h);
    if (aspect > 1.0f) {
        float extra = (right - left) * (aspect - 1.0f) * 0.5f;
        left -= extra;
        right += extra;
    }
    else {
        float extra = (top - bottom) * (1.0f / aspect - 1.0f) * 0.5f;
        bottom -= extra;
        top += extra;
    }
    glOrtho(left, right, bottom, top, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
}

// -----------------------------------------------------
// Display Callback
// -----------------------------------------------------

// Clear the screen, draw the grid, and then draw the blue circle (as an octagon)
// with its edge pixels "rasterized" by Bresenham's algorithm.
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    drawGrid();
    // Draw a circle centered at (0,0) with radius 10 (in grid units).
    drawCircleWithRasterization(0, 0, 10);
    glutSwapBuffers();
}

// -----------------------------------------------------
// Initialization
// -----------------------------------------------------

void init() {
    glClearColor(1.0, 1.0, 1.0, 1.0); // White background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float left = -gridSize - margin;
    float right = gridSize + margin;
    float bottom = -gridSize - margin;
    float top = gridSize + margin;
    gluOrtho2D(left, right, bottom, top);
}

// -----------------------------------------------------
// Main Function
// -----------------------------------------------------

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Blue Circle with Bresenham Rasterisation on Grid");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();
    return 0;
}
