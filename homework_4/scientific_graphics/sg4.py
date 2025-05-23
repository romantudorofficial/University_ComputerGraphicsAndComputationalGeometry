#!/usr/bin/python3
# cython: language_level=3
#This program is (c) Eugen Croitoru, "Al. I. Cuza" University of Iasi, Romania, 2020-2025
#Licensed under GPLv3+
help = '''
if you don't have a C compiler:
pip install PyOpenGL pygame
if you do have a C compiler, accessible system-wide:
pip install PyOpenGL PyOpenGL_accelerate pygame

pyOpenGl reference at:
http://pyopengl.sourceforge.net/documentation/manual-3.0/index.html

pygame reference at e.g.:
https://www.pygame.org/docs/ref/event.html
 
This program works with cython.
How-to: create script named:
setup.py
------------------------------------------
#!/usr/bin/env python3
from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Solver GL',
    ext_modules=cythonize("sg4.pyx"),
    zip_safe=False,
)
------------------------------------------
run the cython-built module:
$ cp sg4.py sg4.pyx -f; python3 setup.py build_ext --inplace; python3 -c 'import sg4; sg4.main()'
'''
print(help)

import sys
import random, time
import math
from math import sin, cos
from ctypes import *
from contextlib import suppress #Amazing! That's why you're the best, Boss, the one and only!

from OpenGL.GL import *
#from OpenGL.GLU import *
#from OpenGL.GLUT import *

import pygame
from pygame.locals import *

def r(radius = 1):
    '''Shorthand function which delivers a random number in [-radius, radius]'''
    return (2 * random.random() - 1) * radius

def randVec(radius = 1, length = 3):
    '''Return a vector (as a list) with the specified length and random values in [-radius, radius]'''
    l = []
    for ii in range(length):
        l.append(r(radius))
    return l


def adjacentFaces(point = [0, 0, 0]):
    '''Get the adjacent faces of a vertex, in a cube.'''
    quad = []
    s = point
    for ii in range(3):
        pl = []
        p = s[:]
        q = s[:]
        for jj in range(3):
            if ii != jj:
                p[jj] = 1 - p[jj]
                q[jj] = p[jj]
                pl.append(p[:])
                p[jj] = 1 - p[jj]
        quad.append(s)
        quad.append(pl[0])
        quad.append(q)
        quad.append(pl[1])
    return(quad)
    
def genCube():
    cube = []
    f = adjacentFaces([0, 0, 0])
    f += adjacentFaces([1, 1, 1])
    for p in f:
        cube.append(p)
    return(cube)

def resizeCube(cube, l = 1, offset = [0, 0, 0]):
    cubeVector = []
    for p in cube:
        for ii in range(3):
            c = offset[ii] + p[ii] * l
            cubeVector.append(c)
    return(cubeVector)
    
def genCol(pointVector, opacity = 1):
    '''Generate the colours for a 3D mesh by simply copying input coordinates.
    For coordinates outside [0, 1], a scaling algorithm should be written.'''
    colVector = []
    for p in pointVector:
        c = p[:]
        c.append(opacity)
        colVector.extend(c)
    return(colVector)

class OglSimpleRenderer():
    '''Requests a drawing surface, and renders and mechanically simulates objects.'''
    def __init__(self, w = 1024, h = 1024, displayFlags = pygame.RESIZABLE | pygame.OPENGL | pygame.DOUBLEBUF, rotation = [-60, 0, 20], scale = 0.5, sceneVolume = [-1, 1, -1, 1, -1, 1]):
        self.w = w
        self.h = h
        self.displayFlags = displayFlags
        self.rotation = rotation
        self.scale = scale
        self.sceneVolume = sceneVolume
        self.viewVolume = self.sceneVolume[:]
        self.viewVolume[-2] *= 2
        self.viewVolume[-1] *= 2
        self.paused = False
        self.shake = False
        pygame.init()
        #init the display size and properties; we specifically request an OpenGL context and a resizable window
        pygame.display.set_mode( (w, h), displayFlags ) # pygame.DOUBLEBUF, pygame.NOFRAME ...
        screen = pygame.display.get_surface()
        if screen.get_flags() & DOUBLEBUF:
            print('Doublebuffered', end = ' ')
        else:
            print('Singlebuffered', end = ' ')
        if screen.get_flags() & OPENGL:
            print('OpenGl')
        else:
            print('No OpenGl')
        #set the window name
        pygame.display.set_caption('SG4 - OpenGL Particle Sim')
        #set the color used to erase the display
        glClearColor(0.0, 0.0, 0.0, 1.0)
        #set the width of lines being drawn
        glLineWidth(3)
        glPointSize(4)
        glEnable(GL_DEPTH_TEST)
        #glEnable(GL_MULTISAMPLE)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        #glEnable(GL_POINT_SMOOTH)
        #glEnable(GL_LINE_SMOOTH)
        #glEnable(GL_POLYGON_SMOOTH)
        #glHint(GL_POINT_SMOOTH, GL_NICEST)
        self.objects = []
        self.clock = pygame.time.Clock()

    def setProjection(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        #Negative z-values wouldn't work for a Frustum, but we're cheating for simplicity -- and because Ortho allows it.
        glOrtho(*self.viewVolume)

    def setModelView(self):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(self.rotation[0], 1, 0, 0)
        glRotatef(self.rotation[1], 0, 1, 0)
        glRotatef(self.rotation[2], 0, 0, 1)
        glScalef(self.scale, self.scale, self.scale)

    def getDt(self):
        '''Get the delta-time between this call and the previous call. Used to compute frame dt.'''
        dtms = self.clock.tick()
        dt = dtms / 1000.0
        return dt
        
    def mouseDown(self, pos, button):
        #Mousewheel up and down 'buttons':
        if button == 4:
            self.scale *= 1.1
        if button == 5:
            self.scale /= 1.1
        
    def mouseMove(self, pos, rel, buttons):
        '''Simpele rotation from mouse drag.'''
        if buttons[2]:
            xrel, yrel = rel
            xrot = xrel / self.w * 180
            yrot = yrel / self.h * 180
            self.rotation[2] += xrot
            self.rotation[0] += yrot

    def keyDown(self, key):
        if key == 'escape':
            sys.exit(0)
        elif key == 'space':
            self.paused = not self.paused
        elif key == 's':
            self.shake = True
        
    def step(self):
        '''An event loop, processing events at every rendering step.'''
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                #self.exitProgram()
                sys.exit(0)
                pass
            elif event.type == pygame.KEYDOWN:
                key = pygame.key.name(event.key)
                self.keyDown(key)
            elif event.type == pygame.KEYUP:
                key = pygame.key.name(event.key)
                #self.keyUp(key)
                pass
            elif event.type == pygame.MOUSEBUTTONDOWN:
                self.mouseDown(event.pos, event.button)
            elif event.type == pygame.MOUSEBUTTONUP:
                #self.mouseUp(event.pos, event.button)
                pass
            elif event.type == pygame.MOUSEMOTION:
                self.mouseMove(event.pos, event.rel, event.buttons)
            elif event.type == pygame.VIDEORESIZE:
                self.resize(event.size)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.display()
        glFlush()
        pygame.display.flip()
        pygame.time.wait(0)

    def loop(self):
        while True:
            self.step()

    def resize(self, size):
        self.w, self.h = size
        print('Window resized to', size)
        glViewport( 0, 0, GLsizei(self.w), GLsizei(self.h) )
        pygame.display.set_mode( (self.w, self.h), self.displayFlags )

    def display(self):
        dt = self.getDt()
        if not self.paused:
            accel = [0, 0, -10]#gravity
            if self.shake:
                self.shake = False
                #The speed we want to add if the accel were to act for a whole second.
                accel[0] += r(10 / dt)
                accel[1] += r(10 / dt)
                #To oppose gravity and down-collision plasticity.
                accel[2] += r(50 / dt)
            for o in self.objects[1:]:  #Quick hack: don't integrate the first object, the larger box.
                self.integrate(o, accel = accel, dt = dt, collisionElasticity = 0.99, collisionFriction = 0.001)
        self.setProjection()
        self.setModelView()
        glMatrixMode(GL_MODELVIEW)
        for o in self.objects:
            self.render(o)

    def integrate(self, obj, accel = [0, 0, 0], dt = 1, collisionElasticity = 0.99, collisionFriction = 0.002):
        '''Integrade the motion of the object, taking into account acceleration, speed, position, and time.
        Also, compute collision with the scene (self.sceneVolume[]).
        collisionElasticity means: how much of the speed is preserved in the direction which caused the collision.
        collisionPlasticity means: how much of the speed is lost in the *other* directions.'''
        _, _, _, _, _, pos, vel, c = obj
        #Homework tasks:
        #1) (1p) Integrate the speed of each object, so that the speed changes position.
        #2) (1p) Also integrate the acceleration, so it influences speed.
        #3) (1p) Do a correct integration with the above, e.g. semi-implicit Euler (also called symplectic Euler)
        #        or the more advanced Runge-Kutta 4. Semi-implicit Euler, for a positive delta-time, requires
        #        speed be updated first, and then position (using the updated speed).
        #4) (1p) Do particle collision with the bounding box.
        #5) (1p) Apply collision elasticity (0.5p) and collision friction (0.5p).
        
    def render(self, obj):
        vertVbo, colVbo, l, cc, primitive, pos, _, _ = obj
        glPushMatrix()
        glTranslate(*pos)
        #activate the drawing of positions
        glEnableClientState(GL_VERTEX_ARRAY)
        #bind the VRAM buffer as active
        glBindBuffer(GL_ARRAY_BUFFER, vertVbo)
        #describe its structure
        glVertexPointer(3, GL_FLOAT, 0, None)

        #same for colours
        glEnableClientState(GL_COLOR_ARRAY)
        glBindBuffer(GL_ARRAY_BUFFER, colVbo)
        #cc is the number of colours, 3 or 4.
        glColorPointer(cc, GL_FLOAT, 0, None)

        #Accelerated draw; the GPU does something like:
        #glBegin(primitive);
        #for(ii = 0, ii < l, ++ii) {
        #  glColour4f(*(colVbo + ii * 4), *(colVbo + ii * 4 + 1)...);
        #  glVertex3f(*(vertVbo + ii * 3), *(vertVbo + ii * 4 + 1)...);
        #glEnd();
        glDrawArrays(primitive, 0, l)
    
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)
        glPopMatrix()

    def addObject(self, vert = [0, 0, 0], col = [1, 1, 1, 1], pos = [0, 0, 0], speed = [0, 0, 0], primitive = GL_QUADS):
        l = len(vert) // 3
        if len(col) // 3 == l:
            cc = 3
        else:
            cc = 4
        #Request a buffer in the VRAM.
        vertVbo = glGenBuffers(1)
        #Bind it as active.
        glBindBuffer(GL_ARRAY_BUFFER, vertVbo)
        #Copy data to it (DYNAMIC_DRAW means the content will not be changed often, but will be
        #read often, so this hints at storing it only on the VRAM. GL_DYNAMIC_DRAW, for instance,
        #would store it both in RAM and VRAM, and monitor the RAM for changes, in a more
        #complex mechanism.
        glBufferData(GL_ARRAY_BUFFER, len(vert) * 4, (c_float * len(vert))(*vert), GL_STATIC_DRAW)
        err = glGetError()
        if err != 0:
            print('Error during buffer allocation', err)
            sys.exit(1)

        colVbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, colVbo)
        glBufferData(GL_ARRAY_BUFFER, len(col) * 4, (c_float * len(col))(*col), GL_STATIC_DRAW)
        err = glGetError()
        if err != 0:
            print('Error during buffer allocation', err)
            sys.exit(1)
        c = [0, 0, 0] #center
        for ii in range(l):
            c[0] += vert[3 * ii]
            c[1] += vert[3 * ii + 1]
            c[2] += vert[3 * ii + 2]
        c = [ c[0] / l, c[1] / l, c[2] / l ]
        self.objects.append( [ vertVbo, colVbo, l, cc, primitive, pos, speed, c ] )



def main():
    print('''Use mouse r-click-drag to rotate, Space to pause-unpause the simulation, and s to shake the box. Mousewheel to zoom.''')
    genCube()
    #half-length of a side
    hl = 10
    speed = hl
    sceneVolume = [-hl, hl, -hl, hl, -hl, hl]
    cubePoints = genCube()
    cubeVert = resizeCube(cubePoints, 2 * hl)
    col = genCol(cubePoints, opacity = 1)
    app = OglSimpleRenderer(sceneVolume = sceneVolume)
    app.addObject(cubeVert, col, pos=[-hl, -hl, -hl], primitive = GL_LINE_STRIP)
    #Add particles, points, which will be rendered as such.
    for ii in range(100):
        app.addObject(pos=randVec(hl), speed = randVec(speed), primitive = GL_POINTS)
    app.loop()

if __name__ == '__main__':
    main()
    
