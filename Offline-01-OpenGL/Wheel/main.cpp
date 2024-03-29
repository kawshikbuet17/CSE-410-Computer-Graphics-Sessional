#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/freeglut.h>
#include <bits/stdc++.h>
using namespace std;

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
double wheelRadius;
double wheelHeight;
double wheelAngle;
double rimAngle;
double angleChange;
double distanceChange;
int drawgrid;
int drawaxes;
double angle;

struct point
{
    double x,y,z;
};

struct point wheelCenter;

void initWheelCenter(){
    wheelCenter.x = 0;
    wheelCenter.y = 0;
    wheelCenter.z = 0;
}

void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);{
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);

            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);

            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }glEnd();
    }
}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);{
            for(i=-20;i<=20;i++){
                //lines parallel to Y-axis
                glVertex3f(i*10, -200, 0);
                glVertex3f(i*10,  200, 0);

                //lines parallel to X-axis
                glVertex3f(-200, i*10, 0);
                glVertex3f( 200, i*10, 0);
            }
        }glEnd();
    }
}

void drawRectangle(double a, double b)
{
    glColor3f(0.6,0.6,0.6);
    glBegin(GL_QUADS);{
        glVertex3f( a, b,0);
        glVertex3f( a,-b,0);
        glVertex3f(-a,-b,0);
        glVertex3f(-a, b,0);
    }glEnd();
}


void drawCylinder(double radius, double height, int slices,int stacks)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0;i<=stacks;i++)
    {
        h=height * (double)i/(double)stacks;
        r=radius;
        for(j=0;j<=slices;j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        double red = (double)i/(double)stacks;
        double green = (double)i/(double)stacks;
        double blue = (double)i/(double)stacks;
        glColor3f(red, green, blue);
        for(j=0;j<slices;j++)
        {
            glBegin(GL_QUADS);{
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }glEnd();
        }
    }
}

void drawWheel(){
    double radius = wheelRadius;
    double height = wheelHeight;
    glTranslatef(wheelCenter.x,wheelCenter.y,wheelCenter.z+wheelRadius);
    glRotatef(wheelAngle, 0, 0, 1);
    glRotatef(rimAngle, 0, 1, 0);
    glRotatef(90, 1, 0, 0);
    drawCylinder(radius, height, 20, 20);
    glRotatef(90, 1, 0, 0);
    drawRectangle(radius, height);
    glRotatef(90, 0, 1, 0);
    drawRectangle(radius, height);
}

void keyboardListener(unsigned char key, int x,int y){
    switch(key){

        case '1':
            drawgrid=1-drawgrid;
            break;
        case 'a':
            wheelAngle -= angleChange;
            break;
        case 'd':
            wheelAngle += angleChange;
            break;
        case 'w':
            wheelCenter.x += distanceChange* cos(wheelAngle * pi / 180);
            wheelCenter.y += distanceChange* sin(wheelAngle * pi /180);
            rimAngle += distanceChange * 360 / (2 * pi * wheelRadius);
            break;
        case 's' :
            wheelCenter.x -= distanceChange* cos(wheelAngle * pi / 180);
            wheelCenter.y -= distanceChange* sin(wheelAngle * pi /180);
            rimAngle -= distanceChange * 360 / (2 * pi * wheelRadius);
            break;
        default:
            break;
    }
}


void specialKeyListener(int key, int x,int y){
    switch(key){
        case GLUT_KEY_DOWN:		//down arrow key
            cameraHeight -= 3.0;
            break;
        case GLUT_KEY_UP:		// up arrow key
            cameraHeight += 3.0;
            break;

        case GLUT_KEY_RIGHT:
            cameraAngle += 0.03;
            break;
        case GLUT_KEY_LEFT:
            cameraAngle -= 0.03;
            break;

        case GLUT_KEY_PAGE_UP:
            break;
        case GLUT_KEY_PAGE_DOWN:
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
    switch(button){
        case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
                drawaxes=1-drawaxes;
            }
            break;

        case GLUT_RIGHT_BUTTON:
            //........
            break;

        case GLUT_MIDDLE_BUTTON:
            //........
            break;

        default:
            break;
    }
}



void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
//    gluLookAt(0,0,200,	0,0,0,	0,1,0);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();
    drawWheel();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate(){
    angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init(){
    //codes for initialization
    drawgrid=1;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;
    wheelRadius = 25;
    wheelHeight = 8;
    wheelAngle = 0;
    rimAngle = 0;
    angleChange = 5;
    distanceChange = 10;

    initWheelCenter();
    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv){
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
