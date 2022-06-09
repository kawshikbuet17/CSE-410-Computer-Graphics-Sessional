#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/freeglut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double length;
double posXangle;
double negXangle;
double pos2Xangle;
double posYangle;
double negYangle;
double angleChange;

struct point
{
    double x,y,z;
};


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
            for(i=-8;i<=8;i++){

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }glEnd();
    }
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);{
        glVertex3f( a, a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }glEnd();
}


void keyboardListener(unsigned char key, int x,int y){
    switch(key){

        case '1':
            posYangle += angleChange;
            if (posYangle > 90){
                posYangle = 90;
            }
            break;
        case '2':
            posYangle -= angleChange;
            if (posYangle < 0){
                posYangle = 0;
            }
            break;
        case '3':
            negYangle -= angleChange;
            if (negYangle < -90){
                negYangle = -90;
            }
            break;
        case '4':
            negYangle += angleChange;
            if (negYangle >= 0){
                negYangle = 0;
            }
            break;
        case '5':
            negXangle += angleChange;
            if (negXangle >= 90){
                negXangle = 90;
            }
            break;
        case '6':
            negXangle -= angleChange;
            if (negXangle <= 0){
                negXangle = 0;
            }
            break;
        case '7':
            posXangle -= angleChange;
            if (posXangle < -90){
                posXangle = -90;
            }
            break;
        case '8':
            posXangle += angleChange;
            if (posXangle >= 0){
                posXangle = 0;
            }
            break;
        case '9':
            pos2Xangle -= angleChange;
            if (pos2Xangle < -90){
                pos2Xangle = -90;
            }
            break;
        case '0':
            pos2Xangle += angleChange;
            if (pos2Xangle >= 0){
                pos2Xangle = 0;
            }
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

void drawCenterSquare(){
    glPushMatrix();
    glColor3f(1,0,0);
    drawSquare(length);
    glPopMatrix();
}

void drawPosXSquare(){
    glPushMatrix();
    glColor3f(1,1,0);
    glTranslatef(length, 0, 0);
    glRotatef(posXangle, 0, 1, 0);
    glTranslatef(length, 0, 0);
    drawSquare(length);
}

void drawNegXSquare(){
    glPushMatrix();
    glColor3f(0,0,1);
    glTranslatef(-length, 0, 0);
    glRotatef(negXangle, 0, 1, 0);
    glTranslatef(-length, 0, 0);
    drawSquare(length);
    glPopMatrix();
}

void drawPos2XSquare(){
    glColor3f(1,0,0);
    glTranslatef(length, 0, 0);
    glRotatef(pos2Xangle, 0, 1, 0);
    glTranslatef(length, 0, 0);
    drawSquare(length);
    glPopMatrix();
}

void drawPosYSquare(){
    glPushMatrix();
    glColor3f(0,1,0);
    glTranslatef(0, length, 0);
    glRotatef(posYangle, 1, 0, 0);
    glTranslatef(0, length, 0);
    drawSquare(length);
    glPopMatrix();
}

void drawNegYSquare(){
    glPushMatrix();
    glColor3f(0,1,1);
    glTranslatef(0, -length, 0);
    glRotatef(negYangle, 1, 0, 0);
    glTranslatef(0, -length, 0);
    drawSquare(length);
    glPopMatrix();
}



void draw6Squares(){
    drawNegYSquare();
    drawPosYSquare();
    drawNegXSquare();
    drawPosXSquare();
    drawPos2XSquare();
    drawCenterSquare();
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

    gluLookAt(100,100,100,	0,0,0,	0,0,2);
//    gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
//    gluLookAt(0,0,200,	0,0,0,	0,1,0);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    //glColor3f(1,0,0);
    draw6Squares();

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
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;
    length = 20;
    posXangle = 0;
    negXangle = 0;
    pos2Xangle = 0;
    posYangle = 0;
    negYangle = 0;
    angleChange = 5;
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
