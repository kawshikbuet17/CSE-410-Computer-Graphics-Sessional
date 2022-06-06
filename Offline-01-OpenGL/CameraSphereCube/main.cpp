#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/freeglut.h>

#include<bits/stdc++.h>
using namespace std;

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double maxLen = 30;
double currentLen = 20;

struct point{
    double x,y,z;
};

struct point pos, u, r, l;

void init_pos_u_r_l(){
    pos.x = 100,
    pos.y = 100,
    pos.z = 0;

    u.x = 0,
    u.y = 0,
    u.z = 1;

    r.x = -1/ sqrt(2),
    r.y = 1/ sqrt(2),
    r.z = 0;

    l.x = -1/ sqrt(2),
    l.y = -1/ sqrt(2),
    l.z = 0;
}

struct point crossProduct(struct point a, struct point b){
    struct point result{};
    result.x = a.y*b.z - b.y*a.z;
    result.y = a.z*b.x - b.z*a.x;
    result.z = a.x*b.y - b.x*a.y;
    return result;
}

struct point normalize(struct point v)
{
    double len = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    v.x /= len;
    v.y /= len;
    v.z /=len;
    return v;
}

struct point rotateVector(struct point v, struct point reference, double rotationAngle){
    struct point result{};
    struct point vperp{};

    vperp = crossProduct(reference, v);

    result.x = v.x * cos(rotationAngle*pi/180) + vperp.x * sin(rotationAngle*pi/180);
    result.y = v.y * cos(rotationAngle*pi/180) + vperp.y * sin(rotationAngle*pi/180);
    result.z = v.z * cos(rotationAngle*pi/180) + vperp.z * sin(rotationAngle*pi/180);
    result = normalize(result);
    return result;
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


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks, string hemisphere, double portion = 1)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0;i<=stacks;i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0;j<=slices;j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi * portion);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi * portion);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        glColor3f(1, 0, 0);   //red color
        for(j=0;j<slices;j++)
        {
            glBegin(GL_QUADS);{
                //upper hemisphere
                if(hemisphere == "upper"){
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                }
                //lower hemisphere
                else if(hemisphere == "lower"){
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
            }glEnd();
        }
    }
}

void draw8Sphere(){
    float trans = currentLen;
    float radius = maxLen - currentLen;
    glColor3f(1, 0, 0);   //red color

    glPushMatrix();
    glTranslatef(trans, trans, trans);
    drawSphere(radius, 20, 20, "upper", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, trans, trans);
    glRotatef(90, 0, 0, 1);
    drawSphere(radius, 20, 20, "upper", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, -trans, trans);
    glRotatef(180, 0, 0, 1);
    drawSphere(radius, 20, 20, "upper", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans, -trans, trans);
    glRotatef(270, 0, 0, 1);
    drawSphere(radius, 20, 20, "upper", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans, trans, -trans);
    drawSphere(radius, 20, 20, "lower", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, trans, -trans);
    glRotatef(90, 0, 0, 1);
    drawSphere(radius, 20, 20, "lower", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, -trans, -trans);
    glRotatef(180, 0, 0, 1);
    drawSphere(radius, 20, 20, "lower", 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans, -trans, -trans);
    glRotatef(270, 0, 0, 1);
    drawSphere(radius, 20, 20, "lower", 0.25);
    glPopMatrix();
}

void drawCylinder(double radius, double height, int slices,int stacks, double portion = 1)
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
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi * portion);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi * portion);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        glColor3f(0, 1, 0);   //red color
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

void draw12Cylinder(){
    float trans = currentLen;
    float radius = maxLen - currentLen;
    glColor3f(1, 0, 0);   //red color

    glPushMatrix();
    glTranslatef(trans, trans, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, trans, 0);
    glRotatef(90, 0, 0, 1);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, -trans, 0);
    glRotatef(180, 0, 0, 1);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans, -trans, 0);
    glRotatef(270, 0, 0, 1);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();



    glPushMatrix();
    glTranslatef(trans, 0, trans);
    glRotatef(90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, trans, trans);
    glRotatef(90, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, 0, trans);
    glRotatef(180, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, -trans, trans);
    glRotatef(270, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();



    glPushMatrix();
    glTranslatef(trans, 0, -trans);
    glRotatef(-90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, trans, -trans);
    glRotatef(90, 0, 0, 1);
    glRotatef(-90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-trans, 0, -trans);
    glRotatef(180, 0, 0, 1);
    glRotatef(-90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, -trans, -trans);
    glRotatef(270, 0, 0, 1);
    glRotatef(-90, 1, 0, 0);
    drawCylinder(radius, currentLen, 20, 20, 0.25);
    glPopMatrix();
}

void drawSquareOfCube(string axis){
    float tx = 0, ty = 0, tz = 0;
    float rx = 0, ry = 0, rz = 0;
    float rangle = 0;

    if(axis == "+z"){
        tz = maxLen;
    }
    else if(axis == "-z"){
        tz = -maxLen;
    }
    else if(axis == "+x"){
        tx = maxLen;
        ry = 1;
        rangle = 90;
    }
    else if(axis == "-x"){
        tx = -maxLen;
        ry = 1;
        rangle = 90;
    }
    else if(axis == "+y"){
        ty = maxLen;
        rx = 1;
        rangle = 90;
    }
    else if(axis == "-y"){
        ty = -maxLen;
        rx = 1;
        rangle = 90;
    }

    glPushMatrix();
    glTranslatef(tx,ty,tz);
    glRotatef(rangle, rx, ry, rz);
    drawSquare(currentLen);
    glPopMatrix();
}

void drawCube(){
    glColor3f(1.0,1.0,1.0);
    drawSquareOfCube("+z");
    drawSquareOfCube("-z");
    drawSquareOfCube("+x");
    drawSquareOfCube("-x");
    drawSquareOfCube("+y");
    drawSquareOfCube("-y");
}

void moveForward(){
    pos.x -= l.x;
    pos.y -= l.y;
    pos.z -= l.z;
}
void moveBackward(){
    pos.x += l.x;
    pos.y += l.y;
    pos.z += l.z;
}
void moveLeft(){
    pos.x -= r.x;
    pos.y -= r.y;
    pos.z -= r.z;
}
void moveRight(){
    pos.x += r.x;
    pos.y += r.y;
    pos.z += r.z;
}
void moveUp(){
    pos.x += u.x;
    pos.y += u.y;
    pos.z += u.z;
}
void moveDown(){
    pos.x -= u.x;
    pos.y -= u.y;
    pos.z -= u.z;
}

void lookLeft(){
    l = rotateVector(l, u, 2);
    r = rotateVector(r, u, 2);
}

void lookRight(){
    l = rotateVector(l, u, -2);
    r = rotateVector(r, u, -2);
}

void lookUp(){
    u = rotateVector(u, r, 2);
    l = rotateVector(l, r, 2);
}

void lookDown(){
    u = rotateVector(u, r, -2);
    l = rotateVector(l, r, -2);
}

void tiltClockwise(){
    u = rotateVector(u, l, -2);
    r = rotateVector(r, l, -2);
}

void tiltAntiClockwise(){
    u = rotateVector(u, l, 2);
    r = rotateVector(r, l, 2);
}

void keyboardListener(unsigned char key, int x,int y){
    switch(key){

        case '1':
            lookLeft();
            break;
        case '2':
            lookRight();
            break;
        case '3':
            lookUp();
            break;
        case '4':
            lookDown();
            break;
        case '5':
            tiltClockwise();
            break;
        case '6':
            tiltAntiClockwise();
            break;
        default:
            break;
    }
}




void specialKeyListener(int key, int x,int y){
    switch(key){
        case GLUT_KEY_DOWN:		//down arrow key
            moveForward();
            break;
        case GLUT_KEY_UP:		// up arrow key
            moveBackward();
            break;
        case GLUT_KEY_RIGHT:
            moveRight();
            break;
        case GLUT_KEY_LEFT:
            moveLeft();
            break;

        case GLUT_KEY_PAGE_UP:
            moveUp();
            break;
        case GLUT_KEY_PAGE_DOWN:
            moveDown();
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            if(currentLen > 0){
                currentLen -= 0.20;
            }
            if(currentLen < 0){
                currentLen = 0;
            }
            break;
        case GLUT_KEY_END:
            if(currentLen <= maxLen){
                currentLen += 0.20;
            }
            if(currentLen > maxLen){
                currentLen = maxLen;
            }
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
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    gluLookAt(pos.x,pos.y,pos.z,	pos.x + l.x,pos.y + l.y,pos.z + l.z,	u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    draw8Sphere();
    drawCube();
    draw12Cylinder();


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
    init_pos_u_r_l();

    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;

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
