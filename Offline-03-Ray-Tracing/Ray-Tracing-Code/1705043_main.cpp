#include<windows.h>
#include<GL/glut.h>
#include "bitmap_image.hpp"
#include "1705043_classes.h"

int windowWidth = 500;
int windowHeight = 500;
double viewAngle = 80.0;
int drawaxes;
int pixelsAlongBothDimensions = 0;
int noOfObjects = 0;
int noOfPointLights = 0;
int bitmapImageCount;

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

Vector crossProduct(Vector a, Vector b){
    Vector result;
    result.x = a.y*b.z - b.y*a.z;
    result.y = a.z*b.x - b.z*a.x;
    result.z = a.x*b.y - b.x*a.y;
    return result;
}

Vector normalize(Vector v){
    double len = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    v.x /= len;
    v.y /= len;
    v.z /=len;
    return v;
}

Vector rotateVector(Vector v, Vector reference, double rotationAngle){
    Vector result;
    Vector vperp;
    vperp = crossProduct(reference, v);
    result.x = v.x * cos(rotationAngle*PI/180) + vperp.x * sin(rotationAngle*PI/180);
    result.y = v.y * cos(rotationAngle*PI/180) + vperp.y * sin(rotationAngle*PI/180);
    result.z = v.z * cos(rotationAngle*PI/180) + vperp.z * sin(rotationAngle*PI/180);
    result = normalize(result);
    return result;
}

void drawAxes()
{
    if(drawaxes == 1)
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

void capture() {
    cout << pos << ": Bitmap Image Capturing" << endl;
    int imageWidth = pixelsAlongBothDimensions;
    int imageHeight = pixelsAlongBothDimensions;

    //initialize bitmap image and set background color
    bitmap_image bitmapImage(imageWidth, imageHeight);

    for(int column=0; column<imageWidth; column++) {
        for(int row=0; row<imageHeight; row++) {
            bitmapImage.set_pixel(column, row, 0, 0, 0);  // color = black
        }
    }

    double planeDistance = windowHeight/(2.0*tan(viewAngle / 2.0 * PI / 180.0));
    Vector topLeft = pos + l * planeDistance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);

    double du = ((double) windowWidth/imageWidth);
    double dv = ((double) windowHeight/imageHeight);

    // Choose middle of the grid cell
    topLeft = topLeft + r * (du / 2.0) - u * (dv / 2.0);

    for(int column=0; column<imageWidth; column++) {
        for(int row=0; row<imageHeight; row++) {
            //calculate curPixel using topleft,r,u,i,j,du,dv
            Vector curPixel = topLeft + r * (column * du) - u * (row * dv);

            //cast ray from eye to (curPixel-eye) direction
            Ray ray(pos, curPixel-pos);

            int nearest = INF;
            double t, tMin=INF;

            for(int i=0; i<objects.size(); i++) {
                Color color;  // color = black
                t = objects[i]->intersect(ray, color, 0);

                if(t>0.0 && t<tMin) {
                    tMin = t;
                    nearest = i;
                }
            }

            /* finding color for current pixel */
            if(nearest != INF) {
                Color color;  // color = black
                tMin = objects[nearest]->intersect(ray, color, 1);
                bitmapImage.set_pixel(column, row, (int) round(color.red*255.0), (int) round(color.green*255.0), (int) round(color.blue*255.0));
            }
        }
    }

    //save image
    //The 1st image you capture after running
    //the program should be named Output_11.bmp, the 2nd image
    //you capture should be named Output_12.bmp and so on.
    bitmapImageCount++;
    string imageName = "Output_1"+to_string(bitmapImageCount)+".bmp";
    bitmapImage.save_image("../" + imageName);
    cout << "Bitmap Image Captured | Image Name - " + imageName << endl;
}

void moveBackward(){
    pos.x -= l.x;
    pos.y -= l.y;
    pos.z -= l.z;
}
void moveForward(){
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

void keyboardListener(unsigned char key, int x, int y) {
	switch(key) {
        case '0':
            capture();
            break;
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

void specialKeyListener(int key, int x, int y) {
	switch(key) {
        case GLUT_KEY_DOWN:		//down arrow key
            moveBackward();
            break;
        case GLUT_KEY_UP:		// up arrow key
            moveForward();
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

void drawObjects(){
    for(int i=0; i<objects.size(); i++) {
        objects[i]->draw();
    }
}

void drawLights(){
    for(int i=0; i<lights.size(); i++) {
        lights[i].draw();
    }
}

void display() {
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
	/* adding axes */
	drawAxes();

	/* adding objects */
    drawObjects();

	/* adding lights */
    drawLights();

	/* ADD this line in the end: if you use double buffer (i.e. GL_DOUBLE) */
	glutSwapBuffers();
}

void animate() {
	/* codes for any changes in Models, Camera */
	glutPostRedisplay();
}

void init() {
    init_pos_u_r_l();
    drawaxes = 1;

	bitmapImageCount = 0;

	/* clearing the screen */
	glClearColor(0, 0, 0, 0);  // color = black

	/* setting up projection here */
	/* loading the PROJECTION matrix */
	glMatrixMode(GL_PROJECTION);

	/* initializing the matrix */
	glLoadIdentity();

    /* setting the camera perspective by providing necessary parameters */
	gluPerspective(viewAngle, 1.0, 1.0, 1000.0);
}

void loadData() {
    ifstream fin;
    fin.open("../scene.txt");
    if(!fin.is_open()) {
        cout << "Error Opening File" << endl;
    }
    fin >> levelOfRecursion >> pixelsAlongBothDimensions;
    fin >> noOfObjects;

    string type;
    Object* temp;

    for(int i=0; i < noOfObjects; i++) {
        fin >> type;
        if(type == "sphere") {
            Vector center;
            double radius;

            fin >> center.x >> center.y >> center.z;
            fin >> radius;

            temp = new Sphere(center, radius);
        }
        else if(type == "triangle") {
            Vector a, b, c;

            fin >> a.x >> a.y >> a.z;
            fin >> b.x >> b.y >> b.z;
            fin >> c.x >> c.y >> c.z;

            temp = new Triangle(a, b, c);
        }
        else if(type == "general") {
            GeneralQuadricSurfaceCoefficient coefficient;
            Vector cubeReferencePoint;
            double length, width, height;

            fin >> coefficient;
            fin >> cubeReferencePoint.x >> cubeReferencePoint.y >> cubeReferencePoint.z;
            fin >> length >> width >> height;

            temp = new GeneralQuadricSurface(coefficient, cubeReferencePoint, length, width, height);
        }
        else {
            cout << type << " Shape Type error" << endl;
        }

        Color color;
        ReflectionCoefficient reflectionCoefficient;
        int shine;

        fin >> color.red >> color.green >> color.blue;
        fin >> reflectionCoefficient.ambient >> reflectionCoefficient.diffuse >> reflectionCoefficient.specular >> reflectionCoefficient.recursive;
        fin >> shine;

        temp->setColor(color);
        temp->setReflectionCoefficient(reflectionCoefficient);
        temp->setShine(shine);

        objects.push_back(temp);
    }
    temp = NULL;

    fin >> noOfPointLights;

    for(int i=0; i < noOfPointLights; i++) {
        Vector position;
        Color color;

        fin >> position.x >> position.y >> position.z;
        fin >> color.red >> color.green >> color.blue;

        lights.push_back(PointLight(position, color, 1.0));
    }
    fin.close();

    /* creating a floor temp and pushing it to objects vector */
    temp = new Floor(1000.0, 20.0, Color());  // color = black

    temp->setColor(Color(1.0, 1.0, 1.0));  // color = white
    temp->setReflectionCoefficient(ReflectionCoefficient(0.25, 0.25, 0.25, 0.25));
    temp->setShine(15);

    objects.push_back(temp);
    temp = NULL;
}

int main(int argc, char **argv){
    glutInit(&argc,argv);
    glutInitWindowSize(windowWidth, windowHeight);
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

    loadData();

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}