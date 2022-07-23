#include<bits/stdc++.h>
using namespace std;
#include "bitmap_image.hpp"

#define pi 2*acos(0.0)
#define INF INT_MAX

class Point{
    public:
        double x, y, z;
        Point(){
            x = y = z = 0.0;
        }
        Point(double x, double y, double z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
        void normalize(){
            double len = sqrt(x*x + y*y + z*z);
            x /= len;
            y /= len;
            z /= len;
        }
};

class Color{
    public:
        double r, g, b;
        Color(){
            r = g = b = 0.0;
        }
        Color(double r, double g, double b){
            this->r = r;
            this->g = g;
            this->b = b;
        }
};

class Triangle{
    public:
        Point points[3];
        Color colors;
        Triangle(Point p1, Point p2, Point p3, Color rgb){
            points[0] = p1;
            points[1] = p2;
            points[2] = p3;

            colors = rgb;
        }
};

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;

vector<vector<double>> matrix;
stack<pair<vector<vector<double>>, bool>> matrixStack;

int Screen_Width, Screen_Height;
double leftLimitOfX, rightLimitOfX, bottomLimitOfY, topLimitOfY, frontLimitOfZ, rearLimitOfZ;

double dx, dy;
double Top_Y, Bottom_Y, Left_X, Right_X;
double **zBuffer;
Color **frameBuffer;


vector<vector<double>> identityMatrix(){
    vector<vector<double>> matrix(4, vector<double>(4, 0));
    for(int i = 0; i < 4; i++){
        matrix[i][i] = 1;
    }
    return matrix;
}

vector<vector<double>> transposeMatrix(vector<vector<double>> matrix){
    vector<vector<double>> transposed(4, vector<double>(4, 0));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}

vector<vector<double>> matrixMultiplication(vector<vector<double>> A, vector<vector<double>> B){
    vector<vector<double>> C(4, vector<double>(4, 0));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            double sum = 0;
            for(int k = 0; k < 4; k++){
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

vector<vector<double>> inputTriangle(ifstream &fin, ofstream &fout, int stage = 1, Point *a = NULL, Point *b = NULL, Point *c = NULL){
    Point p1, p2, p3;
    if(stage == 1){
        fin >> p1.x >> p1.y >> p1.z;
        fin >> p2.x >> p2.y >> p2.z;
        fin >> p3.x >> p3.y >> p3.z;
    }
    else if(stage == 2 or stage == 3){
        p1 = *a;
        p2 = *b;
        p3 = *c;
    }

    //making the points as homogeneous coordinates
    vector<vector<double>> matrix(4, vector<double>(4, 0));
    matrix[0][0] = p1.x;
    matrix[1][0] = p1.y;
    matrix[2][0] = p1.z;
    matrix[3][0] = 1;

    matrix[0][1] = p2.x;
    matrix[1][1] = p2.y;
    matrix[2][1] = p2.z;
    matrix[3][1] = 1;

    matrix[0][2] = p3.x;
    matrix[1][2] = p3.y;
    matrix[2][2] = p3.z;
    matrix[3][2] = 1;

    matrix[0][3] = 1;
    matrix[1][3] = 1;
    matrix[2][3] = 1;
    matrix[3][3] = 1;

    return matrix;
}

vector<vector<double>> translate(ifstream &fin, ofstream &fout){
    double tx, ty, tz;
    fin >> tx >> ty >> tz;

    //make translation matrix vector
    vector<vector<double>> matrix(4, vector<double>(4, 0));
    matrix[0][0] = 1;
    matrix[1][1] = 1;
    matrix[2][2] = 1;
    matrix[3][3] = 1;
    matrix[0][3] = tx;
    matrix[1][3] = ty;
    matrix[2][3] = tz;
    
    return matrix;
}

vector<vector<double>> scale(ifstream &fin, ofstream &fout){
    double sx, sy, sz;
    fin >> sx >> sy >> sz;

    //make scale matrix vector
    vector<vector<double>> matrix(4, vector<double>(4, 0));
    matrix[0][0] = sx;
    matrix[1][1] = sy;
    matrix[2][2] = sz;
    matrix[3][3] = 1;
    return matrix;
}

Point crossProduct(Point p1, Point p2){
    Point p;
    p.x = p1.y * p2.z - p1.z * p2.y;
    p.y = p1.z * p2.x - p1.x * p2.z;
    p.z = p1.x * p2.y - p1.y * p2.x;
    return p;
}

double dotProduct(Point p1, Point p2){
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Point Rodrigues(Point a, Point x, double theta){
    double cosTHETA = cos(theta * pi / 180);
    double sinTHETA = sin(theta * pi / 180);
    
    Point A = Point(x.x * cosTHETA, x.y * cosTHETA, x.z * cosTHETA);

    Point B = Point(a.x * (1 - cosTHETA) * dotProduct(a, x), a.y * (1 - cosTHETA) * dotProduct(a, x), a.z * (1 - cosTHETA) * dotProduct(a, x));

    Point C = crossProduct(a, x);
    C = Point(C.x * sinTHETA, C.y * sinTHETA, C.z * sinTHETA);

    Point p(A.x + B.x + C.x, A.y + B.y + C.y, A.z + B.z + C.z);
    
    return p;
}

vector<vector<double>> rotate(ifstream &fin, ofstream &fout){
    double theta, ax, ay, az;
    fin >> theta >> ax >> ay >> az;

    double sqrtaxayaz = sqrt(ax * ax + ay * ay + az * az);

    Point a = Point((double)ax / sqrtaxayaz, (double)ay / sqrtaxayaz, (double)az / sqrtaxayaz);

    Point x = Point(1, 0, 0);
    Point y = Point(0, 1, 0);
    Point z = Point(0, 0, 1);

    Point c1 = Rodrigues(a, x, theta);
    Point c2 = Rodrigues(a, y, theta);
    Point c3 = Rodrigues(a, z, theta);

    vector<vector<double>> matrix(4, vector<double>(4, 0));
    matrix[0][0] = c1.x;
    matrix[1][0] = c1.y;
    matrix[2][0] = c1.z;

    matrix[0][1] = c2.x;
    matrix[1][1] = c2.y;
    matrix[2][1] = c2.z;

    matrix[0][2] = c3.x;
    matrix[1][2] = c3.y;
    matrix[2][2] = c3.z;

    matrix[0][3] = 0;
    matrix[1][3] = 0;
    matrix[2][3] = 0;
    matrix[3][3] = 1;

    return matrix;
}

void push(stack<pair<vector<vector<double>>, bool>>& matrixStack, vector<vector<double>> matrix, bool flag){
    matrixStack.push({matrix, flag});
}

void pop(stack<pair<vector<vector<double>>, bool>>& matrixStack, ofstream &fout){
    while(!matrixStack.empty()){
        if(matrixStack.top().second){
            matrixStack.top().second = false;
            break;
        }
        else{
            matrixStack.pop();
        }
    }
}

void outputToFile(ofstream &fout, vector<vector<double>> matrix){
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            fout <<fixed <<setprecision(7) << matrix[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

void fileIOtest(ifstream &fin, ofstream &fout){
    double a, b, c;
    while(fin>>a>>b>>c){
        fout << a << " " << b << " " << c << endl;
    }
}

//scale homogeneous matrix
void scaleHomogeneous(vector<vector<double>> &matrix){
    for(int i = 0; i< 4; i++){
        for(int j = 0; j< 4; j++){
            matrix[i][j] /= matrix[3][j];
        }
    }
}

void stage1(ifstream &fin, ofstream &fout){
    push(matrixStack, identityMatrix(), false);
    string command;
    while(fin >> command){
        if(command == "end"){
            break;
        }
        else if(command == "triangle"){
            matrix = inputTriangle(fin, fout);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                scaleHomogeneous(temp);
                outputToFile(fout, transposeMatrix(temp));
            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "translate"){
            matrix = translate(fin, fout);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                push(matrixStack, temp, false);
            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "scale"){
            matrix = scale(fin, fout);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                push(matrixStack, temp, false);

            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "rotate"){
            matrix = rotate(fin, fout);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                push(matrixStack, temp, false);

            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }

        else if(command == "push"){
            matrix = matrixStack.top().first;
            matrixStack.pop();
            push(matrixStack, matrix, true);
        }
        else if(command == "pop"){
            if(!matrixStack.empty()){
                pop(matrixStack, fout);
            }else{
                fout<<"Stack is Empty. Can't pop"<<endl;
            }    
        }
    }
}

void stage2(ifstream &fin, ofstream &fout){
    Point l = Point(lookX-eyeX, lookY-eyeY, lookZ-eyeZ);
    l.normalize();

    Point up = Point(upX, upY, upZ);

    Point r = crossProduct(l, up);
    r.normalize();

    Point u = crossProduct(r, l);

    vector<vector<double>> T = identityMatrix();
    T[0][3] = -eyeX;
    T[1][3] = -eyeY;
    T[2][3] = -eyeZ;
    vector<vector<double>> R = identityMatrix();
    R[0][0] = r.x;
    R[0][1] = r.y;
    R[0][2] = r.z;
    R[1][0] = u.x;
    R[1][1] = u.y;
    R[1][2] = u.z;
    R[2][0] = -l.x;
    R[2][1] = -l.y;
    R[2][2] = -l.z;
    
    //view transformation matrix
    vector<vector<double>> V = matrixMultiplication(T, R);

    Point *p1, *p2, *p3;
    p1 = new Point();
    p2 = new Point();
    p3 = new Point();
    while(fin >> p1->x >> p1->y >> p1->z >> p2->x >> p2->y >> p2->z >> p3->x >> p3->y >> p3->z){
        vector<vector<double>> stage1Tringle =  inputTriangle(fin, fout, 2, p1, p2, p3);
        vector<vector<double>> viewTransformMatrix = matrixMultiplication(V, stage1Tringle);
        scaleHomogeneous(viewTransformMatrix);
        outputToFile(fout, transposeMatrix(viewTransformMatrix));
    }
}

void stage3(ifstream &fin, ofstream &fout){
    double fovX = fovY * aspectRatio;
    double r = near * tan(fovX/2 * pi / 180);
    double t = near * tan(fovY/2 * pi / 180);

    //projection transformation matrix
    vector<vector<double>> P = identityMatrix();
    P[0][0] = near/r;
    P[1][1] = near/t;
    P[2][2] = -(far + near)/(far - near);
    P[2][3] = -2*far*near/(far - near);
    P[3][2] = -1;
    P[3][3] = 0;

    Point *p1, *p2, *p3;
    p1 = new Point();
    p2 = new Point();
    p3 = new Point();
    while(fin >> p1->x >> p1->y >> p1->z >> p2->x >> p2->y >> p2->z >> p3->x >> p3->y >> p3->z){
        vector<vector<double>> stage2Tringle =  inputTriangle(fin, fout, 3, p1, p2, p3);
        vector<vector<double>> projectionTransformMatrix = matrixMultiplication(P, stage2Tringle);
        scaleHomogeneous(projectionTransformMatrix);
        outputToFile(fout, transposeMatrix(projectionTransformMatrix));
    }
}

double **initializeZBuffer(){
    double **zBuffer = new double*[Screen_Height];
    for(int i = 0; i < Screen_Height; i++){
        zBuffer[i] = new double[Screen_Width];
    }
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            zBuffer[i][j] = rearLimitOfZ;
        }
    }
    return zBuffer;
}


Color **initializeFrameBuffer(){
    Color **frameBuffer = new Color*[Screen_Height];
    for(int i = 0; i < Screen_Height; i++){
        frameBuffer[i] = new Color[Screen_Width];
    }
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            frameBuffer[i][j] = Color(0, 0, 0);
        }
    }
    return frameBuffer;
}



void readData(ifstream &fin){
    fin >> Screen_Height >> Screen_Width;
    fin >> leftLimitOfX;
    fin >> bottomLimitOfY;
    fin >> frontLimitOfZ >> rearLimitOfZ;

    rightLimitOfX = -leftLimitOfX;
    topLimitOfY = -bottomLimitOfY;

    fin.close();
}

int findTopScanLine(Triangle *triangle){
    double max_y = max(triangle->points[0].y, max(triangle->points[1].y, triangle->points[2].y));
    
    int top_scanline;

    if(max_y >= Top_Y){
        top_scanline = 0;
    }
    else{
        top_scanline = round((Top_Y - max_y) / dy);
    }
    return top_scanline;
}

int findBottomScanLine(Triangle *triangle){
    double min_y = min(triangle->points[0].y, min(triangle->points[1].y, triangle->points[2].y));
    
    int bottom_scanline;

    if(min_y <= Bottom_Y){
        bottom_scanline = Screen_Height - 1;
    }
    else{
        bottom_scanline = Screen_Height - round((min_y - Bottom_Y) / dy) - 1;
    }
    return bottom_scanline;
}

int findLeftIntersectingColumn(double xa){
    int left_intersecting_column = round((xa - Left_X) / dx);
    left_intersecting_column = max(left_intersecting_column, 0);
    return left_intersecting_column;
}

int findRightIntersectingColumn(double xb){
    int right_intersecting_column = round((xb - Left_X) / dx);
    right_intersecting_column = min(right_intersecting_column, Screen_Width - 1);
    return right_intersecting_column;
}

vector<double> findX(Triangle *triangle, double y){
    //ABC triangle
    //A = triangle->points[0]
    //B = triangle->points[1]
    //C = triangle->points[2]

    //intersect with AB
    // (x - xA)/(xA-xB) = (y - yA)/(yA-yB)
    double xAB = min((y - triangle->points[0].y) * (triangle->points[0].x - triangle->points[1].x) / (triangle->points[0].y - triangle->points[1].y) + triangle->points[0].x, (double)INF);
    //intersect with AC
    // (x - xA)/(xA-xC) = (y - yA)/(yA-yC)
    double xAC = min((y - triangle->points[0].y) * (triangle->points[0].x - triangle->points[2].x) / (triangle->points[0].y - triangle->points[2].y) + triangle->points[0].x, (double)INF);

    //intersect with BC
    // (x - xB)/(xB-xC) = (y - yB)/(yB-yC)
    double xBC = min((y - triangle->points[1].y) * (triangle->points[1].x - triangle->points[2].x) / (triangle->points[1].y - triangle->points[2].y) + triangle->points[1].x, (double)INF);
    vector<double> x;
    x.push_back(xAB);
    x.push_back(xAC);
    x.push_back(xBC);
    return x;
}

vector<double> findZ(Triangle *triangle, double y){
    //ABC triangle
    //A = triangle->points[0]
    //B = triangle->points[1]
    //C = triangle->points[2.
    //intersect with AB
    // (z - zA)/(zA-zB) = (y - yA)/(yA-yB)
    double zAB = min((y - triangle->points[0].y) * (triangle->points[0].z - triangle->points[1].z) / (triangle->points[0].y - triangle->points[1].y) + triangle->points[0].z, (double)INF);
    
    //intersect with AC
    // (z - zA)/(zA-zC) = (y - yA)/(yA-yC)
    double zAC = min((y - triangle->points[0].y) * (triangle->points[0].z - triangle->points[2].z) / (triangle->points[0].y - triangle->points[2].y) + triangle->points[0].z, (double)INF);

    //intersect with BC
    // (z - zB)/(zB-zC) = (y - yB)/(yB-yC)
    double zBC = min((y - triangle->points[1].y) * (triangle->points[1].z - triangle->points[2].z) / (triangle->points[1].y - triangle->points[2].y) + triangle->points[1].z, (double)INF);

    vector<double> z;
    z.push_back(zAB);
    z.push_back(zAC);
    z.push_back(zBC);
    return z;
}

void applyProcedure(ifstream &fin, ofstream &fout){
    fin.open("stage3.txt");
    Point *p1, *p2, *p3;
    p1 = new Point();
    p2 = new Point();
    p3 = new Point();

    while(fin >> p1->x >> p1->y >> p1->z >> p2->x >> p2->y >> p2->z >> p3->x >> p3->y >> p3->z){
        Color rgb = Color(rand()%256, rand()%256, rand()%256);
        Triangle *triangle = new Triangle(*p1, *p2, *p3, rgb);

        int top_scanline = findTopScanLine(triangle);
        int bottom_scanline = findBottomScanLine(triangle);

        for(int row_no = top_scanline; row_no <= bottom_scanline; row_no++){
            double y = Top_Y - row_no * dy;
            vector<double> x = findX(triangle, y);
            vector<double> z = findZ(triangle, y);

            Point A = triangle->points[0];
            Point B = triangle->points[1];
            Point C = triangle->points[2];

            //if intersection point lies between the A and B point
            bool keep_intersection_0 = (x[0] >= A.x && x[0] <= B.x) || (x[0] >= B.x && x[0] <= A.x);

            //if intersection point lies between the C and A point
            bool keep_intersection_1 = (x[1] >= A.x && x[1] <= C.x) || (x[1] >= C.x && x[1] <= A.x);

            //if intersection point lies between the B and C point
            bool keep_intersection_2 = (x[2] >= B.x && x[2] <= C.x) || (x[2] >= C.x && x[2] <= B.x);
            

            double xa = 0, xb = 0, za = 0, zb = 0;
            if(keep_intersection_0 && keep_intersection_1){
                if(x[0] < x[1]){
                    xa = x[0];
                    xb = x[1];
                    za = z[0];
                    zb = z[1];
                }
                else{
                    xa = x[1];
                    xb = x[0];
                    za = z[1];
                    zb = z[0];
                }
            }
            else if(keep_intersection_1 && keep_intersection_2){
                if(x[1] < x[2]){
                    xa = x[1];
                    xb = x[2];
                    za = z[1];
                    zb = z[2];
                }
                else{
                    xa = x[2];
                    xb = x[1];
                    za = z[2];
                    zb = z[1];
                }
            }
            else if(keep_intersection_2 && keep_intersection_0){
                if(x[2] < x[0]){
                    xa = x[2];
                    xb = x[0];
                    za = z[2];
                    zb = z[0];
                }
                else{
                    xa = x[0];
                    xb = x[2];
                    za = z[0];
                    zb = z[2];
                }
            }
            
            int left_intersecting_column = findLeftIntersectingColumn(xa);
            int right_intersecting_column = findRightIntersectingColumn(xb);

            double zp;
            double constantTerm = (zb - za) / (xb - xa);
            for(int col_no = left_intersecting_column; col_no <= right_intersecting_column; col_no++){
                if(col_no == left_intersecting_column){
                    double temp = (Left_X + col_no * dx - xa) * (zb - za) / (xb - xa); 
                    zp = za + temp;
                }else{
                    zp = zp + dx * constantTerm;
                }

                //Compare with Z-buffer and z_front_limit and update if required
                if(zp > frontLimitOfZ && zp < zBuffer[row_no][col_no]){
                    zBuffer[row_no][col_no] = zp;
                    frameBuffer[row_no][col_no] = rgb;
                }
            }
        }
    }
    fin.close();

    //saving outputs
    bitmap_image image(Screen_Width, Screen_Height);

    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            image.set_pixel(j, i, frameBuffer[i][j].r, frameBuffer[i][j].g, frameBuffer[i][j].b);
        }
    }
    image.save_image("output.bmp");

    //output z buffer value to a file
    fout.open("z_buffer.txt");
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            if(zBuffer[i][j] < rearLimitOfZ)
                fout << fixed <<setprecision(6) << zBuffer[i][j] << "\t";
        }
        fout << endl;
    }

    //freeing memory
    for(int i = 0; i < Screen_Height; i++){
        delete[] frameBuffer[i];
    }
    delete[] frameBuffer;

    for(int i = 0; i < Screen_Height; i++){
        delete[] zBuffer[i];
    }
    delete[] zBuffer;
    fout.close();
}



void stage4(ifstream &fin, ofstream &fout){
    fin.open("config.txt");
    readData(fin);
    fin.close();
    dx = (rightLimitOfX - leftLimitOfX) / Screen_Width;
    dy = (topLimitOfY - bottomLimitOfY) / Screen_Height;

    Top_Y = topLimitOfY - dy/2;
    Bottom_Y = -topLimitOfY + dy/2;

    Left_X = leftLimitOfX + dx/2;
    Right_X = leftLimitOfX - dx/2;

    zBuffer = initializeZBuffer();
    frameBuffer = initializeFrameBuffer();

    applyProcedure(fin, fout);
}


int main(){
    ifstream fin;
    ofstream fout;

    fin.open("scene.txt");
    fout.open("stage1.txt");
    fin >> eyeX >> eyeY >> eyeZ;
    fin >> lookX >> lookY >> lookZ;
    fin >> upX >> upY >> upZ;
    fin >> fovY >> aspectRatio >> near >> far;
    stage1(fin, fout);
    fin.close();
    fout.close();

    fin.open("stage1.txt");
    fout.open("stage2.txt");
    stage2(fin, fout);
    fin.close();
    fout.close();

    fin.open("stage2.txt");
    fout.open("stage3.txt");
    stage3(fin, fout);
    fin.close();
    fout.close();

    stage4(fin, fout);

    return 0;
}
