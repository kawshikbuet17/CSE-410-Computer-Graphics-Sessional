#include<bits/stdc++.h>
using namespace std;

#define pi 2*acos(0.0)

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;

vector<vector<double>> matrix;
stack<pair<vector<vector<double>>, bool>> matrixStack;

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

vector<vector<double>> triangle(ifstream &fin, ofstream &fout, int stage = 1, Point *a = NULL, Point *b = NULL, Point *c = NULL){
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
            matrix = triangle(fin, fout);
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
        vector<vector<double>> stage1Tringle =  triangle(fin, fout, 2, p1, p2, p3);
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
        vector<vector<double>> stage2Tringle =  triangle(fin, fout, 3, p1, p2, p3);
        vector<vector<double>> projectionTransformMatrix = matrixMultiplication(P, stage2Tringle);
        scaleHomogeneous(projectionTransformMatrix);
        outputToFile(fout, transposeMatrix(projectionTransformMatrix));
    }
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
}
