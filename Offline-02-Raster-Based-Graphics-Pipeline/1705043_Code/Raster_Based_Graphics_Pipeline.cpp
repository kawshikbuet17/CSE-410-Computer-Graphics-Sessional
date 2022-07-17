#include<bits/stdc++.h>
using namespace std;

#define pi 2*acos(0.0)

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
};

vector<vector<double>> identityMatrix(){
    vector<vector<double>> matrix(4, vector<double>(4, 0));
    for(int i = 0; i < 4; i++){
        matrix[i][i] = 1;
    }
    return matrix;
}

// transpose of a matrix
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

vector<vector<double>> triangle(ifstream &fin, ofstream &fout){
    Point p1, p2, p3;
    fin >> p1.x >> p1.y >> p1.z;
    fin >> p2.x >> p2.y >> p2.z;
    fin >> p3.x >> p3.y >> p3.z;

    fout << "Input Triangle" << endl;
    fout << p1.x << " " << p1.y << " " << p1.z << endl;
    fout << p2.x << " " << p2.y << " " << p2.z << endl;
    fout << p3.x << " " << p3.y << " " << p3.z << endl;
    fout << endl;

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

    fout << "Input Translation" << endl;
    fout << tx << " " << ty << " " << tz << endl;
    fout << endl;

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

    fout << "Input Scale" << endl;
    fout << sx << " " << sy << " " << sz << endl;
    fout << endl;

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

    fout << "Input Rotation" << endl;
    fout << theta << " " << ax << " " << ay << " " << az << endl;
    fout << endl;

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
        fout << "Stack size = " << matrixStack.size() << " | Current Flag is " << matrixStack.top().second << endl;
        if(matrixStack.top().second){
            matrixStack.top().second = false;
            fout << "Popping Stopped" << endl;
            break;
        }
        else{
            fout << "Popping" << endl;
            matrixStack.pop();
        }
    }
    fout << endl;
}

void outputToFile(ofstream &fout, vector<vector<double>> matrix){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            fout << matrix[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main(){
    vector<vector<double>> matrix;
    stack<pair<vector<vector<double>>, bool>> matrixStack;
    push(matrixStack, identityMatrix(), false);

    ifstream fin;
    fin.open("scene.txt");

    ofstream fout;
    fout.open("stage1.txt");

    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;

    fin >> eyeX >> eyeY >> eyeZ;
    fin >> lookX >> lookY >> lookZ;
    fin >> upX >> upY >> upZ;
    fin >> fovY >> aspectRatio >> near >> far;

    string command;
    while(fin >> command){
        fout << "Command = " << command << endl;
        if(!matrixStack.empty()){
            fout<<"Top of stack | Flag " << matrixStack.top().second << " | Stack size = "<< matrixStack.size() <<endl;
            outputToFile(fout, matrixStack.top().first);
        }else{
            fout<<"Stack is Empty"<<endl;
        }

        if(command == "end"){
            break;
        }
        else if(command == "triangle"){
            matrix = triangle(fin, fout);
            fout<<"Triangle"<<endl;
            outputToFile(fout, matrix);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                fout << "After multiplication" << endl;
                outputToFile(fout, temp);
                fout << "Triangle Output format"<<endl;
                outputToFile(fout, transposeMatrix(temp));
            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "translate"){
            matrix = translate(fin, fout);
            fout<<"Translate matrix:"<<endl;
            outputToFile(fout, matrix);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                fout << "After multiplication" << endl;
                outputToFile(fout, temp);
                push(matrixStack, temp, false);
            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "scale"){
            matrix = scale(fin, fout);
            fout<<"Scale matrix:"<<endl;
            outputToFile(fout, matrix);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                fout << "After multiplication" << endl;
                outputToFile(fout, temp);
                push(matrixStack, temp, false);

            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }
        else if(command == "rotate"){
            matrix = rotate(fin, fout);
            fout<<"Rotate matrix:"<<endl;
            outputToFile(fout, matrix);
            if(!matrixStack.empty()){
                vector<vector<double>> temp = matrixMultiplication(matrixStack.top().first, matrix);
                fout << "After multiplication" << endl;
                outputToFile(fout, temp);
                push(matrixStack, temp, false);

            }else{
                fout<<"Error: stack is empty"<<endl;
            }
        }

        else if(command == "push"){
            matrix = matrixStack.top().first;
            matrixStack.pop();
            push(matrixStack, matrix, true);
            fout << "After push top of stack" << endl;
            fout<<"Top of stack | Flag " << matrixStack.top().second << " | Stack size = "<< matrixStack.size() <<endl;
            outputToFile(fout, matrixStack.top().first);
        }
        else if(command == "pop"){
            if(!matrixStack.empty()){
                pop(matrixStack, fout);
            }else{
                fout<<"Stack is Empty. Can't pop"<<endl;
            }    
        }
    }

    fin.close();
    fout.close();
}
