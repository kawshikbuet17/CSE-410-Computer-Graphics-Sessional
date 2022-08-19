#include<bits/stdc++.h>

using namespace std;

#define PI 2*acos(0.0)
#define INF 1e9

class Vector {
public:
    double x, y, z;
    Vector() {
        x = y = z = 0.0;
    }
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void normalize() {
        double len = sqrt(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0));
        x /= len;
        y /= len;
        z /= len;
    }
    Vector operator+(const Vector v) {
        return Vector(x + v.x, y + v.y, z + v.z);
    }

    Vector operator-(const Vector v) {
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    Vector operator*(const double scalar) {
        return Vector(x*scalar, y*scalar, z*scalar);
    }
};

double DOT(Vector v1, Vector v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

Vector CROSS(Vector v1, Vector v2){
    Vector v;
    v.x = v1.y*v2.z - v1.z*v2.y;
    v.y = v1.z*v2.x - v1.x*v2.z;
    v.z = v1.x*v2.y - v1.y*v2.x;
    return v;
}

double ValueOfVector(Vector v){
    return sqrt(pow(v.x, 2.0)+pow(v.y, 2.0)+pow(v.z, 2.0));
}

double distanceBetweenPoints(Vector p1, Vector p2){
    return sqrt(pow(p1.x-p2.x, 2.0)+pow(p1.y-p2.y, 2.0)+pow(p1.z-p2.z, 2.0));
}

class Ray {
public:
    Vector R0;
    Vector Rd;

    Ray(Vector R0, Vector Rd) {
        this->R0 = R0;
        this->Rd = Rd;
        this->Rd.normalize();
    }
};

class Color {
public:
    double red, green, blue;

    Color() {
        red = green = blue = 0;
    }

    Color(double red, double green, double blue) {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
};

class PointLight {
public:
    Vector position;
    Color color;
    double radius;
    int segments;
    int stacks;

    PointLight() {
        radius = 0.0;
        segments = stacks = 0;
    }
    PointLight(Vector position, Color color, double radius, int segments = 30, int stacks = 30) {
        this->position = position;
        this->color = color;
        this->radius = radius;
        this->segments = segments;
        this->stacks = stacks;
    }
    Vector getPosition(){
        return position;
    }
    Color getColor(){
        return color;
    }
    void draw() {
        Vector points[stacks+1][segments+1];
        double height, r;

        /* generating points: segments = segments in plane; stacks = segments in hemisphere */
        for(int i=0; i<=stacks; i++) {
            height = radius*sin(((double)i/(double)stacks)*(PI/2));
            r = radius * cos(((double)i / (double)stacks) * (PI / 2));

            for(int j=0; j<=segments; j++) {
                points[i][j] = Vector(r * cos(((double)j / (double)segments) * 2 * PI), r * sin(((double)j / (double)segments) * 2 * PI), height);
            }
        }

        /* drawing quads using generated points */
        glColor3f(color.red, color.green, color.blue);

        for(int i=0; i<stacks; i++) {
            for(int j=0; j<segments; j++) {
                glBegin(GL_QUADS);
                {
                    /* upper hemisphere */
                    glVertex3f((position+points[i][j]).x, (position+points[i][j]).y, (position+points[i][j]).z);
                    glVertex3f((position+points[i][j+1]).x, (position+points[i][j+1]).y, (position+points[i][j+1]).z);
                    glVertex3f((position+points[i+1][j+1]).x, (position+points[i+1][j+1]).y, (position+points[i+1][j+1]).z);
                    glVertex3f((position+points[i+1][j]).x, (position+points[i+1][j]).y, (position+points[i+1][j]).z);

                    /* lower hemisphere */
                    glVertex3f((position+points[i][j]).x, (position+points[i][j]).y, (position-points[i][j]).z);
                    glVertex3f((position+points[i][j+1]).x, (position+points[i][j+1]).y, (position-points[i][j+1]).z);
                    glVertex3f((position+points[i+1][j+1]).x, (position+points[i+1][j+1]).y, (position-points[i+1][j+1]).z);
                    glVertex3f((position+points[i+1][j]).x, (position+points[i+1][j]).y, (position-points[i+1][j]).z);
                }
                glEnd();
            }
        }
    }
};

class SpotLight {
public:
    Vector position;
    Color color;
    Vector direction;
    double cutoffAngle;
    double radius;
    int segments;
    int stacks;

    SpotLight() {
        radius = 0.0;
        segments = stacks = 0;
    }
    SpotLight(Vector position, Color color, Vector direction, double cutoffAngle, double radius, int segments = 30, int stacks = 30) {
        this->position = position;
        this->color = color;
        this->direction = direction;
        this->cutoffAngle = cutoffAngle;
        this->radius = radius;
        this->segments = segments;
        this->stacks = stacks;
    }
    Vector getPosition(){
        return position;
    }
    Color getColor(){
        return color;
    }
    void draw() {
        Vector points[stacks+1][segments+1];
        double height, r;

        /* generating points: segments = segments in plane; stacks = segments in hemisphere */
        for(int i=0; i<=stacks; i++) {
            height = radius*sin(((double)i/(double)stacks)*(PI/2));
            r = radius * cos(((double)i / (double)stacks) * (PI / 2));

            for(int j=0; j<=segments; j++) {
                points[i][j] = Vector(r * cos(((double)j / (double)segments) * 2 * PI), r * sin(((double)j / (double)segments) * 2 * PI), height);
            }
        }

        /* drawing quads using generated points */
        glColor3f(color.red, color.green, color.blue);

        for(int i=0; i<stacks; i++) {
            for(int j=0; j<segments; j++) {
                glBegin(GL_QUADS);
                {
                    /* upper hemisphere */
                    glVertex3f((position+points[i][j]).x, (position+points[i][j]).y, (position+points[i][j]).z);
                    glVertex3f((position+points[i][j+1]).x, (position+points[i][j+1]).y, (position+points[i][j+1]).z);
                    glVertex3f((position+points[i+1][j+1]).x, (position+points[i+1][j+1]).y, (position+points[i+1][j+1]).z);
                    glVertex3f((position+points[i+1][j]).x, (position+points[i+1][j]).y, (position+points[i+1][j]).z);

                    /* lower hemisphere */
                    glVertex3f((position+points[i][j]).x, (position+points[i][j]).y, (position-points[i][j]).z);
                    glVertex3f((position+points[i][j+1]).x, (position+points[i][j+1]).y, (position-points[i][j+1]).z);
                    glVertex3f((position+points[i+1][j+1]).x, (position+points[i+1][j+1]).y, (position-points[i+1][j+1]).z);
                    glVertex3f((position+points[i+1][j]).x, (position+points[i+1][j]).y, (position-points[i+1][j]).z);
                }
                glEnd();
            }
        }
    }
};


class ReflectionCoefficient {
public:
    double ambient, diffuse, specular, recursive;
    ReflectionCoefficient() {
        ambient = diffuse = specular = recursive = 0.0;
    }
    ReflectionCoefficient(double ambient, double diffuse, double specular, double recursive) {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->recursive = recursive;
    }
};

class Object {
public:
    Color color;
    ReflectionCoefficient reflectionCoefficient;
    int shine;
    Object() {
        shine = 0;
    }
    Color getColor(){
        return color;
    }
    void setColor(Color color) {
        this->color = color;
    }
    void setReflectionCoefficient(ReflectionCoefficient reflectionCoefficient) {
        this->reflectionCoefficient = reflectionCoefficient;
    }
    void setShine(int shine) {
        this->shine = shine;
    }
    virtual void draw() = 0;
    virtual double intersect(Ray*, Color&, int) = 0;
};

Vector pos, u, r, l;
int levelOfRecursion = 0;
vector<Object*> objects;
vector<PointLight> lights;
vector<SpotLight> spotlights;

class Sphere: public Object {
public:
    Vector center;
    double radius;
    int segments;
    int stacks;
    Sphere(Vector center, double radius, int segments = 30, int stacks = 30) {
        this->center = center;
        this->radius = radius;
        this->segments = segments;
        this->stacks = stacks;
    }

    void draw() {
        Vector points[stacks+1][segments+1];
        double height, r;

        /* generating points: segments = segments in plane; stacks = segments in hemisphere */
        for(int i=0; i<=stacks; i++) {
            height = radius*sin(((double)i/(double)stacks)*(PI/2));
            r = radius * cos(((double)i / (double)stacks) * (PI / 2));

            for(int j=0; j<=segments; j++) {
                points[i][j] = Vector(r * cos(((double)j / (double)segments) * 2 * PI), r * sin(((double)j / (double)segments) * 2 * PI), height);
            }
        }

        /* drawing quads using generated points */
        glColor3f(getColor().red, getColor().green, getColor().blue);

        for(int i=0; i<stacks; i++) {
            for(int j=0; j<segments; j++) {
                glBegin(GL_QUADS);
                {
                    /* upper hemisphere */
                    glVertex3f((center+points[i][j]).x, (center+points[i][j]).y, (center+points[i][j]).z);
                    glVertex3f((center+points[i][j+1]).x, (center+points[i][j+1]).y, (center+points[i][j+1]).z);
                    glVertex3f((center+points[i+1][j+1]).x, (center+points[i+1][j+1]).y, (center+points[i+1][j+1]).z);
                    glVertex3f((center+points[i+1][j]).x, (center+points[i+1][j]).y, (center+points[i+1][j]).z);

                    /* lower hemisphere */
                    glVertex3f((center+points[i][j]).x, (center+points[i][j]).y, (center-points[i][j]).z);
                    glVertex3f((center+points[i][j+1]).x, (center+points[i][j+1]).y, (center-points[i][j+1]).z);
                    glVertex3f((center+points[i+1][j+1]).x, (center+points[i+1][j+1]).y, (center-points[i+1][j+1]).z);
                    glVertex3f((center+points[i+1][j]).x, (center+points[i+1][j]).y, (center-points[i+1][j]).z);
                }
                glEnd();
            }
        }
    }
    double intersect(Ray* ray, Color& color, int level) {
        //find intersect value t minimum of at^2+bt+c=0
        double a, b, c, tMin;

        a = DOT(ray->Rd, ray->Rd);
        b = 2.0 * (DOT(ray->R0, ray->Rd) - DOT(ray->Rd, center));
        c = DOT(ray->R0, ray->R0) + DOT(center, center) - DOT(ray->R0, center) * 2.0 - radius * radius;

        double discriminant = b * b - 4.0 * a * c;

        if(discriminant < 0.0) {
            tMin = INF;
        }
        else{
            double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
            tMin = min(t1, t2);
        }

        if(level == 0) {
            return tMin;
        }

        //Illumination with the Phong Lighting Model
        Vector intersectionPoint = ray->R0 + ray->Rd * tMin;
        Color intersectionPointColor = getColor();

        //ambient component of reflected ray
        color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
        color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
        color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;

        //calculate normal at intersectionPoint
        Vector normal = intersectionPoint-center;
        normal.normalize();
        if(distanceBetweenPoints(pos, center) < radius){
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        //for each point light pl in pointLights
        for(int i=0; i<lights.size(); i++) {
            //cast rayl from pl.light_pos to intersectionPoint
            Ray* incidentRay = new Ray(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());


            double t, tMinimum = INF;

            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;  // color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }

            Vector shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd*tMinimum;
            double epsilon = 0.0000001;  // for tuning light effect

            // if intersectionPoint is in shadow, the diffuse
            // and specular components need not be calculated
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }

            //calculate lambertValue using normal, rayl
            //find reflected ray, rayr for rayl
            //calculate phongValue using r, rayr
            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        }
        for (int i = 0; i < spotlights.size(); i++) {
            Ray* incidentRay = new Ray(spotlights[i].getPosition(), spotlights[i].direction);
            double angle = acos(DOT(intersectionPoint, spotlights[i].getPosition())/(ValueOfVector(intersectionPoint)*
                                                                                     ValueOfVector(spotlights[i].getPosition())) * PI / 180);

            if (angle > spotlights[i].cutoffAngle) continue;
            /* checking if intersection point is in shadow */
            double t, tMinimum = INF;
            for (int j = 0; j < objects.size(); j++) {
                Color dummyColor;  // light_color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);
                if (t <= 0) continue;
                if (t < tMinimum) {
                    tMinimum = t;
                }
            }
            Vector shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd * tMinimum;
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }

            //calculate lambertValue using normal, rayl
            //find reflected ray, rayr for rayl
            //calculate phongValue using r, rayr
            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        }

        //if level ≥ recursion_level, return tmin
        if(level >= levelOfRecursion) {
            return tMin;
        }

        //construct reflected ray from intersection point
        Vector reflectionDirection = ray->Rd-normal*(DOT(ray->Rd,normal)*2.0);
        reflectionDirection.normalize();
        Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

        //find tmin from the nearest intersecting object, using
        //intersect() method, as done in the capture() method
        //if found, call intersect(rreflected, colorreflected, level+1)
        int nearest = INF;
        double t, tMinimum=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;  // color = black
            t = objects[i]->intersect(reflectedRay, dummyColor, 0);

            if(t>0.0 and t<tMinimum) {
                tMinimum = t;
                nearest = i;
            }
        }

        // colorreflected will be updated while in the subsequent call
        // update color using the impact of reflection
        Color reflectedColor;
        if(nearest != INF) {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
        }

        color.red += reflectedColor.red*reflectionCoefficient.recursive;
        color.green += reflectedColor.green*reflectionCoefficient.recursive;
        color.blue += reflectedColor.blue*reflectionCoefficient.recursive;
        return tMin;
    }
};

class Triangle: public Object {
public:
    Vector a;
    Vector b;
    Vector c;
    Triangle(Vector a, Vector b, Vector c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
    void draw() {
        glColor3f(getColor().red, getColor().green, getColor().blue);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }
    double intersect(Ray* ray, Color& color, int level) {
        double detA, detBeta, detGhama, detT, tMin;
        double R0x, R0y, R0z;
        double Rdx, Rdy, Rdz;
        R0x = ray->R0.x;
        R0y = ray->R0.y;
        R0z = ray->R0.z;
        Rdx = ray->Rd.x;
        Rdy = ray->Rd.y;
        Rdz = ray->Rd.z;
        detA = (a.x - b.x) * ((a.y - c.y) * Rdz - (a.z - c.z) * Rdy)
               +(a.x - c.x) * ((a.z - b.z) * Rdy - (a.y - b.y) * Rdz)
               +Rdx * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));

        detBeta = (a.x - R0x) * ((a.y - c.y) * Rdz - (a.z - c.z) * Rdy)
                  +(a.x - c.x) * ((a.z - R0z) * Rdy - (a.y - R0y) * Rdz)
                  +Rdx * ((a.y - R0y) * (a.z - c.z) - (a.z - R0z) * (a.y - c.y));

        detGhama = (a.x - b.x) * ((a.y - R0y) * Rdz - (a.z - R0z) * Rdy)
                   +(a.x - R0x) * ((a.z - b.z) * Rdy - (a.y - b.y) * Rdz)
                   +Rdx * ((a.y - b.y) * (a.z - R0z) - (a.z - b.z) * (a.y - R0y));

        detT = (a.x - b.x) * ((a.y - c.y) * (a.z - R0z) - (a.z - c.z) * (a.y - R0y))
               +(a.x - c.x) * ((a.z - b.z) * (a.y - R0y) - (a.y - b.y) * (a.z - R0z))
               +(a.x - R0x) * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));

        tMin = INF;
        if(detA != 0.0) {
            double beta = detBeta / detA;
            double ghama = detGhama / detA;
            double t = detT / detA;
            if(beta + ghama < 1.0 and beta > 0.0 && ghama > 0.0) {
                tMin = t;
            }
        }

        if(level == 0) {
            return tMin;
        }

        //Illumination with the Phong Lighting Model
        Vector intersectionPoint = ray->R0 + ray->Rd * tMin;
        Color intersectionPointColor = getColor();

        //ambient component of reflected ray
        color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
        color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
        color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;

        //calculate normal
        Vector normal = CROSS((b-a),(c-a));
        normal.normalize();
        if(DOT(ray->Rd, normal) > 0.0){
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        //for each point light pl in pointLights
        for(int i=0; i<lights.size(); i++) {
            //cast rayl from pl.light_pos to intersectionPoint
            Ray* incidentRay = new Ray(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

            double t, tMinimum=INF;

            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;  // color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }

            Vector shadowIntersectionPoint = incidentRay->R0+incidentRay->Rd*tMinimum;

            // if intersectionPoint is in shadow, the diffuse
            // and specular components need not be calculated
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) - (1e-7) > distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }

            //calculate lambertValue using normal, rayl
            //find reflected ray, rayr for rayl
            //calculate phongValue using r, rayr
            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        }
        for (int i = 0; i < spotlights.size(); i++) {
            Ray* incidentRay = new Ray(spotlights[i].getPosition(), spotlights[i].direction);
            double angle = acos(DOT(intersectionPoint, spotlights[i].getPosition())/(ValueOfVector(intersectionPoint)*
                                                                                     ValueOfVector(spotlights[i].getPosition())) * PI / 180);

            if (angle > spotlights[i].cutoffAngle) continue;
            /* checking if intersection point is in shadow */
            double t, tMinimum = INF;
            for (int j = 0; j < objects.size(); j++) {
                Color dummyColor;  // light_color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);
                if (t <= 0) continue;
                if (t < tMinimum) {
                    tMinimum = t;
                }
            }
            Vector shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd * tMinimum;
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) - (1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }

            //calculate lambertValue using normal, rayl
            //find reflected ray, rayr for rayl
            //calculate phongValue using r, rayr
            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        }

        //if level ≥ recursion_level, return tmin
        if(level >= levelOfRecursion) {
            return tMin;
        }

        //construct reflected ray from intersection point
        Vector reflectionDirection = ray->Rd-normal*(DOT(ray->Rd, normal)*2.0);
        reflectionDirection.normalize();
        Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

        //find tmin from the nearest intersecting object, using
        //intersect() method, as done in the capture() method
        //if found, call intersect(rreflected, colorreflected, level+1)
        int nearest = INF;
        double t, tMinimum=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;  // color = black
            t = objects[i]->intersect(reflectedRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
                nearest = i;
            }
        }

        // colorreflected will be updated while in the subsequent call
        // update color using the impact of reflection
        Color reflectedColor;  // color = black

        if(nearest != INF) {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
        }

        color.red += reflectedColor.red*reflectionCoefficient.recursive;
        color.green += reflectedColor.green*reflectionCoefficient.recursive;
        color.blue += reflectedColor.blue*reflectionCoefficient.recursive;

        return tMin;
    }
};


class GeneralQuadricSurfaceCoefficient {
public:
    double a, b, c, d, e, f, g, h, i, j;
};

class GeneralQuadricSurface: public Object {
public:
    GeneralQuadricSurfaceCoefficient coefficient;
    Vector cubeReferencePoint;
    double length;
    double width;
    double height;
    GeneralQuadricSurface() {
        coefficient.a = coefficient.b = coefficient.c = coefficient.d = coefficient.e = coefficient.f = coefficient.g = coefficient.h = coefficient.i = coefficient.j = 0.0;
        length = width = height = 0.0;
    }

    GeneralQuadricSurface(GeneralQuadricSurfaceCoefficient coefficient, Vector cubeReferencePoint, double length, double width, double height) {
        this->coefficient = coefficient;
        this->cubeReferencePoint = cubeReferencePoint;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    void draw(){
        ;
    }

    double intersect(Ray* ray, Color& color, int level) {
        double a, b, c, tMin, tMax;

        a = coefficient.a*ray->Rd.x*ray->Rd.x+coefficient.b*ray->Rd.y*ray->Rd.y + coefficient.c*ray->Rd.z*ray->Rd.z + coefficient.d*ray->Rd.x*ray->Rd.y + coefficient.e*ray->Rd.x*ray->Rd.z+coefficient.f*ray->Rd.y*ray->Rd.z;

        b = 2.0*coefficient.a*ray->R0.x*ray->Rd.x+2.0*coefficient.b*ray->R0.y*ray->Rd.y+2.0*coefficient.c*ray->R0.z*ray->Rd.z;
        b += coefficient.d*(ray->R0.x*ray->Rd.y+ray->Rd.x*ray->R0.y);
        b += coefficient.e*(ray->R0.x*ray->Rd.z+ray->Rd.x*ray->R0.z);
        b += coefficient.f*(ray->R0.y*ray->Rd.z+ray->Rd.y*ray->R0.z);
        b += coefficient.g*ray->Rd.x+coefficient.h*ray->Rd.y+coefficient.i*ray->Rd.z;

        c = coefficient.a*ray->R0.x*ray->R0.x+coefficient.b*ray->R0.y*ray->R0.y+coefficient.c*ray->R0.z*ray->R0.z;
        c += coefficient.d*ray->R0.x*ray->R0.y+coefficient.e*ray->R0.x*ray->R0.z+coefficient.f*ray->R0.y*ray->R0.z;
        c += coefficient.g*ray->R0.x+coefficient.h*ray->R0.y+coefficient.i*ray->R0.z+coefficient.j;

        if(a == 0.0) {
            if(b == 0.0){
                tMin = INF;
            }
            else{
                tMin = -c/b;
            }
            tMax = INF;
        }
        else {
            double discriminant = b * b - 4.0 * a * c;

            if(discriminant < 0.0) {
                tMin = tMax = INF;
            }
            else {
                tMax = -b/(2.0*a)+sqrt(discriminant)/(2.0*a);
                tMin = -b/(2.0*a)-sqrt(discriminant)/(2.0*a);
            }
        }

        /* clipping general quadric surface along the dimensions (if necessary) */
        if(tMin < INF) {
            if(tMax < INF) {
                if(tMin > 0.0) {
                    Vector intersectionPoint = ray->R0+ray->Rd*tMin;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                        tMin = INF;
                    }
                }
                if(tMax > 0.0) {
                    Vector intersectionPoint = ray->R0+ray->Rd*tMax;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                        tMax = INF;
                    }
                }
                tMin = (tMin>0.0 && tMin<tMax)? tMin: tMax;
            } else {
                if(tMin > 0.0) {
                    Vector intersectionPoint = ray->R0+ray->Rd*tMin;

                    if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                        tMin = INF;
                    }
                }
            }
        }

        if(level == 0) {
            return tMin;
        }

        Vector intersectionPoint = ray->R0+ray->Rd*tMin;
        Color intersectionPointColor = getColor();

        //ambient component
        color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
        color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
        color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;

        double xNormal, yNormal, zNormal;
        xNormal = 2.0*coefficient.a*intersectionPoint.x+coefficient.d*intersectionPoint.y;
        xNormal += coefficient.e*intersectionPoint.z+coefficient.g;

        yNormal = 2.0*coefficient.b*intersectionPoint.y+coefficient.d*intersectionPoint.x;
        yNormal += coefficient.f*intersectionPoint.z+coefficient.h;

        zNormal = 2.0*coefficient.c*intersectionPoint.z+coefficient.e*intersectionPoint.x;
        zNormal += coefficient.f*intersectionPoint.y+coefficient.i;

        Vector normal(xNormal, yNormal, zNormal);
        normal.normalize();


        for(int i=0; i<lights.size(); i++) {
            Ray* incidentRay= new Ray(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;  // color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }

            Vector shadowIntersectionPoint = incidentRay->R0+incidentRay->Rd*tMinimum;
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }


            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)), reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red*intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
            color.green += lights[i].getColor().green*intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
        }
        for (int i = 0; i < spotlights.size(); i++) {
            Ray* incidentRay = new Ray(spotlights[i].getPosition(), spotlights[i].direction);
            double angle = acos(DOT(intersectionPoint, spotlights[i].getPosition())/(ValueOfVector(intersectionPoint)*
                                                                                     ValueOfVector(spotlights[i].getPosition())) * PI / 180);

            if (angle > spotlights[i].cutoffAngle) continue;
            /* checking if intersection point is in shadow */
            double t, tMinimum = INF;
            for (int j = 0; j < objects.size(); j++) {
                Color dummyColor;  // light_color = black
                t = objects[j]->intersect(incidentRay, dummyColor, 0);
                if (t <= 0) continue;
                if (t < tMinimum) {
                    tMinimum = t;
                }
            }
            Vector shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd * tMinimum;
            double epsilon = 1e-7;  // for tuning light effect
            if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7) > distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
                continue;
            }

            //calculate lambertValue using normal, rayl
            //find reflected ray, rayr for rayl
            //calculate phongValue using r, rayr
            double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
            Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
            double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
            color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

            color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
            color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        }

        if(level >= levelOfRecursion) {
            return tMin;
        }

        Vector reflectionDirection = ray->Rd-normal*(DOT(ray->Rd, normal)*2.0);
        reflectionDirection.normalize();
        Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

        int nearest = INF;
        double t, tMinimum=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;  // color = black
            t = objects[i]->intersect(reflectedRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
                nearest = i;
            }
        }

        Color reflectedColor;

        if(nearest != INF) {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
        }

        color.red += reflectedColor.red*reflectionCoefficient.recursive;
        color.green += reflectedColor.green*reflectionCoefficient.recursive;
        color.blue += reflectedColor.blue*reflectionCoefficient.recursive;

        return tMin;
    }
};

class Floor: public Object {
    double floorWidth;
    double tileWidth;
    Color foregroundColor;

public:
    Floor() {
        floorWidth = tileWidth = 0.0;
    }

    Floor(double floorWidth, double tileWidth, Color foregroundColor) {
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
        this->foregroundColor = foregroundColor;
    }

    void draw();
    double intersect(Ray*, Color&, int);
};

void Floor::draw() {
    for(int i=0, row=(int) floorWidth/tileWidth, column=(int) floorWidth/tileWidth; i<row; i++) {
        for(int j=0; j<column; j++) {
            glColor3f(((i+j)%2 == 0)? getColor().red: foregroundColor.red, ((i+j)%2 == 0)? getColor().green: foregroundColor.green, ((i+j)%2 == 0)? getColor().blue: foregroundColor.blue);
            Vector leftBottomCorner(-floorWidth/2.0+tileWidth*j, -floorWidth/2.0+tileWidth*i, 0.0);
            glBegin(GL_QUADS);
            {
                glVertex3f(leftBottomCorner.x, leftBottomCorner.y, leftBottomCorner.z);
                glVertex3f(leftBottomCorner.x+tileWidth, leftBottomCorner.y, leftBottomCorner.z);
                glVertex3f(leftBottomCorner.x+tileWidth, leftBottomCorner.y+tileWidth, leftBottomCorner.z);
                glVertex3f(leftBottomCorner.x, leftBottomCorner.y+tileWidth, leftBottomCorner.z);
            }
            glEnd();
        }
    }
}

double Floor::intersect(Ray* ray, Color& color, int level) {
    //Normal = (0,0,1)
    Vector normal(0.0, 0.0, 1.0);

    //t = -(D + n·Ro) / n·Rd
    double tMin = INF;
    double D = 0;
    double nR0 = DOT(normal, ray->R0);
    double nRd = DOT(normal, ray->Rd);
    if(nRd != 0.0) {
        tMin = -(D + nR0) / nRd;
    }

    if(level == 0) {
        return tMin;
    }

    Vector intersectionPoint = ray->R0+ray->Rd*tMin;
    Vector referencePosition = intersectionPoint-Vector(-floorWidth/2.0, -floorWidth/2.0, 0.0);
    Color intersectionPointColor;
    if(((int) (floor(referencePosition.x/tileWidth)+floor(referencePosition.y/tileWidth)))%2 == 0){
        intersectionPointColor = getColor();
    }
    else{
        intersectionPointColor = foregroundColor;
    }

    //ambient component
    color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
    color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
    color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;

    for(int i=0; i<lights.size(); i++) {
        Ray* incidentRay = new Ray(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        double t, tMinimum=INF;
        for(int j=0; j<objects.size(); j++) {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
            }
        }

        Vector shadowIntersectionPoint = incidentRay->R0+incidentRay->Rd*tMinimum;

        if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) -(1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
            continue;
        }

        double lambertValue = DOT((incidentRay->Rd*(-1.0)), normal);
        Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
        double phongValue = DOT((ray->Rd*(-1.0)), reflectedRay->Rd);

        color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
        color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
        color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

        color.red += lights[i].getColor().red*intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
        color.green += lights[i].getColor().green*intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
        color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0),shine);
    }
    for (int i = 0; i < spotlights.size(); i++) {
        Ray* incidentRay = new Ray(spotlights[i].getPosition(), spotlights[i].direction);
        double angle = acos(DOT(intersectionPoint, spotlights[i].getPosition())/(ValueOfVector(intersectionPoint)*
                                                                                 ValueOfVector(spotlights[i].getPosition()))) * 180 / PI;

        if (angle > spotlights[i].cutoffAngle) continue;
        /* checking if intersection point is in shadow */
        double t, tMinimum = INF;
        for (int j = 0; j < objects.size(); j++) {
            Color dummyColor;  // light_color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);
            if (t <= 0) continue;
            if (t < tMinimum) {
                tMinimum = t;
            }
        }
        Vector shadowIntersectionPoint = incidentRay->R0 + incidentRay->Rd * tMinimum;
        if(distanceBetweenPoints(intersectionPoint, incidentRay->R0) - (1e-7)> distanceBetweenPoints(shadowIntersectionPoint, incidentRay->R0)){
            continue;
        }

        //calculate lambertValue using normal, rayl
        //find reflected ray, rayr for rayl
        //calculate phongValue using r, rayr
        double lambertValue = DOT((incidentRay->Rd*(-1.0)),normal);
        Ray* reflectedRay = new Ray(intersectionPoint, incidentRay->Rd-normal*(DOT(incidentRay->Rd, normal)*2.0));
        double phongValue = DOT((ray->Rd*(-1.0)),reflectedRay->Rd);

        color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.diffuse * max(lambertValue, 0.0);
        color.green += lights[i].getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse * max(lambertValue, 0.0);
        color.blue += lights[i].getColor().blue*intersectionPointColor.blue * reflectionCoefficient.diffuse * max(lambertValue, 0.0);

        color.red += lights[i].getColor().red * intersectionPointColor.red * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        color.green += lights[i].getColor().green * intersectionPointColor.green * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
        color.blue += lights[i].getColor().blue * intersectionPointColor.blue * reflectionCoefficient.specular * pow(max(phongValue, 0.0), shine);
    }

    if(level >= levelOfRecursion) {
        return tMin;
    }

    Vector reflectionDirection = ray->Rd-normal*(DOT(ray->Rd,normal)*2.0);
    reflectionDirection.normalize();
    Ray* reflectedRay = new Ray(intersectionPoint+reflectionDirection, reflectionDirection);

    int nearest = INF;
    double t, tMinimum=INF;

    for(int i=0; i<objects.size(); i++) {
        Color dummyColor;
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 and t<tMinimum) {
            tMinimum = t;
            nearest = i;
        }
    }

    Color reflectedColor;

    if(nearest != INF) {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    color.red += reflectedColor.red*reflectionCoefficient.recursive;
    color.green += reflectedColor.green*reflectionCoefficient.recursive;
    color.blue += reflectedColor.blue*reflectionCoefficient.recursive;

    return tMin;
}