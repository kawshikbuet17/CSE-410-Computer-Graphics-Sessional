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
    void normalize();
    double computeDistanceBetween(Vector);
    Vector operator+(const Vector);
    Vector operator-(const Vector);
    Vector operator*(const double);
    double operator/(const Vector);
    Vector operator%(const Vector);

    friend ostream& operator<<(ostream&, Vector&);
};

void Vector::normalize() {
    double magnitude = sqrt(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0));
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
}

double Vector::computeDistanceBetween(Vector _vector) {
    return sqrt(pow(x-_vector.x, 2.0)+pow(y-_vector.y, 2.0)+pow(z-_vector.z, 2.0));
}

Vector Vector::operator+(const Vector _vector) {
    return Vector(x+_vector.x, y+_vector.y, z+_vector.z);
}

Vector Vector::operator-(const Vector _vector) {
    return Vector(x-_vector.x, y-_vector.y, z-_vector.z);
}

Vector Vector::operator*(const double scalar) {
    /* vector scalar multiplication */
    return Vector(x*scalar, y*scalar, z*scalar);
}

double Vector::operator/(const Vector _vector) {
    /* vector dot multiplication */
    return x*_vector.x+y*_vector.y+z*_vector.z;
}

Vector Vector::operator%(const Vector _vector) {
    /* vector cross multiplication */
    return Vector(y*_vector.z-z*_vector.y, z*_vector.x-x*_vector.z, x*_vector.y-y*_vector.x);
}

ostream& operator<<(ostream &output, Vector &_vector) {
    output << '(' << _vector.x << ", " << _vector.y << ", " << _vector.z << ')';
    return output;
}

double DOT(Vector v1, Vector v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
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
    Vector position;
    Color color;
    double radius;
    int segments;
    int stacks;

public:
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

    Vector getPosition() const {
        return position;
    }

    Color getColor() const {
        return color;
    }

    void draw();
};

void PointLight::draw() {
    Vector points[stacks+1][segments+1];
    double height, _radius;

    /* generating points: segments = segments in plane; stacks = segments in hemisphere */
	for(int i=0; i<=stacks; i++) {
		height = radius*sin(((double)i/(double)stacks)*(PI/2));
		_radius = radius*cos(((double)i/(double)stacks)*(PI/2));

		for(int j=0; j<=segments; j++) {
            points[i][j] = Vector(_radius*cos(((double)j/(double)segments)*2*PI), _radius*sin(((double)j/(double)segments)*2*PI), height);
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

    Color getColor() const {
        return color;
    }

    void setColor(Color color) {
        this->color = color;
    }

    void setReflectionCoefficient(ReflectionCoefficient reflectionCoefficient) {
        this->reflectionCoefficient = reflectionCoefficient;
    }

    int getShine() const {
        return shine;
    }

    void setShine(int shine) {
        this->shine = shine;
    }

    void computeAmbientLightComponent(Color&, Color);
    void computeReflectionComponents(Ray, Color&, Vector, Color, Vector, PointLight, Ray);
    void computeRecursiveReflectionComponent(Color&, Color);

    virtual void draw() = 0;
    virtual double intersect(Ray, Color&, int) = 0;
};

void Object::computeAmbientLightComponent(Color& color, Color intersectionPointColor) {
    color.red = intersectionPointColor.red*reflectionCoefficient.ambient;
    color.green = intersectionPointColor.green*reflectionCoefficient.ambient;
    color.blue = intersectionPointColor.blue*reflectionCoefficient.ambient;
}

void Object::computeReflectionComponents(Ray ray, Color& color, Vector intersectionPoint, Color intersectionPointColor, Vector normal, PointLight light, Ray incidentRay) {
    double lambertValue = (incidentRay.Rd*(-1.0))/normal;
    Ray reflectedRay(intersectionPoint, incidentRay.Rd-normal*((incidentRay.Rd/normal)*2.0));
    double phongValue = (ray.Rd*(-1.0))/reflectedRay.Rd;

    color.red += light.getColor().red*intersectionPointColor.red*reflectionCoefficient.diffuse*max(lambertValue, 0.0);
    color.green += light.getColor().green*intersectionPointColor.green*reflectionCoefficient.diffuse*max(lambertValue, 0.0);
    color.blue += light.getColor().blue*intersectionPointColor.blue*reflectionCoefficient.diffuse*max(lambertValue, 0.0);

    color.red += light.getColor().red*intersectionPointColor.red*reflectionCoefficient.specular*pow(max(phongValue, 0.0),
                                                                                                    getShine());
    color.green += light.getColor().green*intersectionPointColor.green*reflectionCoefficient.specular*pow(max(phongValue, 0.0),
                                                                                                          getShine());
    color.blue += light.getColor().blue*intersectionPointColor.blue*reflectionCoefficient.specular*pow(max(phongValue, 0.0),
                                                                                                       getShine());
}

void Object::computeRecursiveReflectionComponent(Color& color, Color reflectedColor) {
    color.red += reflectedColor.red*reflectionCoefficient.recursive;
    color.green += reflectedColor.green*reflectionCoefficient.recursive;
    color.blue += reflectedColor.blue*reflectionCoefficient.recursive;
}

Vector pos, u, r, l;

int levelOfRecursion = 0;

vector<Object*> objects;
vector<PointLight> lights;

class Sphere: public Object {
    Vector center;
    double radius;
    int segments;
    int stacks;

public:
    Sphere(Vector center, double radius, int segments = 30, int stacks = 30) {
        this->center = center;
        this->radius = radius;

        this->segments = segments;
        this->stacks = stacks;
    }

    void draw();
    double intersect(Ray, Color&, int);
};

void Sphere::draw() {
    Vector points[stacks+1][segments+1];
    double height, _radius;

    /* generating points: segments = segments in plane; stacks = segments in hemisphere */
	for(int i=0; i<=stacks; i++) {
		height = radius*sin(((double)i/(double)stacks)*(PI/2));
		_radius = radius*cos(((double)i/(double)stacks)*(PI/2));

		for(int j=0; j<=segments; j++) {
            points[i][j] = Vector(_radius*cos(((double)j/(double)segments)*2*PI), _radius*sin(((double)j/(double)segments)*2*PI), height);
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

double Sphere::intersect(Ray ray, Color& color, int level) {
    //find intersect value t minimum of at^2+bt+c=0
    double a, b, c, tMin;

    a = DOT(ray.Rd, ray.Rd);
    b = 2.0 * (DOT(ray.R0, ray.Rd) - DOT(ray.Rd, center));
    c = DOT(ray.R0, ray.R0) + DOT(center, center) - DOT(ray.R0, center) * 2.0 - radius * radius;

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
    Vector intersectionPoint = ray.R0 + ray.Rd * tMin;
    Color intersectionPointColor = getColor();

    //ambient component of reflected ray
    computeAmbientLightComponent(color, intersectionPointColor);

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
        Ray incidentRay(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        // if intersectionPoint is in shadow, the diffuse
        // and specular components need not be calculated
        double t, tMinimum=INF;

        for(int j=0; j<objects.size(); j++) {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.R0 + incidentRay.Rd*tMinimum;
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectionPoint.computeDistanceBetween(incidentRay.R0)-epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.R0)) {
            /* intersection point is, indeed, in shadow */
            continue;
        }

        /* computing diffuse & specular components of reflected ray */
        computeReflectionComponents(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    /* handling recursive reflection */
    if(level >= levelOfRecursion) {
        return tMin;
    }

    /* incorporating concept of evil epsilon to recursive reflection computation */
    Vector reflectionDirection = ray.Rd-normal*((ray.Rd/normal)*2.0);
    reflectionDirection.normalize();
    Ray reflectedRay(intersectionPoint+reflectionDirection, reflectionDirection);

    /* finding nearest intersecting object (if available) */
    int nearest = INT_MAX;
    double t, tMinimum=INF;

    for(int i=0; i<objects.size(); i++) {
        Color dummyColor;  // color = black
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 && t<tMinimum) {
            tMinimum = t;
            nearest = i;
        }
    }

    /* finding color component for reflected ray */
    Color reflectedColor;  // color = black

    if(nearest != INT_MAX) {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    /* computing recursive reflection component of reflected ray */
    computeRecursiveReflectionComponent(color, reflectedColor);

    /* clipping the color values (if necessary) */
    color.red = (color.red > 1.0)? 1.0: ((color.red < 0.0)? 0.0: color.red);
    color.green = (color.green > 1.0)? 1.0: ((color.green < 0.0)? 0.0: color.green);
    color.blue = (color.blue > 1.0)? 1.0: ((color.blue < 0.0)? 0.0: color.blue);

    return tMin;
}

/* Triangle class */
class Triangle: public Object {
    Vector a;
    Vector b;
    Vector c;

public:
    Triangle() {
        /* default constructor */
    }

    Triangle(Vector a, Vector b, Vector c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void draw();
    double intersect(Ray, Color&, int);

    ~Triangle() {
        /* destructor */
    }
};

void Triangle::draw() {
    /* a, b, c - coordinates/position vectors of three corners of the triangle */
    glColor3f(getColor().red, getColor().green, getColor().blue);

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
    }
    glEnd();
}

double Triangle::intersect(Ray ray, Color& color, int level) {
    /* finding intersecting tMin */
    double determinantBase, determinantBeta, determinantGamma, determinantT, tMin;

    determinantBase = (a.x-b.x)*((a.y-c.y)*ray.Rd.z-(a.z-c.z)*ray.Rd.y);
    determinantBase += (a.x-c.x)*((a.z-b.z)*ray.Rd.y-(a.y-b.y)*ray.Rd.z);
    determinantBase += ray.Rd.x*((a.y-b.y)*(a.z-c.z)-(a.z-b.z)*(a.y-c.y));

    determinantBeta = (a.x-ray.R0.x)*((a.y-c.y)*ray.Rd.z-(a.z-c.z)*ray.Rd.y);
    determinantBeta += (a.x-c.x)*((a.z-ray.R0.z)*ray.Rd.y-(a.y-ray.R0.y)*ray.Rd.z);
    determinantBeta += ray.Rd.x*((a.y-ray.R0.y)*(a.z-c.z)-(a.z-ray.R0.z)*(a.y-c.y));

    determinantGamma = (a.x-b.x)*((a.y-ray.R0.y)*ray.Rd.z-(a.z-ray.R0.z)*ray.Rd.y);
    determinantGamma += (a.x-ray.R0.x)*((a.z-b.z)*ray.Rd.y-(a.y-b.y)*ray.Rd.z);
    determinantGamma += ray.Rd.x*((a.y-b.y)*(a.z-ray.R0.z)-(a.z-b.z)*(a.y-ray.R0.y));

    determinantT = (a.x-b.x)*((a.y-c.y)*(a.z-ray.R0.z)-(a.z-c.z)*(a.y-ray.R0.y));
    determinantT += (a.x-c.x)*((a.z-b.z)*(a.y-ray.R0.y)-(a.y-b.y)*(a.z-ray.R0.z));
    determinantT += (a.x-ray.R0.x)*((a.y-b.y)*(a.z-c.z)-(a.z-b.z)*(a.y-c.y));

    if(determinantBase == 0.0) {
        /* ray will not intersect the triangle plane */
        tMin = INF;
    } else {
        /* ray will intersect the triangle plane */
        if(determinantBeta/determinantBase>0.0 && determinantGamma/determinantBase>0.0 && determinantBeta/determinantBase+determinantGamma/determinantBase<1.0) {
            /* intersection point lies within the boundary of the triangle */
            tMin = determinantT/determinantBase;
        } else {
            /* intersection point does not lie within the boundary of the triangle */
            tMin = INF;
        }
    }

    if(level == 0) {
        return tMin;
    }

    /* illuminating with Phong Lighting Model */
    Vector intersectionPoint = ray.R0+ray.Rd*tMin;
    Color intersectionPointColor = getColor();

    /* determining unit normal vector on appropriate side of triangle */
    Vector normal = (b-a)%(c-a);
    normal.normalize();
    normal = ((ray.Rd*(-1.0))/normal > 0.0)? normal: normal*(-1.0);

    /* computing ambient light component of reflected ray */
    computeAmbientLightComponent(color, intersectionPointColor);

    /* computing diffuse & specular reflection components of reflected ray */
    for(int i=0; i<lights.size(); i++) {
        Ray incidentRay(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        /* checking if intersection point is in shadow */
        double t, tMinimum=INF;

        for(int j=0; j<objects.size(); j++) {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.R0+incidentRay.Rd*tMinimum;
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectionPoint.computeDistanceBetween(incidentRay.R0)-epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.R0)) {
            /* intersection point is, indeed, in shadow */
            continue;
        }

        /* computing diffuse & specular components of reflected ray */
        computeReflectionComponents(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    /* handling recursive reflection */
    if(level >= levelOfRecursion) {
        return tMin;
    }

    /* incorporating concept of evil epsilon to recursive reflection computation */
    Vector reflectionDirection = ray.Rd-normal*((ray.Rd/normal)*2.0);
    reflectionDirection.normalize();
    Ray reflectedRay(intersectionPoint+reflectionDirection, reflectionDirection);

    /* finding nearest intersecting object (if available) */
    int nearest = INT_MAX;
    double t, tMinimum=INF;

    for(int i=0; i<objects.size(); i++) {
        Color dummyColor;  // color = black
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 && t<tMinimum) {
            tMinimum = t;
            nearest = i;
        }
    }

    /* finding color component for reflected ray */
    Color reflectedColor;  // color = black

    if(nearest != INT_MAX) {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    /* computing recursive reflection component of reflected ray */
    computeRecursiveReflectionComponent(color, reflectedColor);

    /* clipping the color values (if necessary) */
    color.red = (color.red > 1.0)? 1.0: ((color.red < 0.0)? 0.0: color.red);
    color.green = (color.green > 1.0)? 1.0: ((color.green < 0.0)? 0.0: color.green);
    color.blue = (color.blue > 1.0)? 1.0: ((color.blue < 0.0)? 0.0: color.blue);

    return tMin;
}

/* GeneralQuadricSurfaceCoefficient structure */
struct GeneralQuadricSurfaceCoefficient {
    double a, b, c, d, e, f, g, h, i, j;

    friend ifstream& operator>>(ifstream&, GeneralQuadricSurfaceCoefficient&);
    friend ostream& operator<<(ostream&, GeneralQuadricSurfaceCoefficient&);
};

ifstream& operator>>(ifstream &input, GeneralQuadricSurfaceCoefficient &coefficient) {
    input >> coefficient.a >> coefficient.b >> coefficient.c >> coefficient.d >> coefficient.e;
    input >> coefficient.f >> coefficient.g >> coefficient.h >> coefficient.i >> coefficient.j;
    return input;
}

ostream& operator<<(ostream &output, GeneralQuadricSurfaceCoefficient &coefficient) {
    output << '[' << coefficient.a << ", " << coefficient.b << ", ";
    output << coefficient.c << ", " << coefficient.d << ", ";
    output << coefficient.e << ", " << coefficient.f << ", ";
    output << coefficient.g << ", " << coefficient.h << ", ";
    output << coefficient.i << ", " << coefficient.j << ']';
    return output;
}

/* GeneralQuadricSurface class */
class GeneralQuadricSurface: public Object {
    GeneralQuadricSurfaceCoefficient coefficient;
    Vector cubeReferencePoint;
    double length;
    double width;
    double height;

public:
    GeneralQuadricSurface() {
        coefficient.a = coefficient.b = coefficient.c = coefficient.d = coefficient.e = 0.0;
        coefficient.f = coefficient.g = coefficient.h = coefficient.i = coefficient.j = 0.0;
        length = width = height = 0.0;
    }

    GeneralQuadricSurface(GeneralQuadricSurfaceCoefficient coefficient, Vector cubeReferencePoint, double length, double width, double height) {
        this->coefficient = coefficient;
        this->cubeReferencePoint = cubeReferencePoint;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    void draw() {
        /* draw(): implemented from base class Object */
    }

    double intersect(Ray, Color&, int);

    ~GeneralQuadricSurface() {
        /* destructor */
    }
};

double GeneralQuadricSurface::intersect(Ray ray, Color& color, int level) {
    /* finding intersecting tMin */
    /* reference: http://skuld.bmsc.washington.edu/people/merritt/graphics/quadrics.html */
    double a, b, c, tMin, tMax;

    a = coefficient.a*ray.Rd.x*ray.Rd.x+coefficient.b*ray.Rd.y*ray.Rd.y+coefficient.c*ray.Rd.z*ray.Rd.z;
    a += coefficient.d*ray.Rd.x*ray.Rd.y+coefficient.e*ray.Rd.x*ray.Rd.z+coefficient.f*ray.Rd.y*ray.Rd.z;

    b = 2.0*coefficient.a*ray.R0.x*ray.Rd.x+2.0*coefficient.b*ray.R0.y*ray.Rd.y+2.0*coefficient.c*ray.R0.z*ray.Rd.z;
    b += coefficient.d*(ray.R0.x*ray.Rd.y+ray.Rd.x*ray.R0.y);
    b += coefficient.e*(ray.R0.x*ray.Rd.z+ray.Rd.x*ray.R0.z);
    b += coefficient.f*(ray.R0.y*ray.Rd.z+ray.Rd.y*ray.R0.z);
    b += coefficient.g*ray.Rd.x+coefficient.h*ray.Rd.y+coefficient.i*ray.Rd.z;

    c = coefficient.a*ray.R0.x*ray.R0.x+coefficient.b*ray.R0.y*ray.R0.y+coefficient.c*ray.R0.z*ray.R0.z;
    c += coefficient.d*ray.R0.x*ray.R0.y+coefficient.e*ray.R0.x*ray.R0.z+coefficient.f*ray.R0.y*ray.R0.z;
    c += coefficient.g*ray.R0.x+coefficient.h*ray.R0.y+coefficient.i*ray.R0.z+coefficient.j;

    if(a == 0.0) {
        tMin = (b == 0.0)? INF: -c/b;
        tMax = INF;
    } else {
        double discriminant = b*b-4.0*a*c;

        if(discriminant < 0.0) {
            tMin = tMax = INF;
        } else if(discriminant > 0.0) {
            tMax = -b/(2.0*a)+sqrt(discriminant)/(2.0*a);
            tMin = -b/(2.0*a)-sqrt(discriminant)/(2.0*a);
        } else {
            tMin = -b/(2.0*a);
            tMax = INF;
        }
    }

    /* clipping general quadric surface along the dimensions (if necessary) */
    if(tMin < INF) {
        if(tMax < INF) {
            if(tMin > 0.0) {
                Vector intersectionPoint = ray.R0+ray.Rd*tMin;

                if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                    tMin = INF;
                }
            }
            if(tMax > 0.0) {
                Vector intersectionPoint = ray.R0+ray.Rd*tMax;

                if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                    tMax = INF;
                }
            }
            tMin = (tMin>0.0 && tMin<tMax)? tMin: tMax;
        } else {
            if(tMin > 0.0) {
                Vector intersectionPoint = ray.R0+ray.Rd*tMin;

                if((length!=0.0 && (intersectionPoint.x<cubeReferencePoint.x || intersectionPoint.x>cubeReferencePoint.x+length)) || (width!=0.0 && (intersectionPoint.y<cubeReferencePoint.y || intersectionPoint.y>cubeReferencePoint.y+width)) || (height!=0.0 && (intersectionPoint.z<cubeReferencePoint.z || intersectionPoint.z>cubeReferencePoint.z+height))) {
                    tMin = INF;
                }
            }
        }
    }

    if(level == 0) {
        return tMin;
    }

    /* illuminating with Phong Lighting Model */
    Vector intersectionPoint = ray.R0+ray.Rd*tMin;
    Color intersectionPointColor = getColor();

    /* determining unit normal vector at intersection point on appropriate side of general quadric surface */
    double xNormal, yNormal, zNormal;

    xNormal = 2.0*coefficient.a*intersectionPoint.x+coefficient.d*intersectionPoint.y;
    xNormal += coefficient.e*intersectionPoint.z+coefficient.g;

    yNormal = 2.0*coefficient.b*intersectionPoint.y+coefficient.d*intersectionPoint.x;
    yNormal += coefficient.f*intersectionPoint.z+coefficient.h;

    zNormal = 2.0*coefficient.c*intersectionPoint.z+coefficient.e*intersectionPoint.x;
    zNormal += coefficient.f*intersectionPoint.y+coefficient.i;

    Vector normal(xNormal, yNormal, zNormal);
    normal.normalize();
    normal = ((ray.Rd*(-1.0))/normal > 0.0)? normal: normal*(-1.0);

    /* computing ambient light component of reflected ray */
    computeAmbientLightComponent(color, intersectionPointColor);

    /* computing diffuse & specular reflection components of reflected ray */
    for(int i=0; i<lights.size(); i++) {
        Ray incidentRay(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        /* checking if intersection point is in shadow */
        double t, tMinimum=INF;

        for(int j=0; j<objects.size(); j++) {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.R0+incidentRay.Rd*tMinimum;
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectionPoint.computeDistanceBetween(incidentRay.R0)-epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.R0)) {
            /* intersection point is, indeed, in shadow */
            continue;
        }

        /* computing diffuse & specular components of reflected ray */
        computeReflectionComponents(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    /* handling recursive reflection */
    if(level >= levelOfRecursion) {
        return tMin;
    }

    /* incorporating concept of evil epsilon to recursive reflection computation */
    Vector reflectionDirection = ray.Rd-normal*((ray.Rd/normal)*2.0);
    reflectionDirection.normalize();
    Ray reflectedRay(intersectionPoint+reflectionDirection, reflectionDirection);

    /* finding nearest intersecting object (if available) */
    int nearest = INT_MAX;
    double t, tMinimum=INF;

    for(int i=0; i<objects.size(); i++) {
        Color dummyColor;  // color = black
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 && t<tMinimum) {
            tMinimum = t;
            nearest = i;
        }
    }

    /* finding color component for reflected ray */
    Color reflectedColor;  // color = black

    if(nearest != INT_MAX) {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    /* computing recursive reflection component of reflected ray */
    computeRecursiveReflectionComponent(color, reflectedColor);

    /* clipping the color values (if necessary) */
    color.red = (color.red > 1.0)? 1.0: ((color.red < 0.0)? 0.0: color.red);
    color.green = (color.green > 1.0)? 1.0: ((color.green < 0.0)? 0.0: color.green);
    color.blue = (color.blue > 1.0)? 1.0: ((color.blue < 0.0)? 0.0: color.blue);

    return tMin;
}

/* Floor class */
class Floor: public Object {
    double floorWidth;
    double tileWidth;

    /* considering color from base class Object as background color */
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
    double intersect(Ray, Color&, int);

    ~Floor() {
        /* destructor */
    }
};

void Floor::draw() {
    for(int i=0, row=(int) floorWidth/tileWidth, column=(int) floorWidth/tileWidth; i<row; i++) {
        for(int j=0; j<column; j++) {
            /* drawing square on a plane parallel to x-y plane */
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

double Floor::intersect(Ray ray, Color& color, int level) {
    /* determining unit normal vector on appropriate side of floor */
    Vector normal(0.0, 0.0, 1.0);
    normal = (pos / normal > 0.0) ? normal : normal * (-1.0);

    /* finding intersecting tMin */
    double tMin = INF;

    if(normal/ray.Rd != 0.0) {
        tMin = (-1.0)*(normal/ray.R0)/(normal/ray.Rd);
    }

    if(tMin>0.0 && tMin<INF) {
        /*
            ray intersects the floor plane and is in front of the camera,
                but we need to make sure the intersection point is on the floor
        */
        Vector intersectionPoint = ray.R0+ray.Rd*tMin;

        if(!((intersectionPoint.x>-floorWidth/2.0 && intersectionPoint.x<floorWidth/2.0) && (intersectionPoint.y>-floorWidth/2.0 && intersectionPoint.y<floorWidth/2.0))) {
            /* intersection point is not on the floor */
            tMin = INF;
        }
    }

    if(level == 0) {
        return tMin;
    }

    /* illuminating with Phong Lighting Model */
    Vector intersectionPoint = ray.R0+ray.Rd*tMin;
    Vector referencePosition = intersectionPoint-Vector(-floorWidth/2.0, -floorWidth/2.0, 0.0);
    Color intersectionPointColor = (((int) (floor(referencePosition.x/tileWidth)+floor(referencePosition.y/tileWidth)))%2 == 0)? getColor(): foregroundColor;

    /* computing ambient light component of reflected ray */
    computeAmbientLightComponent(color, intersectionPointColor);

    /* computing diffuse & specular reflection components of reflected ray */
    for(int i=0; i<lights.size(); i++) {
        Ray incidentRay(lights[i].getPosition(), intersectionPoint-lights[i].getPosition());

        /* checking if intersection point is in shadow */
        double t, tMinimum=INF;

        for(int j=0; j<objects.size(); j++) {
            Color dummyColor;  // color = black
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if(t>0.0 && t<tMinimum) {
                tMinimum = t;
            }
        }

        Vector shadowIntersectionPoint = incidentRay.R0+incidentRay.Rd*tMinimum;
        double epsilon = 0.0000001;  // for tuning light effect

        if(intersectionPoint.computeDistanceBetween(incidentRay.R0)-epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.R0)) {
            /* intersection point is, indeed, in shadow */
            continue;
        }

        /* computing diffuse & specular components of reflected ray */
        computeReflectionComponents(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    /* handling recursive reflection */
    if(level >= levelOfRecursion) {
        return tMin;
    }

    /* incorporating concept of evil epsilon to recursive reflection computation */
    Vector reflectionDirection = ray.Rd-normal*((ray.Rd/normal)*2.0);
    reflectionDirection.normalize();
    Ray reflectedRay(intersectionPoint+reflectionDirection, reflectionDirection);

    /* finding nearest intersecting object (if available) */
    int nearest = INT_MAX;
    double t, tMinimum=INF;

    for(int i=0; i<objects.size(); i++) {
        Color dummyColor;  // color = black
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if(t>0.0 && t<tMinimum) {
            tMinimum = t;
            nearest = i;
        }
    }

    /* finding color component for reflected ray */
    Color reflectedColor;  // color = black

    if(nearest != INT_MAX) {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level+1);
    }

    /* computing recursive reflection component of reflected ray */
    computeRecursiveReflectionComponent(color, reflectedColor);

    /* clipping the color values (if necessary) */
    color.red = (color.red > 1.0)? 1.0: ((color.red < 0.0)? 0.0: color.red);
    color.green = (color.green > 1.0)? 1.0: ((color.green < 0.0)? 0.0: color.green);
    color.blue = (color.blue > 1.0)? 1.0: ((color.blue < 0.0)? 0.0: color.blue);

    return tMin;
}