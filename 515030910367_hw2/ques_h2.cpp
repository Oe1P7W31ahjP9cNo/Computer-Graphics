#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#define M_PI 3.141592653589793
#define INFINITY 1e8
template<typename T>
//construct the point object
class Vector3
{
public:
    T x, y, z;
    Vector3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vector3(T xx) : x(xx), y(xx), z(xx) {}
    Vector3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vector3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vector3<T> operator * (const T &f) const
    {
        return Vector3<T>(x * f, y * f, z * f);
    }
    Vector3<T> operator * (const Vector3<T> &v) const
    {
        return Vector3<T>(x * v.x, y * v.y, z * v.z);
    }
    T dot(const Vector3<T> &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    Vector3<T> operator - (const Vector3<T> &v) const
    {
        return Vector3<T>(x - v.x, y - v.y, z - v.z);
    }
    Vector3<T> operator + (const Vector3<T> &v) const
    {
        return Vector3<T>(x + v.x, y + v.y, z + v.z);
    }
    Vector3<T>& operator += (const Vector3<T> &v)
    {
        x += v.x, y += v.y, z += v.z; return *this;
    }
    Vector3<T>& operator *= (const Vector3<T> &v)
    {
        x *= v.x, y *= v.y, z *= v.z; return *this;
    }
    Vector3<T> operator - () const
    {
        return Vector3<T>(-x, -y, -z);
    }
    T length2() const
    {
        return x * x + y * y + z * z;
    }
    T length() const
    {
        return sqrt(length2());
    }
};

typedef Vector3<float> Vector3f;
//construct the sphere object
class Sphere
{
public:
    Vector3f center;
    float radius, radius2;
    Vector3f surface_color, emission_color;      // surface color and emission
    float transparency, reflection;      // surface transparency and reflectivity
        // position, radius, surface color, reflectivity, transparency, emission color
    Sphere(
        const Vector3f &c,
        const float &r,
        const Vector3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vector3f &ec = 0) :
        center(c), radius(r), radius2(r * r), surface_color(sc), emission_color(ec),
        transparency(transp), reflection(refl){}

    // Compute a ray-sphere intersection using the geometric solution
    bool intersect(const Vector3f &rayorig, const Vector3f &raydir, float &t0, float &t1) const
    {
        float tca, thc,d2;
        Vector3f l;
        l = center - rayorig;
        tca = l.dot(raydir);
        if (tca < 0) return false;
        d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
};

#define MAX_RAY_DEPTH 6  //controls the maximum recursion depth

//decide the ratio of reflection
float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

Vector3f trace(const Vector3f &rayorig,const Vector3f &raydir,const std::vector<Sphere> &spheres,const int &depth)
{
    float tnear = INFINITY;
    const Sphere* sphere = NULL;
    // find intersection of this ray with the sphere in the scene
    for (int i = 0; i < spheres.size(); ++i)
    {
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayorig, raydir, t0, t1))
        {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    if (!sphere) return Vector3f(2);
    Vector3f surface_color = 0;
    Vector3f phit = rayorig + raydir * tnear; // point of intersection
    Vector3f nhit = phit - sphere->center; // normal at the intersection point
    nhit.normalize();
    float bias = 1e-4;
    bool inside = false;
    if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH)
    {
        float facingratio = -raydir.dot(nhit);
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        // compute reflection direction
        Vector3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        refldir.normalize();
        Vector3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
        Vector3f refraction = 0;
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->transparency)
        {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior;
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vector3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }
        surface_color = (reflection * fresneleffect +refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surface_color;
    }
    else {
        // it's a diffuse object, no need to raytrace any further
        for (int i = 0; i < spheres.size(); ++i)
        {
            if (spheres[i].emission_color.x > 0)
            {
                Vector3f transmission = 1;
                Vector3f lightDirection = spheres[i].center - phit;
                lightDirection.normalize();
                for (int j = 0; j < spheres.size(); ++j)
                {
                    if (i != j)
                    {
                        float t0, t1;
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1))
                        {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surface_color += sphere->surface_color * transmission *
                std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emission_color;
            }
        }
    }

    return surface_color + sphere->emission_color;
}

void render(const std::vector<Sphere> &spheres)
{
    int width = 640, height = 640;
    Vector3f *image = new Vector3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x, ++pixel)
        {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vector3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vector3f(0), raydir, spheres, 0);
        }
    }
    std::ofstream ofs("./raytracing.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; ++i)
    {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
}

int main()
{
    std::vector<Sphere> spheres;
    //cerate the six sphere
    spheres.push_back(Sphere(Vector3f( 0.0, -10004, -20), 10000, Vector3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres.push_back(Sphere(Vector3f( 0.0,      -1, -20),     3, Vector3f(1.00, 0.92, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vector3f( 5.0,     -1, -15),     2, Vector3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vector3f( 5.0,      0, -25),     1, Vector3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vector3f(-5.5,      0, -15),     3, Vector3f(0.90, 0.50, 0.90), 1, 0.0));
    spheres.push_back(Sphere(Vector3f( 0.0,      6, -20),     3, Vector3f(1.00, 0.22, 0.96), 1, 0.5));
    //  create the light
    spheres.push_back(Sphere(Vector3f( 0.0,     20, -30),     3, Vector3f(0.00, 0.00, 0.00), 0, 0.0, Vector3f(3)));
    render(spheres);
    return 0;
}
