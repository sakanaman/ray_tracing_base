#ifndef BDPT_HPP
#define BDPT_HPP

#include "material.hpp"
#include "IBL.hpp"
#include <cmath>
#include <concurrent.hpp>
#include <wrap_embree.hpp>

class Ray
{
public:
    Ray(){}
    Ray(const Vec3<float>& origin, const Vec3<float>& dir)
    :origin(origin),dir(dir){}
    Ray(const float* _origin, const float* _dir)
    {
        origin = Vec3<float>(_origin);
        dir = Vec3<float>(_dir);
    }

    Vec3<float> origin;
    Vec3<float> dir;
};

//=============================
//       BDPT METHODS
//=============================

// from pbrt
inline uint32_t FloatToBits(float f) {
    uint32_t ui;
    std::memcpy(&ui, &f, sizeof(float));
    return ui;
}

inline float BitsToFloat(uint32_t ui) {
    float f;
    std::memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}

// from pbrt
// Parallel Declarations(I wanna use atomic<float>...)
class AtomicFloat {
public:
    // AtomicFloat Public Methods
    explicit AtomicFloat(float v = 0) { bits = FloatToBits(v); }

    operator float() const { return BitsToFloat(bits); }

    float operator=(float v) {
        bits = FloatToBits(v);
        return v;
    }

    void Add(float v) 
    {
        uint32_t oldBits = bits, newBits;
        do {
            newBits = FloatToBits(BitsToFloat(oldBits) + v);
        } while (!bits.compare_exchange_weak(oldBits, newBits));
    }

private:
    std::atomic<uint32_t> bits;
};

// for t = 1 contribution
class LightImage
{
public:
    LightImage(int width, int height) : width(width), height(height)
    {
        film.resize(width*height, AtomicFloat(0.0f));
    };
    ~LightImage();
    void add(const int i, const int j, const Vec3<float>& color);
private:
    int width, height;
    std::vector<AtomicFloat> film;
};

enum class VertexType
{
    RAMP,
    SURFACE,
    CAMERA,
    UNDEFINED
};

class Vertex
{
public:
    Vertex();
    ~Vertex();
    float getPDFfwd() const
    {
        return pdfFwd;
    }
    float getPDFrev() const
    {
        return pdfRev;
    }
    VertexType getVertexType() const
    {
        return vtype;
    }
    void set_dir_next(const Vec3<float>& _w_next){w_next = _w_next;}
    void set_dir_prev(const Vec3<float>& _w_prev){w_prev = _w_prev;}
protected:
    Vec3<float> w_prev = {0.0f, 0.0f, 0.0f};
    Vec3<float> w_next = {0.0f, 0.0f, 0.0f};
    bool is_delta = false;
    embree::IsectInfo<float> isect;
    float pdfFwd = -1.0f;
    float pdfRev = -1.0f;
    VertexType vtype = VertexType::UNDEFINED;
};

class EyeVertex : public Vertex
{
public:
    EyeVertex();
    ~EyeVertex();
    void createVertex();
};

class LightVertex : public Vertex
{
public:
    LightVertex();
    ~LightVertex();
    void createVertex();
};

class Subpath
{
public:
    Subpath(int maxdepth):maxdepth(maxdepth){};
    ~Subpath();
protected:
    int maxdepth = 7;
};

class LightSubpath : public Subpath
{
public:
    LightSubpath();
    ~LightSubpath();
    // create y(0) vertex
    void init();
    // create y(1), ... y(s-1) vertex
    void build( const Shader<float>& shader,
                const embree::EmbreeManager& emb, RandomManager& rnd_manager, 
                const SceneData<float>& scenedata);
private:
    std::vector<LightVertex> lvs;
    Vec3<float> alpha_l= {1.0f, 1.0f, 1.0f};
};

class EyeSubpath : public Subpath
{
public:
    EyeSubpath();
    ~EyeSubpath();
    // create z(0) vertex
    void init();
    // create z(1), ... z(t-1) vertex
    void build( const Shader<float>& shader,
                const embree::EmbreeManager& emb, RandomManager& rnd_manager, 
                const SceneData<float>& scenedata);
private:
    std::vector<EyeVertex> evs;
    Vec3<float> alpha_e = {1.0f, 1.0f, 1.0f};
};

// calculate C*(s,t)
void calcContrib_unweighted(const int s, const int t, const EyeSubpath& eyepath,
                            const LightSubpath& lightpath,
                            Vec3<float>* contrib);

// calculate W(s,t)
void calcMIS(const int s, const int t, const EyeSubpath& eyepath,
             const LightSubpath& lightpath,
             float* weight);

// task per pixel
Vec3<float> BPT_Test(const float* firstRay_dir, const float* firstRay_origin, 
                     const Shader<float>& shader, LightImage* lightimage,
                     const embree::EmbreeManager& emb, RandomManager& rnd_manager, 
                     const SceneData<float>& scenedata);


#endif