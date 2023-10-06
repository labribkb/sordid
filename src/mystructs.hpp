#ifndef LAYOUT_H
#define LAYOUT_H
#include <omp.h>

#include <vector>
#include <string>
#include <tuple>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
//typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
typedef CGAL::Cartesian<double> Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;
typedef CGAL::Delaunay_mesh_vertex_base_2<Kernel>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<Kernel>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,Tds> CDT;


struct Vec2d
{
    double a;
    double b;
    Vec2d(double x, double y) : a(x), b(y) {}
    Vec2d() = default;

    inline void set(double a, double b)
    {
        this->a = a;
        this->b = b;
    }
};

struct Coord : Vec2d
{
    Coord(const Coord &c) : Vec2d(c.x(), c.y()) {}
    Coord(double x, double y) : Vec2d(x, y) {}
    Coord() = default;

    inline double x() const { return a; }
    inline double y() const { return b; }
    inline void set_x(double x) { a = x; }
    inline void set_y(double y) { b = y; }
};

struct Size : Vec2d
{
    Size(const Size &s) : Vec2d(s.width(), s.height()) {}
    Size(double width, double height) : Vec2d(width, height) {}
    Size() = default;

    inline double width() const { return this->a; }
    inline double height() const { return this->b; }
    inline void set_width(double x) { a = x; }
    inline void set_height(double y) { b = y; }
};
// Compute distances
double FORBID_delta(Coord const &ci, Coord const &cj, Size const &si, Size const &sj, double intersec_width, double intersec_height);
double SIDE2SIDE_delta(Coord const &ci, Coord const &cj, Size const &si, Size const &sj, double intersec_width, double intersec_height);
double PRISM_delta(Coord const &ci, Coord const &cj, Size const &si, Size const &sj, double intersec_width, double intersec_height);

struct Parametrizer
{
    double K = 4;
    double ALPHA = 2;
    double MINIMUM_MOVEMENT = 1e-6;
    int MAX_ITER = 30;
    int MAX_PASSES = 100;
    double SCALE_STEP = 0.1;
    double (*distance_fn)(Coord const &, Coord const &, Size const &, Size const &, double, double) = &FORBID_delta;
    double delta = 0.03;
    double eps = 0.01;
    int seed = 0;
    bool PRIME = false;
    int IMG_WIDTH = 250;
    int N = -1;
    bool MONITOR=false;
};

struct term
{
    int i, j;
    double d, w;
    bool o;
    bool f;
    int face_i;
    int face_j;
    term(int i, int j, double d, double w, bool o, bool f,int face_i,int face_j) : i(i), j(j), d(d), w(w), o(o), f(f),face_i(face_i),face_j(face_j) {}
    term(int i, int j, double d, double w, bool o, bool f) : i(i), j(j), d(d), w(w), o(o), f(f),face_i(0),face_j(0) {}
    term(int i, int j, double d, double w, bool o) : i(i), j(j), d(d), w(w), o(o),face_i(0),face_j(0) {}
    term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w),face_i(0),face_j(0) {}
};
class Shape{
protected: 
    Kernel::Vector_2 recenter;
public:
    void set_recenter(const Kernel::Vector_2 c){this->recenter = c;}
    const Kernel::Vector_2 & recenter_vector(){return this->recenter;}
    virtual double bounding_circle_squared_radius()const =0;
    virtual double area()const=0;
    virtual double left_most_x()const=0;
    virtual double right_most_x()const=0;
    virtual double top_most_y()const=0;
    virtual double bottom_most_y()const=0;
    virtual bool do_intersect(const Shape * other, const Coord & thisPos,const Coord & otherPos)const = 0;
    virtual void save(std::ostream &f,const std::string & sep)const = 0;
    virtual Shape * clone()const = 0;
    virtual void terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const;
    virtual double optimalDist(const Shape * other,bool &overlap,double &a,const Coord & thisPos,const Coord & otherPos)const;
    virtual double maxScaleRatio(const Shape* other,const Coord & thisPos, const Coord & otherPos)const;
    virtual void triangulate();
    virtual Coord center();
    virtual ~Shape() {}

    
};

class PolygonShape: public Shape{
    CGAL::Polygon_with_holes_2<Kernel> geometry;
    double r;
    double a;
    double left;
    double right;
    double top;
    double bottom;

    std::vector<CGAL::Triangle_2<Kernel>> faces;
    std::vector<Kernel::Point_2> faceCenters;
    std::vector<double> faceRadius;
    
public:
    virtual Coord center();
    PolygonShape(const CGAL::Polygon_with_holes_2<Kernel> &p);
    const Coord   getFaceCenters(int i)const ;
    virtual double bounding_circle_squared_radius()const {return this->r;}
    virtual double area()const{return this->a;};
    virtual double left_most_x()const{return this->left;}
    virtual double right_most_x()const{return this->right;}
    virtual double top_most_y()const{return this->top;}
    virtual double bottom_most_y()const{return this->bottom;}
    PolygonShape() = default;
    CGAL::Aff_transformation_2<Kernel> translateShape(const Coord & c) const;
    virtual void save(std::ostream &f,const std::string & sep)const;
    virtual bool do_intersect(const Shape * other, const Coord & thisPos,const Coord & otherPos)const;
    virtual Shape * clone()const;
    PolygonShape(const PolygonShape*other);
    virtual void terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const;
    virtual void triangulate();
    virtual ~PolygonShape() {}

};
class Rectangle : public Shape{
    Size s;
public:
    const Size & size()const{return this->s;}
    virtual double bounding_circle_squared_radius()const {return 0.25*(this->s.width()*this->s.width()+this->s.height()*this->s.height());}
    virtual double area()const{return this->s.width()*this->s.height();};
    virtual double left_most_x()const{return -this->s.width()/2;}
    virtual double right_most_x()const{return this->s.width()/2;}
    virtual double top_most_y()const{return this->s.height()/2;}
    virtual double bottom_most_y()const{return -this->s.height()/2;}
    Rectangle() = default;
    Rectangle(double w,double h):s(w,h){}
    virtual void save(std::ostream &f,const std::string & sep)const;
    virtual bool do_intersect(const Shape * other, const Coord & thisPos,const Coord & otherPos)const;
    virtual Shape * clone()const;
    Rectangle(const Rectangle*other);
    virtual void terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const;
    virtual double maxScaleRatio(const Shape* other,const Coord & thisPos, const Coord & otherPos)const;
    virtual ~Rectangle() {}

};


class INode
{
public:
    Coord coord;
    Size size;
    //CGAL::Polygon_2<Kernel> shape;
    //CGAL::Polygon_with_holes_2<Kernel> shape;
    Shape *shape;
    INode *parent;
    virtual ~INode() {}
    virtual void move(double mvt_x, double mvt_y) = 0;
    virtual void print() = 0;
    virtual void print(std::string prefix) = 0;
    virtual int N() { return -1; }
    virtual int depth() =0;
    virtual int N(bool deep) {return 1;}
    virtual std::tuple<double, double, double, double> getBB() = 0;
    virtual bool overlaps(const INode * other) const;
    virtual double optimalDist(const INode * other,bool & overlap,double& area) const;
    virtual double optimalDist(const INode * other,bool & overlap) const;
    virtual double optimalDist(const INode * other) const;
    virtual void terms(const INode* other,const Parametrizer & params,int i,int j,std::vector<term> & terms)const;
    virtual double maxScaleRatio(const INode* other)const;
};

class Node : public INode
{
private:
    int id;

public:
    Node(const Coord &c, const CGAL::Polygon_with_holes_2<Kernel>&s, int id);
    Node(const Coord &c, const Shape*s, int id);
    Node() = default;
    ~Node();
    void print();
    void print(std::string prefix);
    void move(double mvt_x, double mvt_y);
    int depth() {return 1;}
    std::tuple<double, double, double, double> getBB();
};

class NodesComposite : public INode
{
public:
    std::vector<INode *> nodes;
    NodesComposite(const Coord &c, const Size &s);
    NodesComposite() = default;
    NodesComposite(std::vector<INode *> &nodes);
    NodesComposite(const NodesComposite &other);
    void move(double mvt_x, double mvt_y);
    virtual void add(INode *n);
    void scale(double scaleFactor);
    int N() { return nodes.size(); }
    int N(bool deep);
    void print();
    void print(int n);
    int depth(){return this->N() > 0 ? 1+this->nodes[0]->depth() : 1;};
    std::tuple<double, double, double, double> getBB();
    void print(std::string prefix);
    INode *remove(int i);
};



class Layout : public NodesComposite
{
public:
    Layout(std::vector<Coord> &pos, std::vector<CGAL::Polygon_with_holes_2<Kernel>> &shapes);
    Layout(const Layout &other);
    Layout() = default;
    ~Layout();
    void deepcopy(const Layout &other);
    Layout & operator=(const Layout &other);
    void loadLayout(std::string input_path);
    double maxScaleRatio();
    void move(double mvt_x, double mvt_y){}; // interdire
    bool isCurrentScaleSolvable();
    std::string save(std::string path, std::string extension, double elapsed, bool isPrime);
    void do_triangulate();
};

double maxOfLeft(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double minOfRight(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double minOfTop(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double maxOfBot(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
std::tuple<double, double> nodeRectanglesIntersection(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);


bool contains(Coord &p, Coord &p2, Size &s2);
double vecNorm2D(double vec_x, double vec_y);

double eucl(double x1, double y1, double x2, double y2);
double eucl(Coord const &p1, Coord const &p2);

// Getting overlaps
std::vector<size_t> sort_indexes(const double *v, int size);
bool overlapCheck(INode *ni, INode *nj);
bool scanLineOverlapCheck(NodesComposite &nodes_group);
bool scanLineOverlapCheck(std::vector<INode *> &nodes_group);
void getAllOverlaps(NodesComposite &nodes_group, std::vector<std::tuple<int, int>> &overlaps);



#endif
