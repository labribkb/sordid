#include <iostream>
#include <fstream>
#include <float.h>
#include <CGAL/Bbox_2.h>
//#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/intersections.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/partition_2.h>
#include <CGAL/Polygon_2.h>



#include "mystructs.hpp"

using namespace std;

Node::Node(const Coord &c, const CGAL::Polygon_with_holes_2<Kernel> &s, int id)
{
    this->coord = Coord(c);
    this->shape = new PolygonShape(s);
    CGAL::Bbox_2 bb = s.bbox();
    this->size = Size(bb.xmax()-bb.xmin(),bb.ymax()-bb.ymin());
    this->id = id;
}
Node::Node(const Coord &c, const Shape*s, int id)
{
    this->coord = Coord(c);
    this->shape=s->clone();
    //CGAL::Bbox_2 bb = s.geometry.bbox();
    this->size = Size(s->right_most_x()-s->left_most_x(),s->top_most_y()-s->bottom_most_y());
    this->id = id;
}
Node::~Node(){
    delete this->shape;
}

NodesComposite::NodesComposite(vector<INode *> &nodes)
{
    this->nodes = nodes;
}

NodesComposite::NodesComposite(const Coord &c, const Size &s)
{
    this->coord = Coord(c);
    this->size = Size(s);
}

void Node::move(double mvt_x, double mvt_y)
{
    this->coord.set_x(this->coord.x() + mvt_x);
    this->coord.set_y(this->coord.y() + mvt_y);
}

void NodesComposite::move(double mvt_x, double mvt_y)
{
    this->coord.set_x(this->coord.x() + mvt_x);
    this->coord.set_y(this->coord.y() + mvt_y);
    for (int i = 0; i < this->nodes.size(); ++i)
        this->nodes[i]->move(mvt_x, mvt_y);
}

NodesComposite::NodesComposite(const NodesComposite &other)
{
    this->coord = Coord(other.coord);
    this->size = Size(other.size);
    this->nodes = vector<INode *>(other.nodes);
}

Layout::~Layout()
{
    // cout << "LAYOUT DESTROYER" << endl;
    for (int i = 0; i < this->N(); ++i)
        delete this->nodes[i];
}

INode *NodesComposite::remove(int i)
{
    INode *ret;
    if (i > 0 && i < this->nodes.size())
    {
        ret = this->nodes.at(i);
        this->nodes.erase(this->nodes.begin() + i);
    }
    return ret;
}

Layout::Layout(vector<Coord> &pos, vector<CGAL::Polygon_with_holes_2<Kernel>> &shapes)
{
    for (int i = 0; i < pos.size(); ++i)
        this->nodes.push_back(new Node(pos[i], shapes[i], i));
}

Layout::Layout(const Layout &other)
{
    this->deepcopy(other);
}
void Layout::deepcopy(const Layout &other)
{
    this->nodes.clear();
    this->nodes.reserve(other.nodes.size());
    for (int i = 0; i < other.nodes.size(); ++i)
    {
        this->nodes.push_back(new Node(Coord(other.nodes[i]->coord.x(), other.nodes[i]->coord.y()), other.nodes[i]->shape, i));
    }
}
Layout &Layout::operator=(const Layout &other)
{
    this->deepcopy(other);
    return *this;
}
void NodesComposite::scale(double scaleFactor)
{
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        this->nodes[i]->coord.set(this->nodes[i]->coord.x() * scaleFactor, this->nodes[i]->coord.y() * scaleFactor);
    }
}

void Node::print()
{
    cout << "ID = " << id << "\tX: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t"
         << "W: " << this->size.width() << " ; H: " << this->size.height() << endl;
}

void Node::print(string prefix)
{
    cout << prefix;
    this->print();
}

void NodesComposite::print()
{
    cout << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;
    for (int i = 0; i < this->N(); i++)
    {
        cout << "\t i: " << i << "\t\t";
        this->nodes[i]->print();
    }
    cout << "\n"
         << endl;
}

void NodesComposite::print(string prefix)
{
    cout << prefix << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;
    for (int i = 0; i < this->N(); i++)
    {
        cout << prefix << "i: " << i << "\t";
        this->nodes[i]->print(prefix + "\t\t");
    }
}

void NodesComposite::print(int n)
{
    cout << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;

    int bound = min(n, this->N());
    for (int i = 0; i < bound; i++)
    {
        cout << "\t i: " << i << "\t\t";
        this->nodes[i]->print();
    }
}

// left, right, bot, top
tuple<double, double, double, double> Node::getBB()
{
    return make_tuple(
        this->coord.x() - this->size.width() / 2,
        this->coord.x() + this->size.width() / 2,
        this->coord.y() - this->size.height() / 2,
        this->coord.y() + this->size.height() / 2);
}

// left, right, bot, top
tuple<double, double, double, double> NodesComposite::getBB()
{
    double leftmost = DBL_MAX;
    double rightmost = -DBL_MAX;
    double topmost = -DBL_MAX;
    double botmost = DBL_MAX;

    for (int i = 0; i < this->N(); ++i)
    {
        auto [left, right, bot, top] = this->nodes[i]->getBB();
        leftmost = min(leftmost, left);
        rightmost = max(rightmost, right);
        botmost = min(botmost, bot);
        topmost = max(topmost, top);
    }
    return make_tuple(leftmost, rightmost, botmost, topmost);
}


int NodesComposite::N(bool deep)
{
    if (!deep)
        return this->N();

    int n = 0;
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        n += this->nodes[i]->N(true);
    }
    return n;
}

double Layout::maxScaleRatio()
{
    double padding = 1e-4;
    double maxRatio = 1.;
    double optimalDist, actualDist, ratio, unoverlapRatio;
    double actualX, actualY, desiredWidth, desiredHeight, widthRatio, heightRatio;
    vector<tuple<int, int>> overlaps;
    getAllOverlaps(*this, overlaps);
    INode *nu, *nv;
    for (unsigned int i = 0; i < overlaps.size(); ++i)
    {
        auto [u, v] = overlaps[i];
        nu = this->nodes[u];
        nv = this->nodes[v];
        ratio = nu->maxScaleRatio(nv);
        maxRatio = max(maxRatio, ratio);
    }
   
    return maxRatio;
}

bool Layout::isCurrentScaleSolvable()
{
    double areas_sum = 0;
    double min_x = DBL_MAX;
    double min_y = DBL_MAX;
    double max_x = -DBL_MAX;
    double max_y = -DBL_MAX;

    double left, right, top, bot;
    INode *n;
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        n = this->nodes[i];


	areas_sum+=n->shape->area();

        left = n->coord.x() + n->shape->left_most_x();
        right = n->coord.x() + n->shape->right_most_x();
        top = n->coord.y() + n->shape->top_most_y();
        bot = n->coord.y() + n->shape->bottom_most_y();
        if (left < min_x)
            min_x = left;
        if (right > max_x)
            max_x = right;
        if (bot < min_y)
            min_y = bot;
        if (top > max_y)
            max_y = top;
    }
    double bb_area = (max_x - min_x) * (max_y - min_y);

    return bb_area >= areas_sum;
}

string Layout::save(string path, string extension, double elapsed, bool isPrime)
{
    string sep = " ";
    if (isPrime)
        extension += "p";
    ofstream f(path + "." + extension);
    f << elapsed << endl;
    INode *n;
    for (unsigned int i = 0; i < this->nodes.size(); ++i)
    {
        n = this->nodes[i];
        f << n->coord.x()-n->shape->recenter_vector().x() << sep << n->coord.y()-n->shape->recenter_vector().y() << sep;
        n->shape->save(f,sep);

         f << endl;
    }
    f.close();
    cout << "save : " << path + "." + extension << endl;
    return path + "." + extension;
}

void Layout::loadLayout(string input_path)
{
    ifstream f;
    f.open(input_path);
    int N;
    f >> N;

    vector<Coord> pos;
    pos.reserve(N);
    vector<Shape*> shapes;
    shapes.reserve(N);

    double x, y, w, h;
    for (int i = 0; i < N; i++)
    {
        f >> x >> y;
	size_t nb_vertices;
	f>>nb_vertices;
        if(nb_vertices==0){
            
        }else
        {
	        double xv,yv;
            CGAL::Polygon_2<Kernel> outer;
	        for(size_t j=0;j<nb_vertices;j++){
	            f >> xv >> yv;
	            outer.push_back(Kernel::Point_2(xv,yv));
	        }
            shapes.push_back(new PolygonShape(CGAL::Polygon_with_holes_2<Kernel>(outer)));
        }
        
        pos.emplace_back(Coord(x, y));
    }
    f.close();

    this->nodes.clear();
    this->nodes.resize(pos.size());
    for (int i = 0; i < pos.size(); ++i){
        this->nodes[i] = new Node(pos[i], shapes[i], i);
        delete shapes[i];
    }
}



void NodesComposite::add(INode *n)
{
    this->nodes.push_back(n);
}

double INode::optimalDist(const INode * other) const{
    bool b;
    return this->optimalDist(other,b);
}
double INode::optimalDist(const INode * other,bool & overlap) const{
    double ign;
    return this->optimalDist(other,overlap,ign);
}
PolygonShape::PolygonShape(const CGAL::Polygon_with_holes_2<Kernel> &p){
    this->geometry = p;
    this->a = std::abs(CGAL::to_double(p.outer_boundary().area()));
}
Coord PolygonShape::center(){
    CGAL::Polygon_2<Kernel> outer= this->geometry.outer_boundary();
    double xv,yv;
    CGAL::Min_sphere_of_spheres_d<CGAL::Min_sphere_of_points_d_traits_2<Kernel,Kernel::FT>> mc( outer.begin(), outer.end());
    CGAL::Min_sphere_of_spheres_d<CGAL::Min_sphere_of_points_d_traits_2<Kernel,Kernel::FT>>::Cartesian_const_iterator itc = mc.center_cartesian_begin();

    xv=*itc;
    itc++;
    yv=*itc;
    CGAL::Vector_2<Kernel> v(xv,yv);
    for(auto it =outer.begin();it!=outer.end();++it)
        (*it)-=v;
    this->geometry = CGAL::Polygon_with_holes_2<Kernel>(outer);
    set_recenter(v);


    CGAL::Polygon_2<Kernel> out= this->geometry.outer_boundary();
    this->r = 0.0;
    Kernel::FT r=0.0;
    for (auto it = out.vertices_begin();it!=out.vertices_end();++it){
        r = std::max(r,CGAL::Segment_2<Kernel>(CGAL::Origin(),*it).squared_length());
    }
    this->r = r;
    this->left= out.left_vertex()->x();
    this->right= out.right_vertex()->x();
    this->top= out.top_vertex()->y();
    this->bottom= out.bottom_vertex()->y();
    this->a = std::abs(out.area());
    return Coord(v.x(),v.y());
}
Coord Shape::center(){
    return Coord(0,0);
}
void Shape::triangulate(){
}
void Layout::do_triangulate(){
#pragma omp parallel for 
    for(int i=0;i<this->nodes.size();++i){
        Coord c = this->nodes[i]->shape->center();
        this->nodes[i]->coord.set_x(this->nodes[i]->coord.x()+c.x());
        this->nodes[i]->coord.set_y(this->nodes[i]->coord.y()+c.y());
        this->nodes[i]->shape->triangulate();
    }
}
void PolygonShape::triangulate(){
    const CGAL::Polygon_2<Kernel> & out= this->geometry.outer_boundary();
    CDT cdt;
    cdt.insert_constraint(out.vertices_begin(), out.vertices_end(),true);
    std::vector<Kernel::Point_2> seeds;
    seeds.push_back(Kernel::Point_2(0,0));
    for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end(); ++fit){
        double baryy=0.0;
        double baryx=0.0;
        for(int i=0;i<3;++i){
            baryx+=fit->vertex(i)->point().x();
            baryy+=fit->vertex(i)->point().y();
        }
        baryx/=3;
        baryy/=3;
        
        if(out.bounded_side(Kernel::Point_2(baryx,baryy))==CGAL::ON_BOUNDED_SIDE) 
        {
            
            CGAL::Min_sphere_of_spheres_d<CGAL::Min_sphere_of_points_d_traits_2<Kernel,Kernel::FT>> mc;
            mc.insert(fit->vertex(0)->point());
            mc.insert(fit->vertex(1)->point());
            mc.insert(fit->vertex(2)->point());
            CGAL::Min_sphere_of_spheres_d<CGAL::Min_sphere_of_points_d_traits_2<Kernel,Kernel::FT>>::Cartesian_const_iterator itc = mc.center_cartesian_begin();
            double xv=*itc;
            itc++;
            double yv=*itc;
            Kernel::Point_2 center = Kernel::Point_2(xv,yv);
            this->faceCenters.push_back(center);
            this->faceRadius.push_back(max(
                CGAL::sqrt(Kernel::Segment_2(fit->vertex(0)->point(),center).squared_length()),
                max(		
                CGAL::sqrt(Kernel::Segment_2(fit->vertex(1)->point(),center).squared_length()),
                CGAL::sqrt(Kernel::Segment_2(fit->vertex(2)->point(),center).squared_length())
                )
            ));

            this->faces.emplace_back(fit->vertex(0)->point(),fit->vertex(1)->point(),fit->vertex(2)->point());
        }
    }
}

double INode::optimalDist(const INode * other,bool & overlap,double & overlapArea) const{
    double r1 = (this->shape->bounding_circle_squared_radius());
    double r2 = (other->shape->bounding_circle_squared_radius());
    double sr1 = sqrt(r1);
    double sr2 = sqrt(r2);
    double d = eucl(this->coord,other->coord);
    if(sr1+sr2<=d){
        overlap=false;
        overlapArea=0.0;
        return d;
    }
    if(this->shape->do_intersect(other->shape,this->coord,other->coord))
    {
        double r = sr1+sr2;
        overlap = r>d;
        if(overlap){
                return r;
        }else{
            overlapArea=0.0;
            return d;
        }

    }
    overlap=false;
    overlapArea=0.0;
    return d;
}
bool INode::overlaps(const INode * other) const{
    if(sqrt(this->shape->bounding_circle_squared_radius())+sqrt(other->shape->bounding_circle_squared_radius())<=eucl(this->coord,other->coord)){
        return false;
    }
    return this->shape->do_intersect(other->shape,this->coord,other->coord);
}

const Coord PolygonShape::getFaceCenters(int i)const {
    return Coord(this->faceCenters[i].x(),this->faceCenters[i].y());
}
bool PolygonShape::do_intersect(const Shape * other, const Coord & thisPos,const Coord & otherPos)const{
    const PolygonShape * otherShape = dynamic_cast<const PolygonShape *>(other);
    for(int i=0;i<this->faces.size();++i){
        double ri = this->faceRadius[i];
        auto thisTranslate = this->translateShape(thisPos);
        // const Kernel::Point_2 & ci = this->faceCenters[i].transform(thisTranslate);
        Kernel::Point_2 ci(this->faceCenters[i]);
        ci = ci.transform(thisTranslate);
        if(nullptr!=otherShape){
            for(int j=0;j<otherShape->faces.size();++j){
                double rj = otherShape->faceRadius[j];
                auto otherTranslate = otherShape->translateShape(otherPos);
                // const Kernel::Point_2 & cj = otherShape->faceCenters[j].transform(otherTranslate);
                Kernel::Point_2 cj(otherShape->faceCenters[j]);
                cj = cj.transform(otherTranslate);
                if(ri+rj>=CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(ci,cj)))) {
                    CGAL::Triangle_2<Kernel> thisPolygon = this->faces[i].transform(thisTranslate);
                    CGAL::Triangle_2<Kernel> otherPolygon = otherShape->faces[j].transform(otherTranslate);
                    if(CGAL::do_intersect(thisPolygon,otherPolygon))
                        return true;
                }
            }
        }
    }
    return false;
}
CGAL::Aff_transformation_2<Kernel> PolygonShape::translateShape(const Coord & c) const{
    return CGAL::Aff_transformation_2<Kernel>(CGAL::TRANSLATION,Kernel::Vector_2(c.x(),c.y()));
}
void PolygonShape::save(ostream & f,const string & sep)const{
    const CGAL::Polygon_2<Kernel> & p =this->geometry.outer_boundary();
    f<<p.size();
    for (auto it=p.vertices_begin();it!=p.vertices_end();++it)
        f<<sep<<it->x()+this->recenter.x()<<sep<<it->y()+this->recenter.y();
  
}
Shape * PolygonShape::clone()const{
    return new PolygonShape(this);
}
PolygonShape::PolygonShape(const PolygonShape*other){
    this->geometry = other->geometry;
    this->a = other->a;
    this->r = other->r;
    this->left = other->left;
    this->right = other->right;
    this->top = other->top;
    this->bottom = other->bottom;
    this->faces = other->faces;
    this->faceCenters = other->faceCenters;
    this->faceRadius = other->faceRadius;
    this->recenter = other->recenter;
}
void INode::terms(const INode* other,const Parametrizer & params,int i,int j,vector<term> &terms)const{
    this->shape->terms(other->shape,params,i,j,terms,this->coord,other->coord);
}
void Shape::terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const{
    bool overlap;
    double a;
    double d_ij = this->optimalDist(other,overlap,a,thisPos,otherPos);
    if(overlap)
        terms.push_back(term(i, j, d_ij, pow(d_ij, params.K * params.ALPHA), overlap));
    else
        terms.push_back(term(i, j, d_ij, pow(d_ij, params.ALPHA), overlap));
}
double Shape::optimalDist(const Shape * other,bool &overlap,double &a,const Coord & thisPos,const Coord & otherPos)const {
    if(this->do_intersect(other,thisPos,otherPos)){
        overlap=true;
        return sqrt(this->bounding_circle_squared_radius())+sqrt(other->bounding_circle_squared_radius());
    }else{
        overlap=false;
        return eucl(thisPos,otherPos);
    }
}
void PolygonShape::terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const{
    bool foundOverlap = false;
    vector<term> localTerm;
    if(sqrt(this->bounding_circle_squared_radius())+sqrt(other->bounding_circle_squared_radius())>eucl(thisPos,otherPos)){
        const PolygonShape * otherShape = dynamic_cast<const PolygonShape *>(other);
        for(int fi=0;fi<this->faces.size();++fi){
            double ri = this->faceRadius[fi];
            Kernel::Point_2 ci(this->faceCenters[fi]);
            ci+=CGAL::Vector_2<Kernel>(thisPos.x(),thisPos.y());
            if(nullptr==otherShape){
                
            }else{
                for(int fj=0;fj<otherShape->faces.size();++fj){
                    double rj = otherShape->faceRadius[fj];
                    Kernel::Point_2 cj(otherShape->faceCenters[fj]);
                    cj+=CGAL::Vector_2<Kernel>(otherPos.x(),otherPos.y());

                    auto cij = CGAL::sqrt(CGAL::squared_distance(ci,cj));
                    if(ri+rj>cij){
                        auto thisTranslate = this->translateShape(thisPos);
                        auto otherTranslate = otherShape->translateShape(otherPos);
                        const CGAL::Triangle_2<Kernel> & thisPolygon = this->faces[fi].transform(thisTranslate);
                        const CGAL::Triangle_2<Kernel> & otherPolygon = otherShape->faces[fj].transform(otherTranslate);
                        if(CGAL::do_intersect(thisPolygon,otherPolygon)){
                            double d_ij = (ri+rj);
                            double w_ij = pow(d_ij, params.K * params.ALPHA);
                            localTerm.push_back(term(i, j,d_ij,  w_ij, true,false,fi+1,fj+1));
                            foundOverlap = true;
                        }else{
                            
                            double d_ij = cij;
                            double w_ij = pow(d_ij, params.ALPHA);
                            localTerm.push_back(term(i, j,d_ij,  w_ij, false,false,fi+1,fj+1));
                        }

                    }else{
                        
                        double d_ij = cij;
                        double w_ij = pow(d_ij, params.ALPHA);
                        localTerm.push_back(term(i, j,d_ij,  w_ij, false,false,fi+1,fj+1));
                    }
                }
            }
        }
    }
    if(!foundOverlap){
        double d_ij = eucl(thisPos,otherPos);
        terms.push_back(term(i, j, d_ij, pow(d_ij, params.ALPHA), false));
    }else{

        terms.insert(terms.end(),localTerm.begin(),localTerm.end());

    }

}
double INode::maxScaleRatio(const INode* other)const{
    return this->shape->maxScaleRatio(other->shape,this->coord,other->coord);
}
double Shape::maxScaleRatio(const Shape* other,const Coord & thisPos, const Coord & otherPos)const{
    bool overlap;
    double a;
    return this->optimalDist(other,overlap,a,thisPos,otherPos)/eucl(thisPos,otherPos);
}

void Rectangle::save(std::ostream &f,const std::string & sep)const{
    f<<this->s.width()<<sep<<this->s.height();
}
bool Rectangle::do_intersect(const Shape * other, const Coord & thisPos,const Coord & otherPos)const{
    const Rectangle * otherRect = dynamic_cast<const Rectangle *>(other);
    tuple<double, double> inter = nodeRectanglesIntersection(thisPos, this->s, otherPos, otherRect->s);
    return get<0>(inter) > 0 && get<1>(inter) > 0;
}
Shape * Rectangle::clone()const{
    return new Rectangle(this->s.width(),this->s.height());
}
Rectangle::Rectangle(const Rectangle*other):s(other->s){}
void Rectangle::terms(const Shape* other,const Parametrizer & params,int i,int j,std::vector<term> & terms,const Coord & thisPos, const Coord & otherPos)const{
    const Rectangle * otherRect = dynamic_cast<const Rectangle *>(other);
    tuple<double, double> inter = nodeRectanglesIntersection(thisPos, this->s, otherPos, otherRect->s);
    if(get<0>(inter) > 0 && get<1>(inter) > 0){
        double d_ij = (*params.distance_fn)(thisPos, otherPos, this->s, otherRect->s, get<0>(inter), get<1>(inter));
        double w_ij = pow(d_ij, params.K * params.ALPHA);
        terms.push_back(term(i, j, d_ij, w_ij, true));
    }
    else
    {
        double d_ij = eucl(thisPos, otherPos);
        double w_ij = pow(d_ij, params.ALPHA);
        terms.push_back(term(i, j, d_ij, w_ij, false));
    }
}
double Rectangle::maxScaleRatio(const Shape* other,const Coord & thisPos, const Coord & otherPos)const{
    const Rectangle * otherRect = dynamic_cast<const Rectangle *>(other);
    double actualDist = eucl(thisPos.x(), thisPos.y(), otherPos.x(), otherPos.y());

    double actualX = thisPos.x() - otherPos.x();
    double actualY = thisPos.y() - otherPos.y();
    double padding = 1e-4;
    double desiredWidth = (this->s.width() + otherRect->s.width()) / 2 + padding;
    double desiredHeight = (this->s.height() + otherRect->s.height()) / 2 + padding;
    
    double widthRatio = desiredWidth / actualX;
    double heightRatio = desiredHeight / actualY;
    
    double unoverlapRatio = min(abs(widthRatio), abs(heightRatio));
    actualX *= unoverlapRatio;
    actualY *= unoverlapRatio;
    
    double optimalDist = vecNorm2D(actualX, actualY);
    return  optimalDist / actualDist;
}
