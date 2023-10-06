#include "Monitorer.hpp"
#include <iostream>

void Monitorer::save(const std::string & filepath){
    std::string sep = " ";
    std::ofstream f(filepath);
    f << "stress_pass" + sep;
    for (int i = 0; i < stress_pass.size(); ++i)
    {
        f << stress_pass[i] << sep;
    }

    f << std::endl
      << "overlap_pass" + sep;
    for (int i = 0; i < overlap_pass.size(); ++i)
    {
        f << overlap_pass[i] << sep;
    }

    f << std::endl
      << "stress_iter" + sep;
    for (int i = 0; i < stress_iter.size(); ++i)
    {
        f << stress_iter[i] << sep;
    }

    f << std::endl
      << "overlap_iter" + sep;
    for (int i = 0; i < overlap_iter.size(); ++i)
    {
        f << overlap_iter[i] << sep;
    }

    f << std::endl
      << "scale_pass" + sep;
    for (int i = 0; i < scale_pass.size(); ++i)
    {
        f << scale_pass[i] << sep;
    }

    f << std::endl
      << "scale_iter" + sep;
    for (int i = 0; i < scale_iter.size(); ++i)
    {
        f << scale_iter[i] << sep;
    }

    f << std::endl
      << "pass_length" + sep;
    for (int i = 0; i < pass_length.size(); ++i)
    {
        f << pass_length[i] << sep;
    }

    f.close();
}
void Monitorer::end_pass(int n){
    pass_length.push_back(n);
}
void Monitorer::stress( const Layout & g,const std::vector<term> & terms,bool iter){
    double s=0.0;
    for(const term &t:terms){
        const double &w_ij = t.w;
        const double &delta_ij = t.d;
        const int &i = t.i, &j = t.j;
        double dx = g.nodes[i]->coord.x() - g.nodes[j]->coord.x(); 
        double dy = g.nodes[i]->coord.y() - g.nodes[j]->coord.y();
        if(t.face_i>0){
            const PolygonShape *iShape = dynamic_cast<const PolygonShape*>(g.nodes[i]->shape);
            const auto & ci = iShape->getFaceCenters(t.face_i-1);
            dx+=ci.x();
            dy+=ci.y();
        }
        if(t.face_j>0){
            const PolygonShape *jShape = dynamic_cast<const PolygonShape*>(g.nodes[j]->shape);
            const auto & cj = jShape->getFaceCenters(t.face_j-1);
            dx-=cj.x();
            dy-=cj.y();
       }

        double d_ij = std::sqrt(dx * dx + dy * dy);

        s+= w_ij * std::pow(d_ij - delta_ij, 2.); 
    }
    if(iter)
        stress_iter.push_back(s/(g.nodes.size()*(g.nodes.size()-1)/2));
    else
        stress_pass.push_back(s/(g.nodes.size()*(g.nodes.size()-1)/2));
}
void Monitorer::overlap(Layout & g,bool iter){
    std::vector<std::tuple<int,int>> overlaps;
    getAllOverlaps(g,overlaps);
    int n = overlaps.size();
    if(iter)
        overlap_iter.push_back(n);
    else
        overlap_pass.push_back(n);
}
void Monitorer::scale(bool iter){
    if(iter)
        scale_iter.push_back(curScale);
    else
        scale_pass.push_back(curScale);
}
void Monitorer::setScale(double s){
    curScale = s;
}
