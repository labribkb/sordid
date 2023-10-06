#ifndef MONITORER_H
#define MONITORER_H
#include <vector>
#include <string>
#include "mystructs.hpp"

class Monitorer{
private:
    std::vector<double> stress_pass;
    std::vector<int> overlap_pass;
    std::vector<double> stress_iter;
    std::vector<int> overlap_iter;
    std::vector<double> scale_pass;
    std::vector<double> scale_iter;
    std::vector<int> pass_length;
    double curScale;
public:
    void save(const std::string & filepath);
    void end_pass(int n);
    void stress(const Layout & g,const std::vector<term> & orig_terms,bool iter);
    void overlap( Layout & g,bool iter);
    void scale(bool iter);
    void setScale(double s);
};
#endif
