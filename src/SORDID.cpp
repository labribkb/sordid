#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "SORDID.hpp"
#include "Monitorer.hpp"

using namespace std;
vector<term> layoutToTerms(Layout &g, Layout &init_g, Parametrizer &params)
{
    vector<term> terms;
    vector<vector<term>*> termsMT;
    for(int i=0;i<omp_get_max_threads();++i)
        termsMT.push_back(new vector<term>());
    double d_ij, w_ij; //, d_eucl;
    double xi, xj, yi, yj;
    bool overlap;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < g.N(); i++)
    {
        for (int j = i + 1; j < g.N(); j++)
        {
            if (i != j)
            {
                g.nodes[i]->terms(g.nodes[j],params,i,j,*termsMT[omp_get_thread_num()]);
            }
        }
    }
    for(unsigned int i=0;i<termsMT.size();++i){
        terms.insert(terms.end(),termsMT[i]->begin(),termsMT[i]->end());
        delete termsMT[i];
    }
    return terms;
}

// S_GD2 optim algorithm, adapted from https://github.com/jxz12/s_gd2/blob/master/cpp/s_gd2/layout.cpp
void OPTIMIZATION_PASS(Layout &g, Layout &init_g, vector<term> &terms, const vector<double> &etas, Parametrizer &params,Monitorer & monitor)
{
    rk_state rstate;
    rk_seed(params.seed, &rstate);
    double mvt_sum;
    unsigned int i_eta;
    for (i_eta = 0; i_eta < etas.size(); i_eta++)
    {
        const double eta = etas[i_eta];
        unsigned n_terms = terms.size();
        if (n_terms == 0)
            return;
        fisheryates_shuffle(terms, rstate);

        mvt_sum = 0;
        for (unsigned i_term = 0; i_term < n_terms; i_term++)
        {
            const term &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &w_ij = t.w;
            const double &d_ij = t.d;
            if (true || t.o)
            {
                double mu = eta * w_ij;
                if (mu > 1)
                    mu = 1;

                double dx = g.nodes[i]->coord.x() - g.nodes[j]->coord.x(); // X[i * 2] - X[j * 2];
                double dy = g.nodes[i]->coord.y() - g.nodes[j]->coord.y(); // X[i * 2 + 1] - X[j * 2 + 1];
                if(terms[i_term].face_i>0){
                    const PolygonShape *iShape = dynamic_cast<const PolygonShape*>(g.nodes[i]->shape);
                    const auto & ci = iShape->getFaceCenters(terms[i_term].face_i-1);
                    dx+=ci.x();
                    dy+=ci.y();
                }
                if(terms[i_term].face_j>0){
                    const PolygonShape *jShape = dynamic_cast<const PolygonShape*>(g.nodes[j]->shape);
                    const auto & cj = jShape->getFaceCenters(terms[i_term].face_j-1);
                    dx-=cj.x();
                    dy-=cj.y();
               }

                double mag = sqrt(dx * dx + dy * dy);


                double r = 0;
                if (mag != 0){
                    r = (mu * (mag - d_ij)) / (2 * mag);

                }

                double r_x = r * dx;
                double r_y = r * dy;
                mvt_sum += abs(r_x) + abs(r_y);
                g.nodes[i]->move(-r_x, -r_y);
                g.nodes[j]->move(r_x, r_y);
                if(std::isnan(r_x) || std::isnan(r_y)){
                   cout<<"NAN "<< mag<<" "<<d_ij<<" " <<r_x<< " "<<r_y<<" "<<mu<<" "<<w_ij<<" "<<eta<<endl;
                }
            }
        }
        if (params.MONITOR)
        {
            monitor.stress(g,terms,true);
            monitor.overlap(g,true);
            monitor.scale(true);
        }
        if (mvt_sum < params.MINIMUM_MOVEMENT)
        {
            if (params.MONITOR)
                monitor.end_pass(i_eta + 1);
            return;
        }
        terms = layoutToTerms(g, init_g, params);

    }
    if (params.MONITOR)
        monitor.end_pass(i_eta);
    return;
}

void passInOptim(Layout &g, Layout &init_g, Parametrizer &params,Monitorer & monitor)
{
    vector<term> orig_terms = layoutToTerms(g, init_g, params);
    vector<double> etas = schedule(orig_terms, params.MAX_ITER, params.eps);
    OPTIMIZATION_PASS(g, init_g, orig_terms, etas, params,monitor);
    if (params.MONITOR)
    {
        monitor.stress(g,orig_terms,false);
        monitor.overlap(g,false);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 10)
    {
        cerr << "Wrong parameters" << endl;
        return EXIT_FAILURE;
    }
    string input_layout_path = argv[1];
    Layout g;
    Parametrizer params;
    Monitorer monitor;
    loadParams(params, argv);
    g.loadLayout(input_layout_path);

    int n_passes = 0;
    double scaleFactor;
    auto start = std::chrono::steady_clock::now();
    g.do_triangulate();
    Layout init_g(g);
    
    double maxRatio = g.maxScaleRatio();
    double upperScale = maxRatio;
    double lowerScale = 1.;
    double oldScale = 1.;
    double curScale = 1.;
    monitor.setScale(curScale);
    bool stop = !scanLineOverlapCheck(g); // do not enter if there is no overlap
    Layout best(g);
    bool foundSolution=false;
    // if possible, try to solve the problem in current scale

    if (!stop && g.isCurrentScaleSolvable())
    {
        passInOptim(g, init_g, params,monitor);
        stop = !scanLineOverlapCheck(g);
        if(stop)
            best = g;
        if (params.MONITOR)
            monitor.scale(false);
    }

    vector<tuple<int, int>> overlaps;
    while (!stop)
    {
        curScale = (upperScale + lowerScale) / 2;
        if (params.MONITOR)
            monitor.setScale(curScale);
        cout << "pass " << n_passes << " ; scale : " << curScale << endl;
        scaleFactor = curScale / oldScale;
        oldScale = curScale;
        init_g.scale(scaleFactor);
        if (params.PRIME)
            g = Layout(init_g);
        else
            g.scale(scaleFactor);

        passInOptim(g, init_g, params,monitor);
    
        if (params.MONITOR)
            monitor.scale(false);
        if (scanLineOverlapCheck(g)){
            lowerScale = curScale;
            if (upperScale!=maxRatio && upperScale - lowerScale < params.SCALE_STEP){
                stop = true;
            }
        }
        else // no overlap
        {
            if (upperScale - lowerScale < params.SCALE_STEP){
                stop = true;
            }else
                upperScale = curScale;
            best = Layout(g);
            foundSolution = true;
        }

        ++n_passes;
        if (n_passes >= params.MAX_PASSES){
            if(!foundSolution){
                best.scale(maxRatio);
            }
            break;
	}
    }
    auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    string savepath = best.save(input_layout_path, "sordid", elapsed, params.PRIME);
   if(params.MONITOR)
       monitor.save(savepath +".stats");


    cout << "DONE; saved to " << savepath << endl;

    return EXIT_SUCCESS;
}
