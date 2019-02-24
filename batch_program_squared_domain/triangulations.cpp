#include "triangulations.h"
#include "computations.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <time.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void regular_tri (DrawableTrimesh<> & m,int N)
{
    double deltax=100.0/N;                    // edge length
    double deltay=deltax*sqrt(3)/2;         // height of tris
    int Nh = (int)(100.0/deltay);             // number of rows
    double shift = 0; //commented to stretch tris to y=1        (1.0-Nh*deltay)/2.0;     // brings (0.5,0.5) in the middle
    double offset;                          // offset for odd rows
    double rescale = 100.0/(Nh * deltay);

    std::vector<double> points={0,shift,100,shift,100,Nh*deltay+shift,0,Nh*deltay+shift};
    std::vector<uint>segs={0,1,1,2,2,3,3,0};

    for (int i=0;i<=Nh;++i)                 // y coords
    {
        offset = (i%2!=0) ? deltax/2 : 0;   // shift at odd rows
        if (i%2!=0 && i<Nh) {               // additional point at odd rows
            points.push_back(0);
            points.push_back(i*deltay+shift);
        }
        for(int j=0;j<=N;++j)               // x coords
        {
           if((i==0 && j==0)||(i==Nh && j==N)||(i==Nh && j+offset==0)||(i==0 && j==N))
               continue; // skip corners
           points.push_back(deltax*j+(j==N?0:offset));   // rightmost point within bound
           points.push_back(i*deltay+shift);
        }
    }

    for (int i=1;i<points.size();i+=2) points[i]*=rescale;     // bring y in range [0,1]
    std::string flags;
    triangle_wrap(points,segs, {}, 0, flags.c_str(), m);
    m.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void anisotropic_tri(DrawableTrimesh<> & m, const double max_area, const double min_angle, const double y_anisotropy)
{
    std::string flags = "Qq" + std::to_string(min_angle) + "a" + std::to_string(max_area);
    triangle_wrap({0,0,100*y_anisotropy,0,100*y_anisotropy,100,0,100}, {0,1,1,2,2,3,3,0}, {}, 0, flags.c_str(), m);
    for (int i=0;i<m.num_verts();i++) m.vert(i).x() /= y_anisotropy;
    m.update_bbox();

    m.updateGL();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double my_rand()
{
    double r;
    do r = (double)rand()/RAND_MAX; while(r==0.0 || r==1.0);
    return r;
}

void non_unif_Poisson2D(int n, std::vector<double> & buf)
// non-uniform Poisson disc sampling of n points in [0,1]x[0,1]
// with non-uniform density inv. prop. to distance of point from origin
{
    double scale = n/3;
    int edgesamples = 10*(int)sqrt(n);
    int nsamples = 40000-edgesamples;
    int count = 0;
    double x, y, r, rj, d;

    std::vector<double> bigbuf;
    std::vector<bool> flagbuf(nsamples,true);

    for (int i=0;i<edgesamples;i++) {bigbuf.push_back(my_rand()); bigbuf.push_back(0);}
    for (int i=0;i<edgesamples;i++) {bigbuf.push_back(my_rand()); bigbuf.push_back(1);}
    for (int i=0;i<edgesamples;i++) {bigbuf.push_back(0); bigbuf.push_back(my_rand());}
    for (int i=0;i<edgesamples;i++) {bigbuf.push_back(1); bigbuf.push_back(my_rand());}
    for (int i=0;i<2*nsamples;i++) bigbuf.push_back(my_rand());

    for (int i=0;i<bigbuf.size();i+=2)
    {
        if (!flagbuf[i/2]) continue;
        x = bigbuf[i];y=bigbuf[i+1];
        r = (x*x+y*y)/scale;
        for (int j=0;j<bigbuf.size();j+=2)
        {
            if (i==j) continue;
            if (!flagbuf[j/2]) continue;
            rj = (bigbuf[j]*bigbuf[j]+bigbuf[j+1]*bigbuf[j+1])/scale;
            d = ((x-bigbuf[j])*(x-bigbuf[j])+(y-bigbuf[j+1])*(y-bigbuf[j+1]));
            if (d<std::min(r,rj)) flagbuf[j/2]=false;
        }
    }

    count = 0;
    for (int i=0;i<flagbuf.size();i++)
    {
        if (flagbuf[i]) {buf.push_back(bigbuf[2*i]); buf.push_back(bigbuf[2*i+1]); count++;}
        if (count>=n) break;
    }
    //std::cout << "Sampled " << count << " points out of desired " << n << std::endl;
}


void non_uniform_tri(DrawableTrimesh<> & m,int N)
{
    double min_angle = 2.0;                // this could be a parameter
    std::vector<double> points={0,0,100,0,100,100,0,1};
    std::vector<uint>segs={0,1,1,2,2,3,3,0};

    non_unif_Poisson2D(N*N,points);

    std::string flags; //= "Qq" + std::to_string(min_angle);
    triangle_wrap(points,segs, {}, 0, flags.c_str(), m);
    m.updateGL();
}


void make_domain(DrawableTrimesh<> & m, const double max_area, const double min_angle)
{
    std::string flags = "Qq" + std::to_string(min_angle) + "a" + std::to_string(max_area);
    triangle_wrap({0,0,100,0,100,100,0,100}, {0,1,1,2,2,3,3,0}, {}, 0, flags.c_str(), m);
    m.updateGL();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void make_grid(DrawableTrimesh<> &m,int N)
{
    double delta=100.0/N;
    std::vector<double> points={0,0,0,100,100,100,100,0};
    std::vector<uint>segs={0,1,1,2,2,3,3,0};
    for(int i=0;i<=N;++i)
    {
        for(int j=0;j<=N;++j)
        {
            if((i==0 && j==0)||(i==0 && j==N) || (i==N && j==0) || (i==N && j==N))
            {continue;}
                points.push_back(i*delta);
                points.push_back(j*delta);
        }



    }
    std::string flags;
    triangle_wrap(points,segs, {}, 0, flags.c_str(), m);
    m.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void make_triangulation(DrawableTrimesh<> &m, int mode, int N, double y_anisotropy)
{
   //const double y_anisotropy = 5;   // should be set via GUI
   const double min_angle = 20;     // for Delaunay and anisotropic

   switch(mode)
   {
   case 0:
       regular_tri(m,N);
       break;
   case 1:
       make_domain(m,10000/(1.6*pow(N,2)));
       break;
   case 2:
       non_uniform_tri(m,N);
       break;
   case 3:
       anisotropic_tri(m,10000*y_anisotropy/(1.6*pow(N,2)),min_angle,y_anisotropy);
       break;
   }
}
