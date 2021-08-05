#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>

double f(double x, double y)
{
    return sin(x)*sin(y)-exp(-pow(x+M_PI/2.0, 2)-pow(y-M_PI/2.0,2));
}

double d_rand(const double min_v, const double max_v)
{
    double r = (double)rand() / RAND_MAX;
    r = r * (max_v-min_v) + min_v;
    return r;
}

double trans(double a)
{
    double da = d_rand(-1,1);
    return fabs(a + da) <= 10.0 ? da : (a + da > 10.0 ? 10.0 - a : -10.0 - a);   
}

int main(void)
{
    srand(time(0));
    const int N = 200;
    double **coords = new double*[N];
    for(int i=0; i<N; i++)
    {
        coords[i] = new double[2];
    }
    double *values = new double[N];


    for(int i=0; i<N; i++)
    {
        coords[i][0] = 5.0;
        coords[i][1] = 5.0;
        values[i] = 0.0;
    }

    std::ofstream w0("w0.dat");
    std::ofstream T_file("T.dat");
 
    for(int it=0; it<=20; it++)
    {
        double T = 10.0/(pow(2,it));
        for(int k=0; k<100; k++)
        {
            for(int i=0; i<N; i++)
            {
                 double dx = trans(coords[i][0]);
                 double dy = trans(coords[i][1]);
                 if(f(coords[i][0]+dx, coords[i][1]+dy) < f(coords[i][0], coords[i][1]))
                 {
                     coords[i][0] += dx;
                     coords[i][1] += dy;
                 }
                 else if (d_rand(0,1) < exp(-(f(coords[i][0]+dx, coords[i][1]+dy)-f(coords[i][0],coords[i][1]))/T))
                 {
                     coords[i][0] += dx;
                     coords[i][1] += dy;
                 }
                 values[i] = f(coords[i][0], coords[i][1]); 
            }
            w0 << values[0] << std::endl;    
        }
        if(it==0 || it==7 || it==20)
        {
            std::cout << "it = " << it << ", T = " << T << std::endl;
            for(int i=0; i<N; i++)
            {
                T_file << coords[i][0] << " " << coords[i][1] << std::endl;
            }
            T_file << std::endl << std::endl;
        }
    }

    w0.close();
    T_file.close();

    double temp = values[0];
    int min_index = 0;
    for(int i=1; i<N; i++)
    {  
        if(values[i] < temp)
        {
            temp = values[i];
            min_index = i;
        }
    }

    std::cout << "Minimum funkcji f: " << values[min_index] << " dla (x_min, y_min) = (" << coords[min_index][0] << ", " << coords[min_index][1] << ")" << std::endl; 

    for(int i=0; i<N; i++)
    {
        delete [] coords[i];
    }
    delete [] coords;
    delete [] values;
    
}
