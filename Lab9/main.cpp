#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>

#define frand() ((double)rand())/(RAND_MAX+1.0)

void print_matrix(double **M, int m, int n)
{
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
            printf("%15g", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

void print_vector(double *V, int n)
{
    for(int i=0; i<n; i++)
    {
            printf("%15g", V[i]);
        printf("\n");
    }
    printf("\n");
}

double f(double x, double x0, double x_max, double x_min, double sigma2)
{
    double a = sin((14*M_PI*x)/(x_max-x_min));
    double b = exp(-pow(x-x0, 2)/(2*sigma2));
    double c = exp(-pow(x+x0,2)/(2*sigma2));
    return a*(b+c);
}

double alpha(double *x, double **phi, int j, int n)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for(int i=0; i<n; i++)
    {
        numerator += x[i]*phi[j][i]*phi[j][i];
        denominator += phi[j][i]*phi[j][i];
    }
    return numerator / denominator;
}


double beta(double *x, double **phi, int j, int n)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for(int i=0; i<n; i++)
    {
        numerator += x[i]*phi[j-1][i]*phi[j][i];
        denominator += phi[j-1][i]*phi[j-1][i];
    }
    return numerator / denominator;
}

double c_j(double *y, double **phi, int j, int n)
{
    double result = 0.0;
    for(int i=0; i<n; i++)
    {
        result += y[i]*phi[j][i];
    }
    return result;
}

double s_j(double **phi, int j, int n)
{
    double result = 0.0;
    for(int i=0; i<n; i++)
    {
        result += phi[j][i]*phi[j][i];
    }
    return result;
}

double F(double *y, double **phi, int k, int n, int m)
{
    double result = 0.0;
    for(int j=0; j<=m; j++)
    {
        result += (c_j(y, phi, j, n)/s_j(phi, j, n)) * phi[j][k];
    }
    return result;
}

int main(void)
{   
    srand(time(0));
    const int N = 201;
    double x_max = 4.0;
    double x_min = -4.0;
    double sigma = (x_max - x_min)/16;
    double x_0 = 2.0;
    double h = (x_max - x_min)/(N - 1.0); 
    double C_rand;
    double *x = new double[N];
    double *y_c = new double[N];
    double *y = new double[N];

    for(int i=0; i<N; i++)
    {
        C_rand = (frand()-0.5)/5.0;
        x[i] =  x_min + i*h;
        y_c[i] = f(x[i], x_0, x_max, x_min, sigma*sigma) + C_rand;
    }

    double **phi = new double*[51];
    for(int i=0; i<51; ++i)
    {
        phi[i] = new double[N];
    }    

    for(int k=0; k<N; k++)
    {
        phi[0][k] = 1.0;
    }

    for(int k=0; k<N; k++)
    {
        phi[1][k] = (x[k] - alpha(x, phi, 0, N))*phi[0][k];
    }

    for(int j=1; j<50; j++)
    {
        double a = alpha(x, phi, j, N);
        double b = beta(x, phi, j, N);
        for(int k=0; k<N; k++)
        {
            phi[j+1][k] = (x[k]-a)*phi[j][k]-b*phi[j-1][k];       
        }
    }

    std::ofstream gram("Gram.dat");
    for(int i=0; i<N; i++)
    {
        gram << std::setw(10) << x[i];
        for(int j=0; j<7; j++)
        {
            gram << std::setw(20) << phi[j][i]/phi[j][0];
        }
        gram << std::endl;
    }
    gram.close();

    // Przypadek z szumem
    std::ofstream pkt("pkt.dat");
    for(int i=0; i<N; i++)
    {
        pkt << std::setw(15) << x[i] << std::setw(15) << y_c[i] << std::endl;
    }
    pkt.close();

    std::ofstream approx("approx.dat");
    for(int m=10; m<=50; m+=20)
    {
        for(int k=0; k<N; k++)
        {
            approx << std::setw(15) << x[k] << std::setw(15) 
                   << F(y_c, phi, k, N, m) << std::endl;
        }
        approx << std::endl << std::endl;
    }
    approx.close();

    // Przypadek bez szumu
    std::ofstream pkt_f("pkt_f.dat");
    for(int i=0; i<N; i++)
    {
        pkt_f << std::setw(15) << x[i] << std::setw(15) << y[i] << std::endl;
    }
    pkt_f.close();

    std::ofstream approx_f("approx_f.dat");
    for(int m=10; m<=50; m+=20)
    {
        for(int k=0; k<N; k++)
        {
            approx_f << std::setw(15) << x[k] << std::setw(15) 
                   << F(y, phi, k, N, m) << std::endl;
        }
        approx_f << std::endl << std::endl;
    }
    approx_f.close();

    delete [] x;
    delete [] y_c;
    delete [] y;
    for(int i=0; i<51; ++i){
        delete [] phi[i];
    }
    delete [] phi;
    return 0;
}