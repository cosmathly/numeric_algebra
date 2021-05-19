#include <iostream>
#include <math.h>
#include <stdio.h>
#define maxsize 100
using namespace std;    
double epsilon, a, h;
double epsilon_arr[] = {1.0, 0.1, 0.01, 0.0001}; 
double A[maxsize][maxsize];
double b[maxsize];
double cur_y[maxsize], tmp_y[maxsize], real_y[maxsize], ans_y[maxsize];
const double e = 2.718281828459045; // 自然底数
const double w = 0.5; // 松弛因子，低松弛迭代
int m, n;
int loop_cnt;
double f(double var)
{
       return (1.0-a)*(1.0-pow(e, -1.0*var/epsilon))/(1.0-pow(e, -1.0/epsilon))+a*var;
}
void get_real_y()
{
     for(int i = 1; i <= n; i++)
     real_y[i] = f((double)(i)*h);
}
void initial()
{
     m = 100;
     a = 0.5;
     n = m-1;
     h = 1.0/(double)(m);
     for(int i = 1; i <= n-1; i++)
     b[i] = a*h*h;
     get_real_y();
}
void get_data()
{
     for(int i = 1; i <= n; i++)
     for(int j = 1; j <= n; j++)
     {
         if(i==j) A[i][j] = -1.0*(2.0*epsilon+h);
         else if(i==j+1) A[i][j] = epsilon;
         else if(j==i+1) A[i][j] = epsilon+h;
         else A[i][j] = 0.0;
     }
     b[n] = a*h*h-epsilon-h;
}
void initial_cur_y()
{
     for(int i = 1; i <= n; i++)
     cur_y[i] = 0.0;
}
double abs(double var)
{
       return var>0.0?var:(-1.0*var);
}
double calc_error()
{
     double sum1 = 0.0;
     for(int i = 1; i <= n; i++)
     sum1 += abs(ans_y[i]-real_y[i]);
     double sum2 = 0.0;
     for(int i = 1; i <= n; i++)
     sum2 += abs(real_y[i]);
     return sum1/sum2;
}
void Jacobi()
{
     loop_cnt = 10000;
     initial_cur_y();
     double sum;
     while(loop_cnt>0)
     {
           for(int i = 1; i <= n; i++)
           {
               sum = -1.0*b[i];
               for(int j = 1; j <= n; j++)
               if(j!=i) sum += (A[i][j]*cur_y[j]);
               tmp_y[i] = sum/(-1.0*A[i][i]);
           }
           for(int i = 1; i <= n; i++)
           cur_y[i] = tmp_y[i];
           loop_cnt--;
     }
     for(int i = 1; i <= n; i++)
     ans_y[i] = cur_y[i];
}
void G_S()
{
     loop_cnt = 10000;
     initial_cur_y();
     double sum;
     while(loop_cnt>0)
     {
           for(int i = 1; i <= n; i++)
           {
               sum = -1.0*b[i];
               for(int j = 1; j <= n; j++)
               if(j!=i) sum += (A[i][j]*cur_y[j]);
               cur_y[i] = sum/(-1.0*A[i][i]);
           }
           loop_cnt--;
     }
     for(int i = 1; i <= n; i++)
     ans_y[i] = cur_y[i];
}
void SOR()
{
     loop_cnt = 10000;
     initial_cur_y();
     double sum;
     while(loop_cnt>0)
     {
           for(int i = 1; i <= n; i++)
           {
               sum = -1.0*b[i];
               for(int j = 1; j <= n; j++)
               if(j!=i) sum += (A[i][j]*cur_y[j]);
               cur_y[i] = ((1.0-w)*cur_y[i]-(w*sum/A[i][i]));
           }
           loop_cnt--;
     }
     for(int i = 1; i <= n; i++)
     ans_y[i] = cur_y[i];    
}
void solve()
{ 
     Jacobi();
     cout<<"1) Jacobi:"<<endl;
     for(int i = 1; i <= n; i++)
     cout<<ans_y[i]<<" ";
     cout<<endl;
     cout<<"相对误差:"<<endl;
     cout<<calc_error()<<endl;
     G_S();
     cout<<"2) G-S:"<<endl;
     for(int i = 1; i <= n; i++)
     cout<<ans_y[i]<<" ";
     cout<<endl;
     cout<<"相对误差:"<<endl;
     cout<<calc_error()<<endl;
     SOR();
     cout<<"3) SOR:"<<endl;
     for(int i = 1; i <= n; i++)
     cout<<ans_y[i]<<" ";
     cout<<endl;
     cout<<"相对误差:"<<endl;
     cout<<calc_error()<<endl;
}
int main()
{
    freopen("output.txt", "w", stdout);
    initial();
    for(int i = 0; i < 4; i++)
    {
        epsilon = epsilon_arr[i];
        get_data();
        solve();
    }
    fclose(stdout);
    return 0;
}