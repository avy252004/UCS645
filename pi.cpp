#include <chrono>
#include <iostream>
#include <ratio>
using namespace std;

int main(){
        int num_steps{static_cast<int>(1e9)};
        double pi{};
        double sum{0.0};
        double step{1.0/num_steps};

        cout<<"With No Threading \n";
        auto start_time{chrono::steady_clock::now()};
        for(int i = 0;i<num_steps;++i){
                double x = (i+0.5)*step;
                sum += 4.0/(1.0+x*x);
        }
        pi = step*sum;
        auto end_time{std::chrono::steady_clock::now()};
        auto exec_time{end_time-start_time};
        chrono::duration<double, std::milli> ms{exec_time};
        cout<<"Execution Time:"<<ms.count()<<'\n';
        cout<<"Pi (Serial):"<<pi<<'\n';

        sum = 0.0;
        cout<<"With Threading \n";
        start_time = chrono::steady_clock::now();
        #pragma omp parallel for reduction(+:sum)
        for(int i = 0;i<num_steps;++i){
                double x = (i+0.5)*step;
                sum += 4.0/(1.0+x*x);
        }
        pi = step*sum;
        end_time = chrono::steady_clock::now();
        exec_time = end_time-start_time;
        ms = exec_time;
        cout<<"Execution Time:"<<ms.count()<<'\n';
        cout<<"Pi (Parallel):"<<pi<<'\n';
        return 0;
}
