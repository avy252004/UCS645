#include <chrono>
#include <limits>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

double get_random_number(){
        static thread_local mt19937 mt{random_device{}()};
        static thread_local uniform_real_distribution<double> random_range{0.0,1.0};
        return random_range(mt);
}
int main(){
        int ROW_SIZE{1000};
        int COL_SIZE{1000};
        vector<double> mat_a(ROW_SIZE*COL_SIZE, 0.0);
        vector<double> mat_b(ROW_SIZE*COL_SIZE, 0.0);
        vector<double> mat_b_T(ROW_SIZE * COL_SIZE, 0.0);
        vector<double> res(ROW_SIZE*COL_SIZE, 0.0);


        #pragma omp parallel for collapse(2)
        for(int i = 0;i<ROW_SIZE;++i){
                for(int j = 0;j<COL_SIZE;++j){
                        mat_a[i*COL_SIZE+j] = get_random_number();
                        mat_b[i*COL_SIZE+j] = get_random_number();
                }
        }


        #pragma omp parallel for collapse(2)
        for(int i = 0; i < ROW_SIZE; ++i){
                for(int j = 0; j < COL_SIZE; ++j){
                        mat_b_T[j * ROW_SIZE + i] = mat_b[i * COL_SIZE + j];
                }
        }


        cout<< "Using No Threading \n";
        auto start_time{chrono::steady_clock::now()};
        for(int i = 0;i<ROW_SIZE;++i){
                for(int j = 0;j<COL_SIZE;++j){
                        double sum{0.0};
                        for(int k = 0;k<ROW_SIZE;++k){
                                sum += mat_a[i*COL_SIZE+k]*mat_b_T[j*COL_SIZE+k];
                        }
                        res[i*COL_SIZE+j] = sum;
                }
        }
        auto end_time{chrono::steady_clock::now()};
        auto exec_time{end_time-start_time};
        chrono::duration<double,milli> ms{exec_time};
        cout<<"Execution time:"<<ms.count()<<'\n';


        cout<<"Using 1D Threading \n";
        start_time = chrono::steady_clock::now();
        #pragma omp parallel for
        for(int i = 0;i<ROW_SIZE;++i){
                for(int j = 0;j<COL_SIZE;++j){
                        double sum{0.0};
                        for(int k = 0;k<ROW_SIZE;++k){
                                sum += mat_a[i*COL_SIZE+k]*mat_b_T[j*COL_SIZE+k];
                        }
                        res[i*COL_SIZE+j] = sum;
                }
        }
        end_time = chrono::steady_clock::now();
        exec_time = end_time-start_time;
        ms = exec_time;
        cout<<"Execution time:"<<ms.count()<<'\n';


        cout<<"Using 2D Threading \n";
        start_time = chrono::steady_clock::now();
        #pragma omp parallel for
        for(int i = 0;i<ROW_SIZE;++i){
                #pragma omp parallel for
                for(int j = 0;j<COL_SIZE;++j){
                        double sum{0.0};
                        for(int k = 0;k<ROW_SIZE;++k){
                                sum += mat_a[i*COL_SIZE+k]*mat_b_T[j*COL_SIZE+k];
                        }
                        res[i*COL_SIZE+j] = sum;
                }
        }
        end_time = chrono::steady_clock::now();
        exec_time = end_time-start_time;
        ms = exec_time;
        cout<<"Execution time:"<<ms.count()<<'\n';


        cout<<"Using 2D Threading (collapsed) \n";
        start_time = chrono::steady_clock::now();
        #pragma omp parallel for collapse(2)
        for(int i = 0;i<ROW_SIZE;++i){
                for(int j = 0;j<COL_SIZE;++j){
                        double sum{0.0};
                        for(int k = 0;k<ROW_SIZE;++k){
                                sum += mat_a[i*COL_SIZE+k]*mat_b_T[j*COL_SIZE+k];
                        }
                        res[i*COL_SIZE+j] = sum;
                }
        }
        end_time = chrono::steady_clock::now();
        exec_time = end_time-start_time;
        ms = exec_time;
        cout<<"Execution time:"<<ms.count()<<'\n';
        return 0;
}
