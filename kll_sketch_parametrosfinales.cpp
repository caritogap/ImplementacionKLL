#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
using namespace std;


template <typename T> class KLL_Sketch{
    public:
    vector<vector<T>> sketch;
    vector<int> size;
    vector<int> level;
    vector<pair<T,int>> estimated_vector;
    pair<T,int> reservoir=make_pair(0,0);
    long long int n;
    int k;
    int H;
    int s;
    float c;
    int H1prima;
    int H2prima;
    int count=1;
    KLL_Sketch(float delta, float epsilon, long long int n, float c);
    void reservoir_sample(T element,int w,int h);
    void populate(int H,int H1prima,int H2prima,int k,float c);
    void insert_element(T element);
    int rank(T element);
    T quantile(float q);
    void estimate();
    void compact(int i);
    void insert_into_compactor(int i, T element);
    void clear_sketch(){
        sketch.clear();
        size.clear();
        level.clear();
        estimated_vector.clear();
    }
};


template <typename T> void KLL_Sketch<T>::reservoir_sample(T element,int w,int h){
    float select=1;
    float prob;
    if(h>1){
    if(reservoir.second+w<pow(2,h-1)){
        select=(rand()%10001)/10000.0;
        prob=(((float)w/(float)(reservoir.second+w)));
        reservoir.second=reservoir.second+w;
        if(select<=prob){
            reservoir.first=element;
        }
    }else if((reservoir.second+w)==pow(2,h-1)){
        insert_into_compactor((h-H2prima-1),reservoir.first);
        reservoir.second=0;
    }else if(reservoir.second+w>pow(2,h-1)){
        if(w<reservoir.second){
            reservoir.second=w;
            select=(rand()%10001)/10000.0;
            if(select<=((float)reservoir.second/pow(2,h-1))){
                insert_into_compactor((h-H2prima-1),reservoir.first);
            }
            reservoir.first=element;
            
            
        }else{
            select=(rand()%10001)/10000.0;
            if(select<=((float)reservoir.second/pow(2,h-1))){
                insert_into_compactor((h-H2prima-1),element);
            }
        }
    }
    }else{
        insert_into_compactor(0,element);
    }
}
template <typename T> void KLL_Sketch<T>::estimate(){
    estimated_vector.clear();
    n=0;
    for (int i=0; i<sketch.size(); i++){
        for(int j=0; j<sketch[i].size(); j++){
            estimated_vector.push_back(make_pair(sketch[i].at(j),level[i]-1));
            n=n+pow(2,level[i]-1);
        }
    }
    sort(estimated_vector.begin(),estimated_vector.end());
    
}  
template <typename T> T KLL_Sketch<T>::quantile(float q){
    estimate();
    
    int ran=0;
    for(int i=0; i<estimated_vector.size();i++){

        ran=ran+pow(2,estimated_vector[i].second);

        if(ran>=(q*n)){
            return estimated_vector[i].first;
        }
    }
    return 0;
}

template <typename T> int KLL_Sketch<T>::rank(T element){
    estimate();
    int ran=0;
    for(int i=0; i<estimated_vector.size();i++){
        if(estimated_vector[i].first<=element){
            ran=ran+pow(2,estimated_vector[i].second);
        }else{
            return ran;
        }
    }
    return n;
}

template <typename T> KLL_Sketch<T>::KLL_Sketch(float delta, float epsilon, long long int n, float c){
    H=log2(epsilon*n);
    s=log2(log2(1/delta));
    k=ceil(2*(1/epsilon)*log2(log2(1/delta)));
    H1prima=H-s;
    H2prima=0;
    for(int h=1;h<=H;h++){
        if(ceil(k*pow(c,(H-h)))<=2){
            H2prima++;
        }
    }
    cout<<"H= "<<H<<" s="<<s<<" k="<<k<<endl;
    cout<<"H2prima= "<<H2prima<<"  H1prima= "<<H1prima<<" H="<<H<<endl; 
    populate(H,H1prima,H2prima,k,c);
}

template <typename T> void KLL_Sketch<T>::populate(int H, int H1prima, int H2prima, int k, float c){
    for(int i=1 ;i<=H; i++){
        if(i>H2prima){
            if(i<=H1prima)
            {
                vector<T> compactor;
                sketch.push_back(compactor);
                level.push_back(i);
                size.push_back(ceil(k*pow(c,H-i)));
            }else{
                vector<T> compactor;
                sketch.push_back(compactor);
                level.push_back(i);
                size.push_back(ceil(k));
            }
        }
    }
}

template <typename T> void KLL_Sketch<T>::compact(int i){
    sort(sketch[i].begin(),sketch[i].end());
    int select=rand()%2;
    for(int x=0; x<sketch[i].size(); x++){
        if(select==0 && x%2==0){
            if(i<(H)){
            insert_into_compactor(i+1,sketch[i].at(x));
            }
            } else if(select==1 && x%2==1){
                if(i<(H)){
                    insert_into_compactor(i+1,sketch[i].at(x));
                }
            }
    }
    sketch[i].clear();
}


template <typename T> void KLL_Sketch<T>::insert_element(T element){
    reservoir_sample(element,1,H2prima+1);
}

template <typename T> void KLL_Sketch<T>::insert_into_compactor(int i,T element){
    if(sketch[i].size()==size[i]){
        compact(i);
    }
    sketch[i].push_back(element);
}


int main(int argc, char **argv){
    auto start = std::chrono::high_resolution_clock::now();
    float epsilon,delta,c;
    long long int n;
    epsilon=0.002;
    delta=0.01;
    c=0.66666;
    n=4000000000;
                    
    srand(3312455);
    std::vector<uint32_t> v;
    std::vector<uint32_t> v2=v;
    std::string line;
    std::ifstream myfile (argv[1]);
    std::ofstream outputfile(argv[2]);
    int totalelems=0;
    if (myfile.is_open()) {
        while ( getline (myfile,line) ) {
		std::stringstream s(line);
		uint32_t x;
		s>>x;
		v.push_back(x);
		totalelems++;
    	}
    	myfile.close();
    }
    auto end_read = std::chrono::high_resolution_clock::now();
    KLL_Sketch<uint32_t> kll_sketch(delta,epsilon,n,c);
    for (int i = 0; i < v.size(); i++) {
      kll_sketch.insert_element(v[i]); // mean=0, stddev=1
      
    }
    
    auto end_insert= std::chrono::high_resolution_clock::now();

    
    printf("%d \n",kll_sketch.quantile(0));
    cout<<"Rank 0: "<<kll_sketch.rank(0)<<endl;
    cout<<"Rank 300: "<<kll_sketch.rank(300)<<endl;
    cout<<"Rank 1000: "<<kll_sketch.rank(1000)<<endl;
    cout<<"Rank 5000: "<< kll_sketch.rank(5000)<<endl;
    cout<<"Rank 15000: "<< kll_sketch.rank(15000)<<endl;


    
    //cout<<"Quantile 0: "<<v[0]<<endl;
    printf("%d \n",kll_sketch.quantile(0));
    //cout<<"Quantile 10: "<<v[(int)v.size()/10-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.10));
     //cout<<"Quantile 20: "<<v[(int)(2*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.20));
     //cout<<"Quantile 30: "<<v[(int)(3*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.30));
     //cout<<"Quantile 40: "<<v[(int)(4*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.40));
    //cout<<"Quantile 50: "<<v[(int)(5*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.50));
    //cout<<"Quantile 60: "<<v[(int)(6*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.60));
    //cout<<"Quantile 70: "<<v[(int)(7*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.70));
    //cout<<"Quantile 80: "<<v[(int)(8*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.80));
    //cout<<"Quantile 90: "<<v[(int)(9*v.size()/10)-1]<<endl;
    printf("%d \n",kll_sketch.quantile(0.90));
    //cout<<"Quantile 100: "<< v[v.size()-1]<<endl;
    printf("%d \n",kll_sketch.quantile(1));

    auto end_query= std::chrono::high_resolution_clock::now();

    sort(v.begin(),v.end());
cout<<"Quantile 0: "<<v[0]<<endl;
    cout<<"Quantile 10: "<<v[ceil(v.size()/10)-1]<<endl;
     cout<<"Quantile 20: "<<v[ceil(2*v.size()/10)-1]<<endl;
     cout<<"Quantile 30: "<<v[ceil(3*v.size()/10)-1]<<endl;
     cout<<"Quantile 40: "<<v[ceil(4*v.size()/10)-1]<<endl;
    cout<<"Quantile 50: "<<v[ceil(5*v.size()/10)-1]<<endl;
    cout<<"Quantile 60: "<<v[ceil(6*v.size()/10)-1]<<endl;
    cout<<"Quantile 70: "<<v[ceil(7*v.size()/10)-1]<<endl;
    cout<<"Quantile 80: "<<v[ceil(8*v.size()/10)-1]<<endl;
    cout<<"Quantile 90: "<<v[ceil(9*v.size()/10)-1]<<endl;
    cout<<"Quantile 100: "<< v[v.size()-1]<<endl;
    
    auto time_read = std::chrono::duration_cast<std::chrono::milliseconds>(end_read-start);
    auto time_insert = std::chrono::duration_cast<std::chrono::milliseconds>(end_insert-end_read);
    auto time_query =std::chrono::duration_cast<std::chrono::milliseconds>(end_query-end_insert);
    cout<<time_read.count()<<" "<<time_insert.count()<<" "<<time_query.count()<<endl;
    cout<<time_read.count()+time_insert.count()+time_query.count()<<endl;

    outputfile<<"\\"<<"documentclass{article}"<<endl;
    outputfile<<"\\usepackage{longtable}"<<endl;
    outputfile<<"\\usepackage{lscape}"<<endl;
    outputfile<<"\\begin{document}"<<endl;
    outputfile<<"\\begin{longtable}[c]{|c|c|c|}"<<endl;
    outputfile<<"\\hline"<<endl;
    outputfile<<"\\textbf{Number} & \\textbf{Error rank} & \\textbf{Error quantile} \\\\"<<endl; 
    outputfile<<"\\hline"<<endl;
    uint32_t rank_queries[11];
    uint32_t min=v[0],max=v[v.size()-1];
    
    for(int i=0; i<11; i++){
        if(i==10){
            rank_queries[10]=max;
        }else{
        rank_queries[i]=((max-min)/10)*i+min;
        }
    }

    float r[11],r_est[11];
    
    for(int x=0; x<11; x++){
        vector<uint32_t> sequence={rank_queries[x]};
        vector<uint32_t>::iterator pos=find_end(v.begin(),v.end(),sequence.begin(),sequence.end());
        r[x]=pos-v.begin();
    }  
    float q[11],q_est[11];
    int ksk,Hsk,H1sk,H2sk,capacity;
    for(int x=0; x<11; x++){
        if(x==0){
            q[0]=v[0];
        }else{
            q[x]=v[ceil(x*v.size()/10)-1];
        }
    }   
   for(int i=0; i<11; i++){
    r_est[i]=kll_sketch.rank(rank_queries[i]);
    q_est[i]=kll_sketch.quantile(0.1*i);
   }
    for(int i=0; i<11; i++){
    outputfile<<i<<" & ";
    if(r[i]!=0){
    outputfile<<fixed<<abs(r_est[i]-r[i])/r[i]<<setprecision(4);
    }else{
        outputfile<<fixed<<abs(r_est[i]-r[i])<<setprecision(4);
    }
    outputfile<<" & "<< abs(q_est[i]-q[i])/q[i]<<"\\\\"<<endl;
   }
       
    outputfile<<"\\hline"<<endl;
     outputfile<<"\\end{longtable}"<<endl;
     outputfile<<"\\end{document}"<<endl;




}

