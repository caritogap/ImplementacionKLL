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
#include "lfsreqs.cpp"
using namespace std;


template <typename T> class KLL_SketchES2{
    public:
    LFSRES lfsr_insert;
    vector <LFSRES> lfsr_compact;
    uint32_t random_number;
    vector<vector<T>> sketch;
    vector<int> size;
    vector<int> level;
    vector<pair<T,int>> estimated_vector;
    pair<T,int> reservoir=make_pair(0,0);
    long long int n,inserted_elements=0;;
    int k;
    int H;
    int s;
    float c;
    int H1prima;
    int H2prima;
    int count=1;


    KLL_SketchES2(float delta, float epsilon, long long int n, float c,int seed);
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


template <typename T> void KLL_SketchES2<T>::reservoir_sample(T element,int w,int h){
    float select=1;
    uint32_t prob;
    random_number=lfsr_insert.new_random_number();
    if(h>1){
        
    if(reservoir.second+w<=pow(2,h-1)){
        select=random_number;
        prob=(uint32_t)(((float)w/(float)(reservoir.second+w))*(1<<13));
        reservoir.second=reservoir.second+w;
        if(select<prob){
            reservoir.first=element;
        }
    }
    
    if((reservoir.second)==pow(2,h-1)){
        insert_into_compactor((h-H2prima-1),reservoir.first);
        reservoir.second=0;
    }else if(reservoir.second+w>pow(2,h-1)){
        if(w<reservoir.second){
            reservoir.second=w;
            select=random_number;
            if(select<=(int)((float)reservoir.second/pow(2,h-1))*(1<<13)){
                insert_into_compactor((h-H2prima-1),reservoir.first);
            }
            reservoir.first=element;
            
            
        }else{
            select=(lfsr_insert.new_random_number()%10001)/10000.0;
            
            if(select<=((float)reservoir.second/pow(2,h-1))){
                insert_into_compactor((h-H2prima-1),element);
            }
        }
    }
    }else{
        insert_into_compactor(0,element);
    }
}
template <typename T> void KLL_SketchES2<T>::estimate(){
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
template <typename T> T KLL_SketchES2<T>::quantile(float q){
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

template <typename T> int KLL_SketchES2<T>::rank(T element){
    estimate();
    int ran=0;
    for(int i=0; i<estimated_vector.size();i++){
        if(estimated_vector[i].first<=element){
            ran=ran+pow(2,estimated_vector[i].second);
        }else{
            return ran;
        }
    }
    return ran;
}

template <typename T> KLL_SketchES2<T>::KLL_SketchES2(float delta, float epsilon, long long int n, float c,int seed){
    LFSRES a(seed);
    LFSRES b0(seed*324);
    LFSRES b1(seed*400);
    LFSRES b2(seed*430);
    LFSRES b3(seed*525);
    LFSRES b4(seed*978);
    LFSRES b5(seed*766);
    LFSRES b6(seed*945);
    LFSRES b7(seed*546);
    LFSRES b8(seed*876);
    LFSRES b9(seed*23);
    
    LFSRES b10(seed*574);
    LFSRES b11(seed*964);
    LFSRES b12(seed*123);
    LFSRES b13(seed*253);
    LFSRES b14(seed*877);
    LFSRES b15(seed*924);
    
    lfsr_compact.push_back(b0);
    lfsr_compact.push_back(b1);
    lfsr_compact.push_back(b2);
    lfsr_compact.push_back(b3);
    lfsr_compact.push_back(b4);
    lfsr_compact.push_back(b5);
    lfsr_compact.push_back(b6);
    lfsr_compact.push_back(b7);
    lfsr_compact.push_back(b8);
    lfsr_compact.push_back(b9);
    lfsr_compact.push_back(b10);
    lfsr_compact.push_back(b11);
    lfsr_compact.push_back(b12);
    lfsr_compact.push_back(b13);
    lfsr_compact.push_back(b14);
    lfsr_compact.push_back(b15);

    lfsr_insert=a;
    

    H=ceil(log2(epsilon*n));
    s=ceil(log2(log2(1/(delta*epsilon))));
    k=ceil(2*(1/epsilon)*log2(log2(1/(delta*epsilon))));
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

template <typename T> void KLL_SketchES2<T>::populate(int H, int H1prima, int H2prima, int k, float c){
    for(int i=1 ;i<=H; i++){
        if(i>H2prima && ceil(k*pow(c,H-i))>2){
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

template <typename T> void KLL_SketchES2<T>::compact(int i){
 
   
    sort(sketch[i].begin(),sketch[i].end());
    int select=lfsr_compact[i].new_random_number();

    select=select%2;
    int start;
    if(sketch[i].size()%2==0){
        if(select==1){
            start=sketch[i].size()-1;

        }else{
            start=sketch[i].size()-2;
        }
    }else{
        if(select==0){
            start=sketch[i].size()-1;

        }else{
            start=sketch[i].size()-2;
        }
    }

    for(int x=start; x>=0; x--){
        
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


template <typename T> void KLL_SketchES2<T>::insert_element(T element){
    inserted_elements++;
    reservoir_sample(element,1,H2prima+1);
    
}

template <typename T> void KLL_SketchES2<T>::insert_into_compactor(int i,T element){

        if(sketch[i].size()==size[i]){
            compact(i);
        }
    
    

        sketch[i].push_back(element);


    
}

int main(int argc, char **argv){
    KLL_SketchES2<uint32_t> kll_sketch(0.01,0.005,200000000,0.66666,2138840);

    std::vector<uint32_t> v;


    std::ifstream myfile (argv[1]);
    std::ifstream myfile2 (argv[2]);
    std::ofstream outputqueries(argv[3]);
    std::ofstream outputresults(argv[4]);

    string line;
    int elins=0;
    if (myfile.is_open()) {
    while ( getline (myfile,line)) {
		//float x = std::stof( line );
		//float x = atof( line.c_str() );
		std::stringstream s(line);
		//printf("x %s %f %d \n", line.c_str(), x, *ending);
		uint32_t x;
		s>>x;
		//printf("x %s %lf \n", line.c_str(), x);
		v.push_back(x);
        elins++;
    	}
        myfile.close();
        
    }

    for (int i = 0; i < v.size(); i++) {

        kll_sketch.insert_element(v[i]); // mean=0, stddev=1
        //cout<<v[i]<<endl;
     }
    
    //     cout<<"Rank 0: "<<kll_sketch.rank(0)<<endl;
    //     cout<<"Rank 300: "<<kll_sketch.rank(300)<<endl;
    //     cout<<"Rank 1000: "<<kll_sketch.rank(1000)<<endl;
    //     cout<<"Rank 5000: "<< kll_sketch.rank(5000)<<endl;
    //     cout<<"Rank 15000: "<< kll_sketch.rank(15000)<<endl;

    cout<<"FINAL"<<endl;

    for(int x=0; x<kll_sketch.sketch.size(); x++){
            cout<<"Compactor "<<x<<endl;
            sort(kll_sketch.sketch[x].begin(),kll_sketch.sketch[x].end());
            for(int y=0; y<kll_sketch.sketch[x].size(); y++){
                cout<<std::setw(10)<<kll_sketch.sketch[x].at(y)<<endl;
            }
     }


    
    auto end_insert= std::chrono::high_resolution_clock::now();


    auto end_query= std::chrono::high_resolution_clock::now();

    sort(v.begin(),v.end());



    outputresults<<"\\"<<"documentclass{article}"<<endl;
    outputresults<<"\\usepackage{longtable}"<<endl;
    outputresults<<"\\usepackage{lscape}"<<endl;
    outputresults<<"\\begin{document}"<<endl;
    outputresults<<"\\begin{longtable}[c]{|c|c|}"<<endl;
    outputresults<<"\\hline"<<endl;
    outputresults<<"\\textbf{Number} & \\textbf{Error rank}\\\\"<<endl; 
    outputresults<<"\\hline"<<endl;
    uint32_t rank_queries[11];
    uint32_t min=v[0],max=v[v.size()-1];
    
    string quer;
     
    for(int i=0; i<11; i++){
        getline (myfile2,quer);
        std::stringstream qu(quer);
        // if(i==10){
        //     rank_queries[10]=max;
        // }else{
        // rank_queries[i]=((max-min)/10)*i+min;
        // }
        qu>>rank_queries[i];
        outputqueries<<rank_queries[i]<<endl;
    }

    double r[11],r_est[11];
    cout<<v.size()<<endl;
    for(int x=0; x<11; x++){
        r[x]=0;
        for(int y=0; y<v.size(); y++){
                
                if(v[y]<=rank_queries[x]){
                    if(y==(v.size()-1)){
                        r[x]=v.size();
                    }
                }else{
                    cout<<y<<" "<<v[y]<<endl;
                    r[x]=y;
                    cout<<x<<" "<<fixed<<r[x]<<endl;
                    break;
                }
            
        }
        
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
    outputqueries<<fixed<<r_est[i]<<endl;
    q_est[i]=kll_sketch.quantile(0.1*i);
   }
    for(int i=0; i<11; i++){
    outputresults<<i<<" & ";
    outputqueries<<fixed<<r[i]<<endl;
    if(r[i]!=0){
    outputresults<<fixed<<abs(r_est[i]-r[i])/kll_sketch.inserted_elements<<setprecision(4);
    }else{
        outputresults<<fixed<<abs(r_est[i]-r[i])/kll_sketch.inserted_elements<<setprecision(4);
    }
    //outputresults<<" & "<< abs(q_est[i]-q[i])
    outputresults<<"\\\\"<<endl;
   }
       
    outputresults<<"\\hline"<<endl;
    outputresults<<"\\end{longtable}"<<endl;
    outputresults<<"\\end{document}"<<endl;

}


