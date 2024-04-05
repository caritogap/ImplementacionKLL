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
#include <map>
using namespace std;


template <typename T> class KLL_Sketch_d{
    public:
    vector<vector<T>> sketch;
    vector<int> size,sizes;
    vector<int> level;
    int k;
    int H;
    int s;
    float c;
    int H1prima;
    int H2prima;
    float delta,epsilon;
    long long int n;
    int maxH=0;
    bool added_sampler=false,final_compact=false;
    int sampler_height=0;
    pair<T,int> reservoir_insertions=make_pair(0,0);

    map<T,int> item_sorted_map;
    long long int total_weight;

    long long int inels=0;

    vector<pair<int,T>> elements_inserted;
    vector<pair<int,T>> elements_sampler_pos;
    vector<pair<int,T>> elements_sampler_neg;
    vector<pair<int,T>> elements_compacted;

    vector<bool> compacting;
    KLL_Sketch_d(float delta, float epsilon, long long int n, float c);
    void populate(int H,int H1prima,int H2prima,int k,float c);

    void insert_element(T element);
    pair<T,int> reservoir_sample(T element,int w,int h);
    void insert_into_compactor(int i, T element);
    
    void compact(int i);
    int grow();

    void estimate();
    int rank(T element);
    T quantile(float q);
    
    //KLL_Sketch_d<T> merge(KLL_Sketch_d a);
    void clear_sketch(){
        sketch.clear();
        size.clear();
        level.clear();
        item_sorted_map.clear();
        KLL_Sketch_d(delta,epsilon,n,c);
    }
};

template <typename T> pair <T,int> KLL_Sketch_d<T>::reservoir_sample(T element,int w,int h){
    float select=1;
    float prob;
        if(h>0){
            if(reservoir_insertions.second+w>pow(2,h)){
                if(w<reservoir_insertions.second){
                    
                    select=(rand()%10001)/10000.0;
                    int last_weight=reservoir_insertions.second;
                    reservoir_insertions.second=w;
                    reservoir_insertions.first=element;
                    if(select<=((float)last_weight/pow(2,h))){
                        elements_sampler_pos.clear();
                        return make_pair(reservoir_insertions.first,1);
                    }
                    
                }else{
                    select=(rand()%10001)/10000.0;
                    if(select<=((float)w/pow(2,h))){
                        for(int x=0; x<elements_sampler_pos.size(); x++){
                        }
                        elements_sampler_pos.clear();
                        return make_pair(element,1);
                    }
                }
            }

            
            select=(rand()%10001)/10000.0;
            prob=(((float)w/(float)(reservoir_insertions.second+w)));
            reservoir_insertions.second=reservoir_insertions.second+w;
            if(select<=prob){
                reservoir_insertions.first=element;
            }

            if((reservoir_insertions.second)==pow(2,h)){
                reservoir_insertions.second=0;
                for(int x=0; x<elements_sampler_pos.size(); x++){
                }
                elements_sampler_pos.clear();
                return make_pair(reservoir_insertions.first,1);
            }
           
        }else{
               return make_pair(element,1);
        }
    return make_pair(element,-1);
}


template <typename T> void KLL_Sketch_d<T>::estimate(){
    
    int weight;
    total_weight=0;
    item_sorted_map.clear();
    for (int i=0; i<sketch.size(); i++){
        for(int j=0; j<sketch[i].size(); j++){
            weight=pow(2,level[i+sampler_height]-1);            
            total_weight=total_weight+weight;
            item_sorted_map[sketch[i].at(j)]+=weight;
        }
    } 
}  

template <typename T> T KLL_Sketch_d<T>::quantile(float q){
    estimate();
    int ran=0;

    for(auto& i: item_sorted_map){
        ran=ran+i.second;
        if(ran>=(q*total_weight)){
            return i.first;
        }
    }
    return 0;
}

template <typename T> int KLL_Sketch_d<T>::rank(T element){
    estimate();
    int ran=0;
    for(auto& i: item_sorted_map){
        if(i.first<=element){
            ran=ran+i.second;
        }else{
            if(ran>=0){
                return ran;
            }else{
                return 0;
            }
            
        }
    }
        return total_weight;
    
}

template <typename T> KLL_Sketch_d<T>::KLL_Sketch_d(float delta, float epsilon, long long int n, float c){
    H=ceil(log2(epsilon*n));
    s=ceil(log2(log2(1/(epsilon*delta))));
    k=ceil((1/epsilon)*log2(log2(1/(epsilon*delta))));
    H1prima=H-s;
    H2prima=0;
    for(int h=0;h<=log2(n);h++){
        if(ceil(k*pow(c,h))>2){
            if(h<s){
                sizes.push_back(ceil(k));
            }else{
                sizes.push_back(ceil(k*pow(c,h)));
            }
            
            maxH++;
        }else{
            sizes.push_back(2);
        }
    }
    cout<<"maxH: "<<maxH<<endl;
    this->n=n;
    this->delta=delta;
    this->epsilon=epsilon;
    this->c=c;
    cout<<"H= "<<H<<" s="<<s<<" k="<<k<<endl;
    cout<<"H2prima= "<<H2prima<<"  H1prima= "<<H1prima<<" H="<<H<<endl; 
}

template <typename T> void KLL_Sketch_d<T>::populate(int H, int H1prima, int H2prima, int k, float c){
    for(int i=1 ;i<=H; i++){
        if(i>H2prima){
            if(i<=H1prima)
            {
                vector<pair<int,T>> compactor;
                sketch.push_back(compactor);
                level.push_back(i);
                size.push_back(ceil(k*pow(c,H-i)));
            }else{
                vector<pair<int,T>> compactor;
                sketch.push_back(compactor);
                level.push_back(i);
                size.push_back(ceil(k));
            }
        }
    }
    

}
template <typename T> int KLL_Sketch_d<T>::grow(){
    if(sketch.size()>=maxH){
        sampler_height++;
        level.push_back(level.size()+1);
        return 1;
    }else{
        size.insert(size.begin(),sizes[sketch.size()]);
        vector<T> compactor;
        sketch.push_back(compactor);
        level.push_back(level.size()+1);
        return 2;
    }
}

template <typename T> void KLL_Sketch_d<T>::compact(int i){
    if(sketch[i].size()==0){
        return;
    }
    sort(sketch[i].begin(),sketch[i].end());
    compacting.resize(sketch.size(),0);
    compacting[i]=true;
    vector<bool> erased;
    vector<pair<int,T>> temp;
    vector<pair<int,T>> push;
    erased.resize(sketch[i].size(),0);
    
    push.clear();
    int select=rand()%2;

    if(sketch.size()==i+1){
        int cs=grow();
        if(cs==2){
        }else{
            added_sampler=true;
            size.insert(size.begin(),sizes[sketch.size()]);
            vector<T> compactor;
            sketch.push_back(compactor);

            
        }
    }
    
                for(int x=select; x<sketch[i].size(); x=x+2){
                    if(i<H){
                        insert_into_compactor(i+1,sketch[i].at(x));
                    }
                }
    sketch[i].clear();
    compacting[i]=false;
    if(count(compacting.begin(),compacting.end(),true)==0 && final_compact==false){
        if(added_sampler==true){

        final_compact=true;
        compact(0);
        
        vector <T> left=sketch[0];
        sketch.erase(sketch.begin());
        size.erase(size.begin());
        
        for(int i=0; i<left.size();i++){
            pair<T,int> element_to_insert_2=reservoir_sample(left[i],pow(2,sampler_height-1),sampler_height);
                if(element_to_insert_2.second!=-1){
                    insert_into_compactor(0,element_to_insert_2.first);
                }
                

        }
            added_sampler=false;
            final_compact=false;
        }   
        
    }
        
    return;
}


template <typename T> void KLL_Sketch_d<T>::insert_element(T element){

    inels++;
    if(sketch.size()==0){
            grow();
    }

    for(int i=0; i<sketch.size(); i++){
                if(sketch[i].size()>size[i]){
                        compact(i);
                }
    }
    pair<T,int> element_to_insert=reservoir_sample(element,1,sampler_height);
    if(element_to_insert.second!=-1){
        insert_into_compactor(0,element_to_insert.first);
    }
    
    
}

template <typename T> void KLL_Sketch_d<T>::insert_into_compactor(int i,T element){
    if(sketch[i].size()>=size[i]){
        compact(i);
    }
    sketch[i].push_back(element);
}


int main(){
    KLL_Sketch_d<uint32_t> kll_sketch(0.01,0.005,200000000,0.66666);
    std::vector<uint32_t> v;
    std::ifstream myfile ("201911021400_length.txt");
    std::ofstream outputqueries("201911021400_queriesd.txt");
    std::ofstream outputresults("201911021400_resultsd.txt");


    string line;
    int elins=0;
    if (myfile.is_open()) {
    while ( getline (myfile,line)) {
		std::stringstream s(line);
		uint32_t x;
		s>>x;
		v.push_back(x);
        elins++;
    	}
        myfile.close();
        
    }

    for (int i = 0; i < v.size(); i++) {

        kll_sketch.insert_element(v[i]); // mean=0, stddev=1
     }
    

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
    
    for(int i=0; i<11; i++){
        if(i==10){
            rank_queries[10]=max;
        }else{
        rank_queries[i]=((max-min)/10)*i+min;
        }
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
    outputresults<<fixed<<abs(r_est[i]-r[i])/kll_sketch.inels<<setprecision(4);
    }else{
        outputresults<<fixed<<abs(r_est[i]-r[i])/kll_sketch.inels<<setprecision(4);
    }
    outputresults<<"\\\\"<<endl;
   }
       
    outputresults<<"\\hline"<<endl;
    outputresults<<"\\end{longtable}"<<endl;
    outputresults<<"\\end{document}"<<endl;

}

