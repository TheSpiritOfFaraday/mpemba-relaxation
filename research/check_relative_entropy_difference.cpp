#define main relative_entropy_coverage_main
#include "physical_box_relative_entropy.cpp"
#undef main
int main(int argc,char**argv){
    try{
        if(argc!=10)throw std::runtime_error("expected normalized time and eight parameters");
        Box b;b[8]=decimal(argv[1]);
        for(int k=0;k<8;k++){b[k]=decimal(argv[k+2]);if(k>=6)b[k]=b[k]-I(1);}
        Model m(b);J F=m.difference(b[8]);
        std::cout<<std::setprecision(17)<<"{\"value\":["<<F.v.l<<","<<F.v.h<<"],\"gradient\":[";
        for(int k=0;k<9;k++){if(k)std::cout<<",";std::cout<<"["<<F.d[k].l<<","<<F.d[k].h<<"]";}
        std::cout<<"]}"<<std::endl;return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}
}
