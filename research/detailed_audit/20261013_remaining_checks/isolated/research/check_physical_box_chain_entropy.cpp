// Read-only point/derivative adapter for an independent high-precision audit.
#define main coverage_main
#include "physical_box_chain_entropy.cpp"
#undef main
int main(int argc,char**argv){
    try{
        if(argc!=10)throw std::runtime_error("expected normalized time and eight parameters");
        Box b;
        b[8]=decimal(argv[1]);
        for(int k=0;k<8;k++){b[k]=decimal(argv[k+2]);if(k>=6)b[k]=b[k]-I(1);}
        Model m(b);auto D=m.eval(b[8]);
        std::cout<<std::setprecision(17)<<"{\"slow\":["<<m.spectral->slow.v.l<<","<<m.spectral->slow.v.h<<"],\"KL\":[";
        for(int j=0;j<2;j++){
            if(j)std::cout<<",";
            std::cout<<"{\"value\":["<<D[j].v.l<<","<<D[j].v.h<<"],\"gradient\":[";
            for(int k=0;k<9;k++){if(k)std::cout<<",";std::cout<<"["<<D[j].d[k].l<<","<<D[j].d[k].h<<"]";}
            std::cout<<"]}";
        }
        std::cout<<"]}"<<std::endl;return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}
}
