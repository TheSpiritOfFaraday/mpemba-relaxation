#define main temperature_coverage_main
#include "physical_temperature_reduced_combined.cpp"
#undef main
int main(int argc,char**argv){
    try{
        if(argc!=10)throw std::runtime_error("expected normalized time and eight parameters");
        Box box;box[8]=decimal(argv[1]);for(int k=0;k<8;k++){box[k]=decimal(argv[k+2]);if(k>=6)box[k]=box[k]-I(1);}
        Model model(box);std::cout<<std::setprecision(17)<<"[";
        for(int b=0;b<2;b++){
            if(b)std::cout<<",";std::cout<<"[";auto d=model.thermal_derivatives(box[8],b);
            for(int n=0;n<3;n++){
                if(n)std::cout<<",";std::cout<<"{\"value\":["<<d[n].v.l<<","<<d[n].v.h<<"],\"gradient\":[";
                for(int k=0;k<9;k++){if(k)std::cout<<",";std::cout<<"["<<d[n].d[k].l<<","<<d[n].d[k].h<<"]";}
                std::cout<<"]}";
            }
            std::cout<<"]";
        }
        std::cout<<"]"<<std::endl;
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}
}
