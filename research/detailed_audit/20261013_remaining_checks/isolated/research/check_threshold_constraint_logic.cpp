// Positive and deliberately false controls for the conditional certificate.
#define main constrained_coverage_main
#include "physical_temperature_reduced_combined.cpp"
#undef main
int main(){
    try{
        Box good=initial(0);good[8]=I(1);
        if(!verify_cell(good,0).good)throw std::runtime_error("nominal positive control unresolved");
        Box equal=good;equal[7]=equal[6];
        if(verify_cell(equal,0).good)throw std::runtime_error("false time-saving claim accepted for equal preparations");
        if(verify_cell(equal,1).good)throw std::runtime_error("false distance claim accepted for equal preparations");
        equal[8]=decimal(".02");
        if(!verify_cell(equal,0).good)throw std::runtime_error("outside-band implication unresolved");
        std::cout<<"{\"positive_control\":true,\"equal_preparation_time_claim_rejected\":true,\"equal_preparation_distance_claim_rejected\":true,\"outside_band_control\":true}"<<std::endl;
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}
}
