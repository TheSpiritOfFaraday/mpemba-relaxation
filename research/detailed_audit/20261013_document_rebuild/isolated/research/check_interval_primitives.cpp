// Boundary regressions against 512-bit reference arithmetic, not a universal proof.
#define main physical_coverage_main
#include "physical_temperature_reduced_combined.cpp"
#undef main
int main(){
 try{
  int checks=0;mpfr_t a,b,r,z,t;for(auto p:{a,b,r,z,t})mpfr_init2(p,512);
  auto enclosed=[&](I x){if(mpfr_cmp_d(r,x.l)<0||mpfr_cmp_d(r,x.h)>0)throw std::runtime_error("reference escaped enclosure");checks++;};
  std::vector<double> values={0.,1.,-1.,.1,-.1,std::nextafter(1.,0.),std::numeric_limits<double>::min(),std::numeric_limits<double>::denorm_min()};
  for(double x:values)for(double y:values){
   mpfr_set_d(a,x,MPFR_RNDN);mpfr_set_d(b,y,MPFR_RNDN);
   mpfr_add(r,a,b,MPFR_RNDN);enclosed(I(x)+I(y));
   mpfr_sub(r,a,b,MPFR_RNDN);enclosed(I(x)-I(y));
   mpfr_mul(r,a,b,MPFR_RNDN);enclosed(I(x)*I(y));
   if(y!=0){mpfr_div(r,a,b,MPFR_RNDN);enclosed(I(x)/I(y));}
  }
  for(double x:{-700.,-1.,0.,1.,700.}){mpfr_set_d(a,x,MPFR_RNDN);mpfr_exp(r,a,MPFR_RNDN);enclosed(expi(I(x)));}
  for(double x:{std::numeric_limits<double>::denorm_min(),.1,1.,2.,100.}){
   mpfr_set_d(a,x,MPFR_RNDN);mpfr_log(r,a,MPFR_RNDN);enclosed(logi(I(x)));mpfr_sqrt(r,a,MPFR_RNDN);enclosed(sqrti(I(x)));
  }
  for(double x:{0.,1e-12,std::nextafter(1e-8,0.),1e-8,std::nextafter(1e-8,INFINITY),1.,100.}){
   mpfr_set_d(z,x,MPFR_RNDN);mpfr_sqrt(a,z,MPFR_RNDN);
   for(int which=0;which<3;which++){
    if(x==0){mpfr_set_ui(r,1,MPFR_RNDN);if(which==2)mpfr_div_ui(r,r,6,MPFR_RNDN);}
    else if(which==0)mpfr_cosh(r,a,MPFR_RNDN);
    else if(which==1){mpfr_sinh(r,a,MPFR_RNDN);mpfr_div(r,r,a,MPFR_RNDN);}
    else {mpfr_cosh(r,a,MPFR_RNDN);mpfr_sinh(t,a,MPFR_RNDN);mpfr_div(t,t,a,MPFR_RNDN);mpfr_sub(r,r,t,MPFR_RNDN);mpfr_div(r,r,z,MPFR_RNDN);mpfr_div_ui(r,r,2,MPFR_RNDN);}
    enclosed(entire(I(x),which));
   }
  }
  I range=phii(I(-.9,.5));if(range.l!=0)throw std::runtime_error("phi minimum not zero");checks++;
  for(double x:{-.9,-.1,0.,.2,.5}){mpfr_set_d(a,x,MPFR_RNDN);mpfr_add_ui(b,a,1,MPFR_RNDN);mpfr_log(r,b,MPFR_RNDN);mpfr_mul(r,r,b,MPFR_RNDN);mpfr_sub(r,r,a,MPFR_RNDN);enclosed(range);}
  for(int k=0;k<5;k++){
   bool rejected=false;try{
    if(k==0)(void)(I(1)/I(-1,1));
    if(k==1)(void)logi(I(0,1));
    if(k==2)(void)sqrti(I(-1,1));
    if(k==3)(void)cap(I(0,1),I(2,3));
    if(k==4)(void)phii(I(-1,0));
   }catch(const std::exception&){rejected=true;}
   if(!rejected)throw std::runtime_error("invalid enclosure accepted");checks++;
  }
  for(auto p:{a,b,r,z,t})mpfr_clear(p);
  std::cout<<"{\"passed\":true,\"checks\":"<<checks<<",\"reference_bits\":512,\"scope\":\"Primitive and entire-function boundary regression points; analytic enclosure proof remains separate\"}"<<std::endl;
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
}
