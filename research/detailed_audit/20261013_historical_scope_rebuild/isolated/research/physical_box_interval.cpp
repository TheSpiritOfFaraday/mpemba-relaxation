// Experimental interval/automatic-derivative engine for the uncertainty box.
// Build: clang++ -O2 -std=c++17 -ffp-contract=off -I/opt/homebrew/include
//   research/physical_box_interval.cpp -L/opt/homebrew/lib -lmpfr -lgmp
//   -o tmp/physical_box_interval
//
// Basic IEEE binary64 operations are widened by one nextafter in each
// direction. Transcendentals/decimal input use MPFR with directed rounding.
// No fast-math. The executable rejects unresolved signs; it is NOT a
// robustness certificate unless its entire coverage run completes.
#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpfr.h>
constexpr int ND=9;
double dn(double x){return std::nextafter(x,-INFINITY);}
double up(double x){return std::nextafter(x,INFINITY);}
struct I {
    double l,h;
    I(double x=0):l(x),h(x){}
    I(double a,double b):l(a),h(b){if(!(a<=b))throw std::runtime_error("invalid interval");}
};
I operator+(I a,I b){return I(dn(a.l+b.l),up(a.h+b.h));}
I operator-(I a){return I(-a.h,-a.l);}
I operator-(I a,I b){return a+(-b);}
I operator*(I a,I b){
    double z[4]={a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h};
    return I(dn(*std::min_element(z,z+4)),up(*std::max_element(z,z+4)));
}
I operator/(I a,I b){
    if(b.l<=0&&b.h>=0)throw std::runtime_error("division through zero");
    return a*I(dn(1/b.h),up(1/b.l));
}
I cap(I a,I b){return I(std::max(a.l,b.l),std::min(a.h,b.h));}
double mid(I a){return a.l+(a.h-a.l)/2;}
double rad(I a){return (a.h-a.l)/2;}
double mag(I a){return std::max(std::abs(a.l),std::abs(a.h));}
double trans(double x,bool exponential,mpfr_rnd_t rnd){
    mpfr_t a;mpfr_init2(a,128);mpfr_set_d(a,x,MPFR_RNDN);
    if(exponential)mpfr_exp(a,a,rnd);else mpfr_log(a,a,rnd);
    double y=mpfr_get_d(a,rnd);mpfr_clear(a);return y;
}
I expi(I a){return I(trans(a.l,true,MPFR_RNDD),trans(a.h,true,MPFR_RNDU));}
I logi(I a){
    if(a.l<=0)throw std::runtime_error("nonpositive log argument");
    return I(trans(a.l,false,MPFR_RNDD),trans(a.h,false,MPFR_RNDU));
}
I decimal(const char*s){
    mpfr_t a;mpfr_init2(a,128);
    mpfr_set_str(a,s,10,MPFR_RNDD);double lo=mpfr_get_d(a,MPFR_RNDD);
    mpfr_set_str(a,s,10,MPFR_RNDU);double hi=mpfr_get_d(a,MPFR_RNDU);
    mpfr_clear(a);return I(lo,hi);
}
// phi(z)=(1+z)log(1+z)-z. Monotone on each side of zero.
I phi_point(double x){
    I z(x);return (I(1)+z)*logi(I(1)+z)-z;
}
I phii(I z){
    if(z.l<=-1)throw std::runtime_error("unresolved probability positivity");
    I a=phi_point(z.l),b=phi_point(z.h);
    double low=z.l<=0&&z.h>=0?0:std::min(a.l,b.l);
    return I(std::max(0.,low),std::max(a.h,b.h));
}
struct J{
    I v;
    std::array<I,ND> d{};
    J(double x=0):v(x){}
    J(I x):v(x){}
    static J var(I x,int k){J a(x);a.d[k]=I(1);return a;}
};
J operator+(J a,J b){J c(a.v+b.v);for(int k=0;k<ND;k++)c.d[k]=a.d[k]+b.d[k];return c;}
J operator-(J a){J c(-a.v);for(int k=0;k<ND;k++)c.d[k]=-a.d[k];return c;}
J operator-(J a,J b){return a+(-b);}
J operator*(J a,J b){J c(a.v*b.v);for(int k=0;k<ND;k++)c.d[k]=a.d[k]*b.v+a.v*b.d[k];return c;}
J operator/(J a,J b){J c(a.v/b.v);for(int k=0;k<ND;k++)c.d[k]=(a.d[k]-c.v*b.d[k])/b.v;return c;}
J expj(J a){J c(expi(a.v));for(int k=0;k<ND;k++)c.d[k]=c.v*a.d[k];return c;}
J phij(J a){J c(phii(a.v));I g=logi(I(1)+a.v);for(int k=0;k<ND;k++)c.d[k]=g*a.d[k];return c;}
using V=std::array<J,4>;
using Mat=std::array<std::array<J,3>,3>;
Mat identity(){Mat m{};for(int i=0;i<3;i++)m[i][i]=J(1);return m;}
Mat mul(const Mat&a,const Mat&b){
    Mat c{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)c[i][j]=c[i][j]+a[i][k]*b[k][j];return c;
}
Mat basis(){
    const char* entries[3][3]={
        {"0.8096550290462836","-0.9995127825310614","-0.0013154541523163457"},
        {"-0.4384607226250695","-0.0038635126342375535","-0.7070467305505806"},
        {"-0.39014219030395","-0.03097209755797053","0.707165603234566"}};
    Mat K{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)K[i][j]=J(decimal(entries[i][j]));
    return K;
}
Mat inverse(Mat a){
    Mat b=identity();
    for(int k=0;k<3;k++){
        J pivot=a[k][k];
        for(int j=0;j<3;j++){a[k][j]=a[k][j]/pivot;b[k][j]=b[k][j]/pivot;}
        for(int i=0;i<3;i++)if(i!=k){
            J x=a[i][k];
            for(int j=0;j<3;j++){a[i][j]=a[i][j]-x*a[k][j];b[i][j]=b[i][j]-x*b[k][j];}
        }
    }
    return b;
}
double norm(const Mat&a,int derivative=-1){
    double n=0;
    for(int i=0;i<3;i++){double s=0;for(int j=0;j<3;j++)s=up(s+mag(derivative<0?a[i][j].v:a[i][j].d[derivative]));n=std::max(n,s);}
    return n;
}
Mat expm(Mat A){
    int scaling=0;
    while(norm(A)>.125){
        for(auto&row:A)for(auto&x:row)x=x/J(2);
        scaling++;if(scaling>20)throw std::runtime_error("exponential scaling limit");
    }
    Mat total=identity(),term=identity();
    const int order=16;
    for(int k=1;k<=order;k++){
        term=mul(term,A);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            term[i][j]=term[i][j]/J(k);total[i][j]=total[i][j]+term[i][j];
        }
    }
    I r(norm(A)),power(1);
    for(int k=1;k<=order;k++)power=power*r/I(k);
    I derivative_tail=expi(r)*power;
    I value_tail=derivative_tail*r/I(order+1);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        total[i][j].v=total[i][j].v+I(-value_tail.h,value_tail.h);
        for(int k=0;k<ND;k++){
            I b=derivative_tail*I(norm(A,k));
            total[i][j].d[k]=total[i][j].d[k]+I(-b.h,b.h);
        }
    }
    for(int k=0;k<scaling;k++)total=mul(total,total);
    return total;
}
V gibbs(const V&E,J T){
    V q{};J sum=0;
    for(int i=0;i<4;i++){q[i]=expj(-E[i]/T);sum=sum+q[i];}
    for(auto&a:q)a=a/sum;return q;
}
using Box=std::array<I,9>; // gaps, barriers, temperature excesses, KL threshold
void tighten(J&j,I at_center,const Box&box){
    I enclosure=at_center;
    for(int k=0;k<8;k++)enclosure=enclosure+j.d[k]*(box[k]-I(mid(box[k])));
    j.v=cap(j.v,enclosure);
}
struct Model{
    V pi,initial[2];Mat R{},K{},Kinv{};
    Model(const Box&box,bool use_mean=true){
        std::array<J,8>x;
        for(int k=0;k<8;k++)x[k]=J::var(box[k],k);
        V E={J(0),x[0],x[0]+x[1],x[0]+x[1]+x[2]};
        std::array<std::array<J,4>,4>M{};
        for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){
            J h=(i==0&&j==3)?x[3]:(i==1&&j==2)?x[4]:x[5];
            M[i][j]=expj(-h);M[j][i]=expj(-h-(E[j]-E[i]));
        }
        for(int j=0;j<4;j++)for(int i=0;i<4;i++)if(i!=j)M[j][j]=M[j][j]-M[i][j];
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)R[i][j]=M[i][j]-M[i][3];
        K=basis();Kinv=inverse(K);R=mul(mul(Kinv,R),K);
        pi=gibbs(E,J(1));initial[0]=gibbs(E,J(1)+x[6]);initial[1]=gibbs(E,J(1)+x[7]);
        if(use_mean){
            Box c=box;for(int k=0;k<8;k++)c[k]=I(mid(box[k]));
            Model point(c,false);
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)tighten(R[i][j],point.R[i][j].v,box);
            for(int i=0;i<4;i++){
                tighten(pi[i],point.pi[i].v,box);
                for(int b=0;b<2;b++)tighten(initial[b][i],point.initial[b][i].v,box);
            }
        }
    }
    std::array<J,2> eval(I time)const{
        J t=J::var(time,8);Mat A=R;
        for(auto&row:A)for(auto&z:row)z=z*t;
        Mat P=mul(mul(K,expm(A)),Kinv);std::array<J,2>D{};
        for(int b=0;b<2;b++){
            V dev{};
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)dev[i]=dev[i]+P[i][j]*(initial[b][j]-pi[j]);
            dev[3]=-(dev[0]+dev[1]+dev[2]);
            for(int i=0;i<4;i++)D[b]=D[b]+pi[i]*phij(dev[i]/pi[i]);
        }
        return D;
    }
};
Box center(const Box&b){Box c;for(int k=0;k<9;k++)c[k]=I(mid(b[k]));return c;}
I mean_value(I nominal,const J&box_eval,const Box&box){
    I value=nominal;
    for(int k=0;k<8;k++)value=value+box_eval.d[k]*(box[k]-I(mid(box[k])));
    return cap(value,box_eval.v);
}
double central_root(const Model&point,double delta,int b){
    double lo=0,hi=600;
    for(int i=0;i<34;i++){
        double t=(lo+hi)/2;I D=point.eval(I(t))[b].v;
        if(mid(D)>delta)lo=t;else hi=t;
    }
    return (lo+hi)/2;
}
struct Root{I value;std::array<I,9> d;};
Root root(const Box&box,const Model&family,const Model&point,int b){
    double c=central_root(point,mid(box[8]),b);
    double lo=.65*c,hi=1.35*c+1;
    auto bound=[&](double t){
        auto d=family.eval(I(t))[b];auto dc=point.eval(I(t))[b];
        return mean_value(dc.v,d,box);
    };
    for(int i=0;i<6&&bound(lo).l<=box[8].h;i++)lo*=.7;
    for(int i=0;i<6&&bound(hi).h>=box[8].l;i++)hi*=1.3;
    if(!(bound(lo).l>box[8].h&&bound(hi).h<box[8].l))throw std::runtime_error("root bracket unresolved");
    // Separate monotone searches bound the earliest and latest roots over
    // the whole family before differentiating on their enclosing interval.
    double left=lo,right=hi;
    for(int i=0;i<22;i++){
        double m=(left+right)/2;
        if(bound(m).l>box[8].h)left=m;else right=m;
    }
    double earliest=left;
    left=lo;right=hi;
    for(int i=0;i<22;i++){
        double m=(left+right)/2;
        if(bound(m).h<box[8].l)right=m;else left=m;
    }
    lo=earliest;hi=right;
    I T(lo,hi);
    for(int i=0;i<6;i++){
        auto d=family.eval(T)[b];
        if(d.d[8].h>=0)break;
        double m=mid(T);
        I next=cap(T,I(m)-(bound(m)-box[8])/d.d[8]);
        if(next.h-next.l>.999*(T.h-T.l))break;
        T=next;
    }
    auto d=family.eval(T)[b];
    if(d.d[8].h>=0)throw std::runtime_error("root derivative unresolved");
    Root r;r.value=T;
    for(int k=0;k<8;k++)r.d[k]=-d.d[k]/d.d[8];
    r.d[8]=I(1)/d.d[8];
    return r;
}
struct Margin{I value;std::array<I,9>d;};
Margin margin(const Box&box){
    Box c=center(box);Model family(box),point(c);
    Root cold=root(box,family,point,0),hot=root(box,family,point,1);
    Root cc=root(c,point,point,0),hc=root(c,point,point,1);
    Margin m;m.value=decimal(".95")*cold.value-hot.value;
    I centered=decimal(".95")*cc.value-hc.value;
    for(int k=0;k<9;k++){
        m.d[k]=decimal(".95")*cold.d[k]-hot.d[k];
        centered=centered+m.d[k]*(box[k]-I(mid(box[k])));
    }
    m.value=cap(m.value,centered);return m;
}
Box initial(double radius){
    const char*params[]={"0.8308942048583045","0.11635274930182392","2.729136958784003","2.11406392019699","2.822657578561883","5.640080424208037","3.8694077398491555","9.472085161836484"};
    Box b;
    for(int k=0;k<8;k++){
        I x=decimal(params[k]);if(k>=6)x=x-I(1);
        b[k]=x*(I(1)+I(-radius,radius));
    }
    b[8]=decimal(".001");return b;
}
int main(int argc,char**argv){
    try{
        double radius=argc>1?std::stod(argv[1]):0;
        Box box=initial(radius);
        if(argc>2)box[8]=decimal(argv[2]);
        std::cout<<std::setprecision(17);
        Model family(box);
        auto test=family.eval(I(60));
        std::cout<<"KL at t=60: cold ["<<test[0].v.l<<","<<test[0].v.h<<"] hot ["<<test[1].v.l<<","<<test[1].v.h<<"]\n";
        Margin m=margin(box);
        std::cout<<"time margin ["<<m.value.l<<","<<m.value.h<<"]\n";
        for(int k=0;k<9;k++)std::cout<<"derivative "<<k<<" ["<<m.d[k].l<<","<<m.d[k].h<<"]\n";
        return 0;
    }catch(const std::exception&e){std::cerr<<"UNRESOLVED: "<<e.what()<<"\n";return 2;}
}
