// Spectral interval/automatic-derivative coverage engine for the uncertainty box.
// Build: clang++ -O2 -std=c++17 -ffp-contract=off -I/opt/homebrew/include
//   research/physical_box_spectral.cpp -L/opt/homebrew/lib -lmpfr -lgmp
//   -o tmp/physical_box_spectral
//
// Basic IEEE binary64 operations are widened by one nextafter in each
// direction. Transcendentals/decimal input use MPFR with directed rounding.
// No fast-math. The executable rejects unresolved signs; it is NOT a
// robustness certificate unless its entire coverage run completes.
// cover RADIUS MAXNODES MODE: modes 0=time saving, 1=distance window,
// 2=both preparations above .01 at .02 tau, 3=cold below .001 at 3 tau.
// Time in this engine is u=t/tau, not the attempt-frequency clock.
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
#include <memory>
#include <chrono>
#include <fstream>
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
I sqrti(I a){
    if(a.l<0)throw std::runtime_error("negative sqrt interval");
    mpfr_t x;mpfr_init2(x,128);
    mpfr_set_d(x,a.l,MPFR_RNDN);mpfr_sqrt(x,x,MPFR_RNDD);double lo=mpfr_get_d(x,MPFR_RNDD);
    mpfr_set_d(x,a.h,MPFR_RNDN);mpfr_sqrt(x,x,MPFR_RNDU);double hi=mpfr_get_d(x,MPFR_RNDU);
    mpfr_clear(x);return I(lo,hi);
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
bool center_enabled=false;
std::array<I,8> center_deviation{};
struct J{
    I v,c;
    std::array<I,ND> d{};
    J(double x=0):v(x),c(x){}
    J(I x):v(x),c(x){}
    static J var(I x,int k){J a(x);a.d[k]=I(1);if(k<8)a.c=I(mid(x));return a;}
};
J finish(J a){
    a.c=cap(a.c,a.v);
    if(center_enabled){
        I bound=a.c;
        for(int k=0;k<8;k++)bound=bound+a.d[k]*center_deviation[k];
        a.v=cap(a.v,bound);
    }
    return a;
}
J operator+(J a,J b){J c(a.v+b.v);c.c=a.c+b.c;for(int k=0;k<ND;k++)c.d[k]=a.d[k]+b.d[k];return finish(c);}
J operator-(J a){J c(-a.v);c.c=-a.c;for(int k=0;k<ND;k++)c.d[k]=-a.d[k];return finish(c);}
J operator-(J a,J b){return a+(-b);}
J operator*(J a,J b){J c(a.v*b.v);c.c=a.c*b.c;for(int k=0;k<ND;k++)c.d[k]=a.d[k]*b.v+a.v*b.d[k];return finish(c);}
J operator/(J a,J b){J c(a.v/b.v);c.c=a.c/b.c;for(int k=0;k<ND;k++)c.d[k]=(a.d[k]-c.v*b.d[k])/b.v;return finish(c);}
J expj(J a){J c(expi(a.v));c.c=expi(a.c);for(int k=0;k<ND;k++)c.d[k]=c.v*a.d[k];return finish(c);}
J phij(J a){J c(phii(a.v));c.c=phii(a.c);I g=logi(I(1)+a.v);for(int k=0;k<ND;k++)c.d[k]=g*a.d[k];return finish(c);}
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
J trace(const Mat&A){return A[0][0]+A[1][1]+A[2][2];}
J det(const Mat&A){return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);}
// For z>=0, both functions and their first derivatives have positive
// Taylor coefficients. Directed MPFR evaluates their endpoint bounds.
double entire_point(double z,int which,mpfr_rnd_t rnd){
    if(z<1e-8){
        I Z(z),term(1),s(which==2?1./6:1);
        if(which==2){
            // f1'(z)=sum_{n>=0}(n+1)z^n/(2n+3)!.
            term=I(1)/I(6);s=term;
            for(int n=1;n<=8;n++){term=term*Z*I(n+1)/I(n)/I(2*n+2)/I(2*n+3);s=s+term;}
        }else{
            for(int n=1;n<=8;n++){term=term*Z/I(2*n)/I(2*n+(which==0?-1:1));s=s+term;}
        }
        // For z<1e-8 all omitted positive terms sum to less than 1e-80.
        s=s+I(0,1e-80);return rnd==MPFR_RNDD?s.l:s.h;
    }
    // Compute an interval at high precision, with cancellation included.
    mpfr_t x,lo,hi,a,b,c;for(auto p:{x,lo,hi,a,b,c})mpfr_init2(p,192);
    mpfr_set_d(x,z,MPFR_RNDN);
    mpfr_sqrt(lo,x,MPFR_RNDD);mpfr_sqrt(hi,x,MPFR_RNDU);
    if(which==0){mpfr_cosh(a,lo,MPFR_RNDD);mpfr_cosh(b,hi,MPFR_RNDU);}
    else if(which==1){mpfr_sinh(a,lo,MPFR_RNDD);mpfr_div(a,a,hi,MPFR_RNDD);mpfr_sinh(b,hi,MPFR_RNDU);mpfr_div(b,b,lo,MPFR_RNDU);}
    else{
        mpfr_cosh(a,lo,MPFR_RNDD);mpfr_sinh(c,hi,MPFR_RNDU);mpfr_div(c,c,lo,MPFR_RNDU);mpfr_sub(a,a,c,MPFR_RNDD);mpfr_div(a,a,x,MPFR_RNDD);mpfr_div_ui(a,a,2,MPFR_RNDD);
        mpfr_cosh(b,hi,MPFR_RNDU);mpfr_sinh(c,lo,MPFR_RNDD);mpfr_div(c,c,hi,MPFR_RNDD);mpfr_sub(b,b,c,MPFR_RNDU);mpfr_div(b,b,x,MPFR_RNDU);mpfr_div_ui(b,b,2,MPFR_RNDU);
    }
    double out=mpfr_get_d(rnd==MPFR_RNDD?a:b,rnd);
    for(auto p:{x,lo,hi,a,b,c})mpfr_clear(p);return out;
}
I entire(I z,int which){return I(entire_point(z.l,which,MPFR_RNDD),entire_point(z.h,which,MPFR_RNDU));}
J entirej(J z,int which){
    z.v=cap(z.v,I(0,INFINITY));z.c=cap(z.c,z.v);J ans(entire(z.v,which));ans.c=entire(z.c,which);
    I grad=which==0?entire(z.v,1)/I(2):entire(z.v,2);
    for(int k=0;k<ND;k++)ans.d[k]=grad*z.d[k];return finish(ans);
}
struct Spectral{
    Mat P,Q,NQ,H{};std::array<J,2> v{};J slow,a,d2;
    Spectral(const Mat&R,bool evaluate_center=true){
        J tr=trace(R),de=det(R),minor=0;
        for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)minor=minor+R[i][i]*R[j][j]-R[i][j]*R[j][i];
        auto characteristic=[&](J x){
            J determinant=(x-R[1][1])*(x-R[2][2])-R[1][2]*R[2][1];
            J coupling=R[0][1]*((x-R[2][2])*R[1][0]+R[1][2]*R[2][0])+R[0][2]*(R[2][1]*R[1][0]+(x-R[1][1])*R[2][0]);
            return (x-R[0][0])*determinant-coupling;
        };
        auto poly=[&](I x){return characteristic(J(x)).v;};
        auto deriv=[&](I x){
            return (x-R[1][1].v)*(x-R[2][2].v)-R[1][2].v*R[2][1].v+(x-R[0][0].v)*(I(2)*x-R[1][1].v-R[2][2].v)-R[0][1].v*R[1][0].v-R[0][2].v*R[2][0].v;
        };
        I L(-.015,-.001);
        if(!(poly(I(L.l)).h<0&&poly(I(L.h)).l>0&&deriv(L).l>0))throw std::runtime_error("slow eigenvalue bracket unresolved");
        for(int it=0;it<20;it++){
            double c=mid(L);I next=cap(L,I(c)-poly(I(c))/deriv(L));
            if(next.h-next.l>.9999*(L.h-L.l))break;L=next;
        }
        for(int it=0;it<20;it++){
            I determinant=(L-R[1][1].v)*(L-R[2][2].v)-R[1][2].v*R[2][1].v;
            I coupling=R[0][1].v*((L-R[2][2].v)*R[1][0].v+R[1][2].v*R[2][0].v)+R[0][2].v*(R[2][1].v*R[1][0].v+(L-R[1][1].v)*R[2][0].v);
            I next=cap(L,R[0][0].v+coupling/determinant);
            if(next.h-next.l>.9999*(L.h-L.l))break;L=next;
        }
        slow=J(L);I denom=deriv(L);J ch=characteristic(J(L));
        for(int k=0;k<ND;k++)slow.d[k]=-ch.d[k]/denom;
        if(evaluate_center){
            Mat Rc{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)Rc[i][j]=J(R[i][j].c);
            Spectral at_center(Rc,false);slow.c=at_center.slow.v;
        }
        slow=finish(slow);

        J determinant=(slow-R[1][1])*(slow-R[2][2])-R[1][2]*R[2][1];
        std::array<J,3> right={J(1),((slow-R[2][2])*R[1][0]+R[1][2]*R[2][0])/determinant,(R[2][1]*R[1][0]+(slow-R[1][1])*R[2][0])/determinant};
        std::array<J,3> left={J(1),(R[0][1]*(slow-R[2][2])+R[0][2]*R[2][1])/determinant,(R[0][1]*R[1][2]+R[0][2]*(slow-R[1][1]))/determinant};
        J H11=R[1][1]-right[1]*R[0][1],H12=R[1][2]-right[1]*R[0][2];
        J H21=R[2][1]-right[2]*R[0][1],H22=R[2][2]-right[2]*R[0][2];
        H[0][0]=H11;H[0][1]=H12;H[1][0]=H21;H[1][1]=H22;v[0]=right[1];v[1]=right[2];
        a=(H11+H22)/J(2);J diff=(H11-H22)/J(2);d2=diff*diff+H12*H21;d2.v=cap(d2.v,I(0,INFINITY));
        if((a.v+sqrti(d2.v)).h>=slow.v.l)throw std::runtime_error("slowest-mode isolation unresolved");
        J normalization=0;for(int i=0;i<3;i++)normalization=normalization+left[i]*right[i];
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            P[i][j]=right[i]*left[j]/normalization;
            Q[i][j]=J(i==j)-P[i][j];
        }
        Mat N=R;for(int i=0;i<3;i++)N[i][i]=N[i][i]-a;NQ=mul(N,Q);
    }
    Mat eval(J u)const{
        J t=-u/slow;Mat A{};
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)A[i][j]=H[i][j]*t;
        J delta=(A[0][0]-A[1][1])/J(2);
        J kernel=expj((A[0][0]+A[1][1])/J(2))*entirej(delta*delta,1);
        Mat Z{};Z[0][0]=expj(A[0][0]);Z[1][1]=expj(A[1][1]);
        Z[0][1]=A[0][1]*kernel;Z[1][0]=A[1][0]*kernel;
        // Dyson expansion around the diagonal, retaining orders 0 and 1.
        // The n>=2 tail is bounded in the induced infinity norm.
        I c(std::max(mag(A[0][1].v),mag(A[1][0].v)));
        I em=expi(I(std::max(A[0][0].v.h,A[1][1].v.h))),ec=expi(c);
        I second=ec-I(1)-c,first=ec-I(1);
        double tail=std::max(0.,(em*second).h);
        I cc(std::max(mag(A[0][1].c),mag(A[1][0].c)));
        I emc=expi(I(std::max(A[0][0].c.h,A[1][1].c.h)));
        double center_tail=std::max(0.,(emc*(expi(cc)-I(1)-cc)).h);
        for(int i=0;i<2;i++)for(int j=0;j<2;j++){
            Z[i][j].v=Z[i][j].v+I(-tail,tail);
            Z[i][j].c=Z[i][j].c+I(-center_tail,center_tail);
            for(int k=0;k<ND;k++){
                double dD=std::max(mag(A[0][0].d[k]),mag(A[1][1].d[k]));
                double dC=std::max(mag(A[0][1].d[k]),mag(A[1][0].d[k]));
                double bound=std::max(0.,(em*(I(dD)*second+I(dC)*first)).h);
                Z[i][j].d[k]=Z[i][j].d[k]+I(-bound,bound);
            }
        }
        std::array<std::array<J,3>,2> W={{{-v[0],J(1),J(0)},{-v[1],J(0),J(1)}}};
        Mat out{};J es=expj(-u);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            out[i][j]=es*P[i][j];
            for(int a=0;a<2;a++)for(int b=0;b<2;b++)out[i][j]=out[i][j]+Q[i][a+1]*Z[a][b]*W[b][j];
        }
        return out;
    }

};
V gibbs(const V&E,J T){
    V q{};J sum=0;
    for(int i=0;i<4;i++){q[i]=expj(-E[i]/T);sum=sum+q[i];}
    for(auto&a:q)a=a/sum;return q;
}
using Box=std::array<I,9>; // gaps, barriers, temperature excesses, KL threshold
struct CenterGuard{
    bool old_enabled;std::array<I,8> old_deviation;
    CenterGuard(const Box&box):old_enabled(center_enabled),old_deviation(center_deviation){
        center_enabled=true;for(int k=0;k<8;k++)center_deviation[k]=box[k]-I(mid(box[k]));
    }
    ~CenterGuard(){center_enabled=old_enabled;center_deviation=old_deviation;}
};
void tighten(J&j,I at_center,const Box&box){
    I enclosure=at_center;
    for(int k=0;k<8;k++)enclosure=enclosure+j.d[k]*(box[k]-I(mid(box[k])));
    j.v=cap(j.v,enclosure);j.c=cap(j.c,at_center);
}
struct Model{
    Box domain;
    V pi,initial[2];J basin_A,qA,qB;std::array<J,3> modal[2]{};Mat R{},K{},Kinv{};std::shared_ptr<Spectral> spectral;
    Model(const Box&box,bool use_mean=true):domain(box){
        CenterGuard guard(box);
        std::array<J,8>x;
        for(int k=0;k<8;k++)x[k]=J::var(box[k],k);
        V E={J(0),x[0],x[0]+x[1],x[0]+x[1]+x[2]};
        std::array<std::array<J,4>,4>M{};
        for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){
            J h=(i==0&&j==3)?x[3]:(i==1&&j==2)?x[4]:x[5];
            M[i][j]=expj(-h);M[j][i]=expj(-h-(E[j]-E[i]));
        }
        for(int j=0;j<4;j++)for(int i=0;i<4;i++)if(i!=j)M[j][j]=M[j][j]-M[i][j];
        K=basis();Kinv=inverse(K);
        // Transform each constant edge matrix before multiplying its rate.
        // This preserves cancellations of intrabasin flux in the slow row.
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)if(i!=j){
            Mat edge{};
            for(int a=0;a<3;a++)for(int b=0;b<3;b++)edge[a][b]=J((int(a==i)-int(a==j))*(j==3?-1:int(b==j)));
            Mat C=mul(mul(Kinv,edge),K);
            for(int a=0;a<3;a++)for(int b=0;b<3;b++)R[a][b]=R[a][b]+M[i][j]*C[a][b];
        }
        pi=gibbs(E,J(1));initial[0]=gibbs(E,J(1)+x[6]);initial[1]=gibbs(E,J(1)+x[7]);
        J tail=expj(-E[3]),second=expj(-E[1]),third=expj(-E[2]);
        basin_A=(J(1)+tail)/(J(1)+second+third+tail);
        qA=tail/(J(1)+tail);J ratioB=expj(-x[1]);qB=ratioB/(J(1)+ratioB);
        for(int b=0;b<2;b++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)modal[b][i]=modal[b][i]+Kinv[i][j]*(initial[b][j]-pi[j]);
        if(use_mean){
            Box c=box;for(int k=0;k<8;k++)c[k]=I(mid(box[k]));
            Model point(c,false);
            tighten(basin_A,point.basin_A.v,box);tighten(qA,point.qA.v,box);tighten(qB,point.qB.v,box);
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)tighten(R[i][j],point.R[i][j].v,box);
            for(int b=0;b<2;b++)for(int i=0;i<3;i++)tighten(modal[b][i],point.modal[b][i].v,box);
            for(int i=0;i<4;i++){
                tighten(pi[i],point.pi[i].v,box);
                for(int b=0;b<2;b++)tighten(initial[b][i],point.initial[b][i].v,box);
            }
        }
        spectral=std::make_shared<Spectral>(R);
    }
    std::array<J,2> eval(I time)const{
        CenterGuard guard(domain);
        J t=J::var(time,8);
        Mat P=spectral->eval(t);std::array<J,2>D{};
        for(int b=0;b<2;b++){
            V dev{};std::array<J,3> out{};
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)out[i]=out[i]+P[i][j]*modal[b][j];
            J deltaA=0,contrastA=0,contrastB=0;
            for(int j=0;j<3;j++){
                J ca=-(K[1][j]+K[2][j]),fourth=-(K[0][j]+K[1][j]+K[2][j]);
                deltaA=deltaA+ca*out[j];
                contrastA=contrastA+(fourth-qA*ca)*out[j];
                contrastB=contrastB+(K[2][j]+qB*ca)*out[j];
            }
            auto binaryKL=[&](J q,J deviation){return q*phij(deviation/q)+(J(1)-q)*phij(-deviation/(J(1)-q));};
            J massA=basin_A+deltaA,massB=J(1)-massA;
            D[b]=binaryKL(basin_A,deltaA)+massA*binaryKL(qA,contrastA/massA)+massB*binaryKL(qB,contrastB/massB);
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
    double lo=0,hi=5;
    for(int i=0;i<34;i++){
        double t=(lo+hi)/2;I D=point.eval(I(t))[b].v;
        if(mid(D)>delta)lo=t;else hi=t;
    }
    return (lo+hi)/2;
}
struct Root{I value;std::array<I,9> d;};
Root root(const Box&box,const Model&family,const Model&point,int b){
    double c=central_root(point,mid(box[8]),b);
    double lo=.65*c,hi=1.35*c+.01;
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
        b[k]=x*(I(1)+I(-up(radius),up(radius)));
    }
    b[8]=decimal(".001");return b;
}
struct Check{bool good=false;std::array<double,9> score{};std::string reason;};
Check verify_cell(const Box&box,int mode){
    Check result;Model family(box),point(center(box));
    auto D=family.eval(box[8]);auto Dc=point.eval(I(mid(box[8])));
    auto tightened=[&](J value,I center_value){
        I enclosure=center_value;
        for(int k=0;k<9;k++)enclosure=enclosure+value.d[k]*(box[k]-I(mid(box[k])));
        value.v=cap(value.v,enclosure);return value;
    };
    J cold=tightened(D[0],Dc[0].v),hot=tightened(D[1],Dc[1].v);
    if(mode==0&&(cold.v.l>decimal(".01").h||cold.v.h<decimal(".001").l)){result.good=true;result.reason="outside threshold band";return result;}
    if(mode==2){result.good=cold.v.l>decimal(".01").h&&hot.v.l>decimal(".01").h;return result;}
    if(mode==3){result.good=cold.v.h<decimal(".001").l;return result;}
    J F;int selected_kind=0;double selected_multiplier=0;
    if(mode==0){
        auto H=family.eval(decimal(".95")*box[8])[1];
        H.d[8]=H.d[8]*decimal(".95");
        auto Hc=point.eval(decimal(".95")*I(mid(box[8])))[1];
        J original=D[0]-H;I at_center=Dc[0].v-Hc.v;
        F=tightened(original,at_center);
        // On the relevant band, both Dc-.01 and .001-Dc are <=0.
        // A nonnegative multiplier therefore yields a sufficient margin.
        // Its value is only a search heuristic; every chosen margin is
        // evaluated with outward interval arithmetic on the whole box.
        auto upper_constraint=[&](double multiplier){
            if(!(multiplier>=0&&multiplier<=32&&std::isfinite(multiplier)))return;
            J trial=tightened((J(1)+J(multiplier))*D[0]-H-J(multiplier)*J(decimal(".01")),(I(1)+I(multiplier))*Dc[0].v-Hc.v-I(multiplier)*decimal(".01"));
            if(trial.v.l>F.v.l){F=trial;selected_kind=1;selected_multiplier=multiplier;}
        };
        auto lower_constraint=[&](double multiplier){
            if(!(multiplier>=0&&multiplier<=32&&std::isfinite(multiplier)))return;
            J trial=tightened((J(1)-J(multiplier))*D[0]-H+J(multiplier)*J(decimal(".001")),(I(1)-I(multiplier))*Dc[0].v-Hc.v+I(multiplier)*decimal(".001"));
            if(trial.v.l>F.v.l){F=trial;selected_kind=2;selected_multiplier=multiplier;}
        };
        for(double x:{.025,.05,.1,.2,.4,.6,1.,2.,4.}){upper_constraint(x);lower_constraint(x);}
        double slope=mid((Dc[0]-Hc).d[8])/mid(Dc[0].d[8]);
        upper_constraint(-slope);lower_constraint(slope);
        for(int k=0;k<ND;k++){
            double ratio=mid(original.d[k])/mid(D[0].d[k]);
            upper_constraint(-ratio);lower_constraint(ratio);
            double central=mid((Dc[0]-Hc).d[k])/mid(Dc[0].d[k]);
            upper_constraint(-central);lower_constraint(central);
        }
    }else F=tightened(J(decimal(".9"))*D[0]-D[1],decimal(".9")*Dc[0].v-Dc[1].v);
    if(mode==0&&F.v.l<0){
        Box face=box;int collapsed=0;
        for(int k=0;k<9;k++)if(box[k].l<box[k].h){
            if(F.d[k].l>=0){face[k]=I(box[k].l);collapsed++;}
            else if(F.d[k].h<=0){face[k]=I(box[k].h);collapsed++;}
        }
        if(collapsed){
            auto expr=[&](J c,J h){
                if(selected_kind==1)return (J(1)+J(selected_multiplier))*c-h-J(selected_multiplier)*J(decimal(".01"));
                if(selected_kind==2)return (J(1)-J(selected_multiplier))*c-h+J(selected_multiplier)*J(decimal(".001"));
                return c-h;
            };
            Model fm(face),fc(center(face));
            J d=fm.eval(face[8])[0],h=fm.eval(decimal(".95")*face[8])[1];h.d[8]=h.d[8]*decimal(".95");
            J dc=fc.eval(I(mid(face[8])))[0],hc=fc.eval(decimal(".95")*I(mid(face[8])))[1];
            J ff=expr(d,h);I centered=expr(dc,hc).v;
            for(int k=0;k<9;k++)centered=centered+ff.d[k]*(face[k]-I(mid(face[k])));
            I bound=cap(ff.v,centered);
            F.v=I(std::max(F.v.l,bound.l),F.v.h);
        }
    }
    result.good=F.v.l>=0&&(mode==0||cold.v.l>=decimal(".001").h);
    for(int k=0;k<9;k++)result.score[k]=mag(F.d[k])*(box[k].h-box[k].l);
    return result;
}
int cover(double radius,long maxnodes,int mode,I start_u=decimal(".02"),I end_u=I(3)){
    Box base=initial(radius);
    base[8]=mode==0?I(start_u.l,end_u.h):mode==1?I(decimal(".4").l,decimal(".8").h):mode==2?start_u:end_u;
    std::vector<Box> pending={base};long seen=0,passed=0,errors=0;int maxdepth=0;
    auto start=std::chrono::steady_clock::now();
    const double weights[9]={1.,.15,.7,.6,.4,1.,.4,.4,1.};
    while(!pending.empty()&&seen<maxnodes){
        Box b=pending.back();pending.pop_back();seen++;Check c;
        bool wide_time=b[8].h-b[8].l>.04;
        if(!wide_time){try{c=verify_cell(b,mode);}catch(const std::exception&e){errors++;c.reason=e.what();}}
        if(c.good){passed++;continue;}
        int split=-1;double best=-1;
        for(int k=0;k<9;k++)if(mid(b[k])>b[k].l&&mid(b[k])<b[k].h){
            double fallback=weights[k]*(b[k].h-b[k].l)/(base[k].h-base[k].l);
            double score=c.score[k]>0?c.score[k]:fallback;
            if(score>best){best=score;split=k;}
        }
        if(wide_time)split=8;
        if(split<0)throw std::runtime_error("unresolved point");
        double m=mid(b[split]);if(m==b[split].l||m==b[split].h)throw std::runtime_error("subdivision precision exhausted");
        Box upper=b;b[split].h=m;upper[split].l=m;pending.push_back(upper);pending.push_back(b);
        if(seen%1000==0){long double remaining=0;for(const Box& cell:pending){long double vol=1;for(int k=0;k<9;k++)if(base[k].h>base[k].l)vol*=((long double)cell[k].h-cell[k].l)/((long double)base[k].h-base[k].l);remaining+=vol;}std::cout<<"remaining_volume="<<(double)remaining<<" u=["<<b[8].l<<","<<b[8].h<<"] ";double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"progress cells="<<seen<<" passed="<<passed<<" pending="<<pending.size()<<" enclosure_errors="<<errors<<" seconds="<<sec<<std::endl;}
    }
    bool done=pending.empty();double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"complete\":"<<(done?"true":"false")<<",\"mode\":"<<mode<<",\"radius\":"<<radius<<",\"start_u\":["<<start_u.l<<","<<start_u.h<<"],\"end_u\":["<<end_u.l<<","<<end_u.h<<"],\"cells\":"<<seen<<",\"passed\":"<<passed<<",\"pending\":"<<pending.size()<<",\"enclosure_errors\":"<<errors<<",\"seconds\":"<<sec<<"}"<<std::endl;
    return done?0:2;
}
int main(int argc,char**argv){
    try{
        std::cout<<std::setprecision(17);
        if(argc>1&&std::string(argv[1])=="cover")return cover(std::stod(argv[2]),std::stol(argv[3]),std::stoi(argv[4]),argc>6?decimal(argv[5]):decimal(".02"),argc>6?decimal(argv[6]):I(3));
        double radius=argc>1?std::stod(argv[1]):0;
        Box box=initial(radius);if(argc>2)box[8]=decimal(argv[2]);
        Model family(box);auto test=family.eval(I(.4));
        std::cout<<"KL at t/tau=.4: cold ["<<test[0].v.l<<","<<test[0].v.h<<"] hot ["<<test[1].v.l<<","<<test[1].v.h<<"]\n";
        Margin m=margin(box);std::cout<<"time margin / tau ["<<m.value.l<<","<<m.value.h<<"]\n";
        for(int k=0;k<9;k++)std::cout<<"derivative "<<k<<" ["<<m.d[k].l<<","<<m.d[k].h<<"]\n";
        return 0;
    }catch(const std::exception&e){std::cerr<<"UNRESOLVED: "<<e.what()<<"\n";return 2;}
}
