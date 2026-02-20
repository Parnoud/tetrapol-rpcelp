#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <assert.h>
#include <math.h>
using namespace std;

double QBL[8]={0.3,0.45,0.60,0.74,0.86,0.98,1.14,1.23};
double sto_gain[32]={0,5,11,19,27,35,43,51,59,71,87,103,119,143,175,207,239,287,351,415,479,575,703,831,959,1151,1407,1663,1919,2303,2815,3583};

double ls[10][2]={
		  //  {-2.0643047, 0.38465041} ,
		  {-1.625, 0.38465041} ,
		  //{-0.7255134, 2.02982235} ,
		  {-0.7255134, 1.625},
		  {-1.3538182, 0.62888535} ,
		  {-0.4090979, 1.36354402} ,
		  {-0.7099015, 0.57737373} ,
		  {-0.4737691, 0.82117651} ,
		  {-0.8141674, 0.47146487} ,
		  {-0.5001262, 0.71935197} ,
		  {-0.7204825, 0.47273149} ,
		  {-0.3206912, 0.55081164}};

void sat(double &x,double M)
{
  if(x>M)x=M;
  if(x<-M)x=-M;
}

struct val_t {
  val_t(string n) { name=n;}
  string name;
  vector<int> vb;
  int len() const {
    return vb.size();
  }
  int get(const vector<bool> &v) const {
    int r=0;
    for(auto &it:vb) {
      r=r*2;
      if(v[it]) r++;
    }
    return r;
  }
};

vector<val_t*> vals,valslars;
val_t *valsltp[3][2];
val_t *valsex[3][3];

double ffdec=1.;

double declar(int u,int i,int nb) {
  assert(i>=0 && i<10);
  assert(u>=0);
  if(u>=(1<<nb))
    printf("error %d >= 2^%d\n",i,nb);
  double r=ffdec*ls[i][0]+(ffdec*ls[i][1]-ffdec*ls[i][0])*u/((1<<nb)-1);
  return r;
}

struct lar_t {
  vector<double> t;
};

double refl2lar_approx(double x)
{
  double a=fabs(x);
  double r;
  if(a<0.675) {
    r=a;
  } else if(a<.95) {
    r=2*a-0.675;
  } else
    r=8*a-6.375;
  if(x<0) r*=-1;
  return r;
}

double lar2refl_approx(double lar) {
  double a=fabs(lar);
  double r;
  if(a<0.675) {
    r=a;
  } else if(a<1.225) {
    r=0.5*a+0.3375;
  } else
    r=0.125*a+0.796875;

  if(lar<0) r*=-1;
  return r;
}

double lar2refl_eval(double lars) {
  lars=pow(10.,lars);
  return (lars-1)/(lars+1);
}

lar_t lar2refl(const lar_t &lar,double ex=1) {
  lar_t r;
  for(auto &it:lar.t) {
    if(ex)
      r.t.push_back(lar2refl_eval(it));
    else 
      r.t.push_back(lar2refl_approx(it));
  }
  return r;
}

lar_t operator*(double x,const lar_t &a) {
  lar_t r;
  r.t=a.t;
  assert(r.t.size()==10);
  for(int i=0;i<r.t.size();i++)
    r.t[i]*=x;
  return r;
}

lar_t operator+(const lar_t &b,const lar_t &a) {
  lar_t r;
  r.t=a.t;
  assert(a.t.size()==10 && b.t.size()==10);
  for(int i=0;i<r.t.size();i++)
    r.t[i]+=b.t[i];
  return r;
}

struct ltp_t {
  int ilag;
  int igain;
  int dec;
  int iexgain;
  int iph;
  int N;
  ltp_t(int n=0):N(n) {}

  int vp() const {
    int t[4]={3,4,4,6};
    return t[dec];
  }
  int Q() const {
    if(dec==0) return N==48?6:7;
    int t[4]={7,4,3,1};
    return t[dec];
    }
  int D() const {
    int t[4]={8,12,15,56};
    return t[dec];
  }
  int ph() const {
    int n;
    if(N==48) n=9;
    else n=10;
    return (iph>>(n-vp()))&((1<<(vp()))-1);
  }
  int sigs() const {
    int n;
    if(N==48) n=9;
    else n=10;
    return iph&((1<<(n-vp()))-1);
  }
  double gain() const {
    if(ilag==255) return 0;
    else return QBL[igain];
  }
  double exgain() const {
    return sto_gain[iexgain]/100000.;
  }
};

struct vfr_t {
  vector<bool> v;
  vfr_t(const vector<bool> &vv) {
    v=vv;
  }

  lar_t lar() const {
    lar_t l;
    assert(valslars.size()==10);
    for(int i=0;i<10;i++) {
      auto &it(valslars[i]);
      l.t.push_back(declar(it->get(v),i,it->len()));
    }
    return l;
  }

  vector<ltp_t> ltp() const {
    vector<ltp_t> r;
    for(int i=0;i<3;i++) {
      ltp_t u(i==1?48:56);
      u.ilag=valsltp[i][0]->get(v);
      u.igain=valsltp[i][1]->get(v);
      u.dec=valsex[i][0]->get(v);
      u.iexgain=valsex[i][1]->get(v);
      u.iph=valsex[i][2]->get(v);
      r.push_back(u);
    }
    
    return r;
  }
};

vector<lar_t> interpol(const lar_t &a, const lar_t &old)
{
  vector<lar_t> v;
  v.push_back(0.875*old+0.125*a);
  v.push_back(0.500*old+0.500*a);
  v.push_back(0.125*old+0.875*a);
  return v;
}

struct rpcelp_decode_t {
  double tab[160*3];
  double tab2[160*3];

  void clear() {
    memset(tab,0,sizeof(tab));
    memset(tab2,0,sizeof(tab));
    old.t.resize(0);
    old.t.resize(10);
    memset(v,0,sizeof(v));
  }
  
  rpcelp_decode_t() {
    clear();
  }
  
  lar_t correctlar(const lar_t &lar,double x) {
    lar_t r=lar;
    for(int i=1;i<=10;i++)
      r.t[i-1]/=1.001*pow(x,i);
    return r;
  }

  const double *decode(const vector<bool> &v)
  {
    vfr_t vfr(v);
    return s(vfr);
  }

  const double *decode(const unsigned char *bf) {
    vector<bool> v;
    for(int i=0;i<160;i++)
      v.push_back((bf[i/8]>>(i%8))&1);
    return decode(v);
  }

  ltp_t altp;
  lar_t old;
  
  const double *s(const vfr_t &vfr) {
    memmove(tab,tab+160,sizeof(double)*160*2);
    memmove(tab2,tab2+160,sizeof(double)*160*2);
    memset(tab+320,0,sizeof(double)*160);
    memset(tab2+320,0,sizeof(double)*160);

    auto lar=vfr.lar();
    auto ltp=vfr.ltp();
    auto lars=interpol(lar,old);
    old=lar;

    altp=ltp[0];
    lt(320,56,altp);
    ex(320,56,altp);
    
    altp=ltp[1];
    lt(320+56,48,altp);
    ex(320+56,48,altp);

    altp=ltp[2];
    lt(320+56+48,56,altp);
    ex(320+56+48,56,altp);
    
    st(320,56,lars[0]);
    st(320+56,48,lars[1]);
    st(320+56+48,56,lars[2]);
    
    return tab2+160*2;
  }
  
  double v[11];

  double interm13(double *t) 
  {
    static const double c[32] = {
				 -49,    66,   -96,   142,  -207,   294,  -407,   553,  -739,
				 981, -1304,  1758, -2452,  3688, -6669, 27072, 13496, -5287,  3179,
				 -2182,  1587, -1185,   893,  -672,   500,  -366,   263,  -183,   125,
				 -84,    59,   -47 };
    double s=0;
    for(int i=0;i<32;i++)
      s+=t[i-15]*c[i]/32768;
    return s;
  }

  double inter13(double *t)
  {
    static const double c[32] = {
				 -47,    59,   -84,   125,  -183,   263,  -366,   500,  -672,   893,
				 -1185,  1587, -2182,  3179, -5287, 13496, 27072, -6669,  3688, -2452,
				 1758, -1304,   981,  -739,   553,  -407,   294,  -207,   142,   -96,
				 66,   -49};
    
    double s=0;
    for(int i=0;i<32;i++)
      s+=t[i-16]*c[i]/32768;
    return s;
  }

  void ex(int st,int len,const ltp_t &ltp) { //excitation
    int D=ltp.D();
    int ph=ltp.ph();
    ph=ph%D;
    int s=ltp.sigs();
    int Q=ltp.Q();
    auto g1=ltp.gain();
    auto g2=ltp.exgain();
    for(int i=0;i<len;i++)
      tab[st+i]*=g1;
    for(int i=0;i<Q && ph<len;i++) {
      if( (s>>(Q-1-i)) & 1)
	tab[st+ph]+=+g2;
      else 
	tab[st+ph]+=-g2;
      ph+=D;
    }
  }

  void lt(int st,int len,const ltp_t &ltp) { //long time prediction
    auto ilag=ltp.ilag;
    auto g=ltp.gain();
    int fr=(ilag+60)%3;
    int T0=(ilag+61)/3;
    for(int i=st;i<st+len;i++) {
      int j=i-T0;
      if(fr==0)
	tab[i]+=tab[j];
      else if(fr==1)
         tab[i]+=inter13(&tab[i-T0]);
      else {
	assert(fr==2);
	tab[i]+=interm13(&tab[i-T0]);
      }
    }
  }

  double ffb=1.05; 
  
  void st(int st,int len,const lar_t &ll_) { //short time prediction
    auto lar=lar2refl(ll_,0);
    lar=correctlar(lar,ffb); // avoid overflow due to LAR cooef errors
    
    for(int i=st;i<st+len;i++) {
      double sri=tab[i];

      for(int t=1;t<=10;t++) {
	sri-=lar.t[10-t]*v[10-t];
	v[11-t]=v[10-t]+lar.t[10-t]*sri;
	sat(v[11-t],1000);
      }
      tab2[i]=sri;
      v[0]=sri;
      sat(v[0],1000);
    }
  }
};

const char *bits_raw[]={
			   "LAR01", "1", "0", "23", "22", "21",
			   "LAR02", "2", "27", "26", "25", "24",
			   "LAR03", "3", "30", "29", "28",
			   "LAR04", "4", "31", /***/ "33", "32",
			   "LAR05", "5",  "36", "35", "34",
			   "LAR06", "39", "38", "37",/***/
			   "LAR07", "42", "41", "40",
			   "LAR08", "45", "44", "43",
			   "LAR09", "47", "46",/***/ "48",
			   "LAR10", "51", "50", "49",
			   "LTP_lag1", "6", "58", "57", "56", "55", "54", "53", "52", /***/    // 52-58 56-sp
			   "LTP_gain1", "7", "60", "59",
			   "E_gain1", "8", "63", "62", "61", /***/ "64",
			   "E_dec1", "10", "9",
			   "E_sph1",  "74", "73", "72", "71", "70", "69", "68", "67", "66", "65",  //68sp

			   "LTP_lag2", "11", "79", "78", "77", "76", "75", /***/ "81", "80",   // 75-81  77sp
			   "LTP_gain2", "12", "83", "82",
			   "E_gain2", "13", "87", "86", "85", "84", /***/
			   "E_dec2", "15", "14",
			   "E_sph2", "95", "94", "93", "92", "91", "90", "89", "88", /***/ "96", //89sp

			   "LTP_lag3", "16", "103", "102", "101", "100", "99", "98", "97", /***/ // 97-103 101sp
			   "LTP_gain3", "17", "105", "104",
			   "E_gain3", "18", "109", "108", "107", "106",
			   "E_dec3", "20", "19",
			   "E_sph3", "111", "110", "119", "118", "117", "116", "115", "114", "113", "112", //115sp
	    NULL
};

map<string,list<int> > bits;

val_t *val_search(const string &a)
{
  for(auto &it:vals)
    if(it->name==a) return it;
  cout<<a<<endl;
  assert(0);
  return NULL;
}


void init() {
  for(auto &it:vals)
    it->vb.clear();

  string rr="";
  for(int i=0;bits_raw[i];i++) {
    if(bits_raw[i][0]>='0' && bits_raw[i][0]<='9') {
      int u=(atoi(bits_raw[i]));
      for(auto &it:vals)
	if(it->name==rr)
	  it->vb.push_back(u); 
    } else {
      rr=(bits_raw[i]);
      vals.push_back(new val_t(rr));
    }
  }

  
  valslars.clear();
  for(int i=0;i<10;i++) {
    char bf[100];
    snprintf(bf,99,"LAR%02d",i+1);
    auto t=val_search(bf);
    assert(t);
    valslars.push_back(t);
  }

  valsltp[0][0]=val_search("LTP_lag1");
  valsltp[0][1]=val_search("LTP_gain1");
  valsltp[1][0]=val_search("LTP_lag2");
  valsltp[1][1]=val_search("LTP_gain2");
  valsltp[2][0]=val_search("LTP_lag3");
  valsltp[2][1]=val_search("LTP_gain3");

  valsex[0][0]=val_search("E_dec1");
  valsex[0][1]=val_search("E_gain1");
  valsex[0][2]=val_search("E_sph1");
  valsex[1][0]=val_search("E_dec2");
  valsex[1][1]=val_search("E_gain2");
  valsex[1][2]=val_search("E_sph2");
  valsex[2][0]=val_search("E_dec3");
  valsex[2][1]=val_search("E_gain3");
  valsex[2][2]=val_search("E_sph3");
  
  rr="";
  for(int i=0;bits_raw[i];i++) {
    if(bits_raw[i][0]>='0' && bits_raw[i][0]<='9')
      bits[rr].push_back(atoi(bits_raw[i]));
    else
      rr=(bits_raw[i]);
  }

}

int h2i(char c) {
  if(c>='0' && c<='9') return c-'0';
  if(c>='a' && c<='f') return c-'a'+10;
  return -1;
}


string a(char *c)
{
  char bf[15];
  memset(bf,0,sizeof(bf));
  for(int u=0;u<30;u++) {
    int z=h2i(c[u]);
    if(z<0) return "";
    bf[u/2]|=(z<<(4*(1-(u%2))));
  }
  return bf;
}


int main(int ac, char **av) {
  double volume=100000;
  if(ac>1) volume=atof(av[1]);
  init();
  char bf[1000];
  rpcelp_decode_t decode;
  while(fgets(bf,9999,stdin)) {
    const double *r;

#if 0 //json file
    for(int i=0;bf[i];i++) {
      auto x=a(bf+i);
      if(x.empty()) continue;
      r=decode.decode((unsigned char *)x.c_str());
      break;
    }
#endif

#if 1 //text file with one  binary frame by line
    vector<bool> v;
    for(int i=0;i<120;i++)
      v.push_back(bf[i]-'0');
    r=decode.decode(v);
#endif

    short s[160];
    for(int i=0;i<160;i++)
      s[i]=r[i]*volume;
    fwrite(s,2,160,stdout);
  }
  return 0;
}
