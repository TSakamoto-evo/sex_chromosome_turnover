#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<cmath>

int sign(double a){
  if(a >= 0){
    return(1);
  }else{
    return(-1);
  }
}

double cubic(double a3, double a2, double a1, double a0){
  double b2 = a2 / a3;
  double b1 = a1 / a3;
  double b0 = a0 / a3;

  double p = b1 - std::pow(b2, 2.0) / 3.0;
  double q = b0 - 1.0/3.0 * b1 * b2 + 2.0/27.0 * std::pow(b2, 3.0);
  double det = std::pow(q/2.0, 2.0) + std::pow(p/3.0, 3.0);

  double x, y;
  if(det > 0){
    y = std::pow( std::abs(-q/2.0 + std::sqrt(det)), 1.0/3.0 ) *
          sign(-q/2.0 + std::sqrt(det)) +
          std::pow( std::abs(-q/2.0 - std::sqrt(det)), 1.0/3.0 ) *
          sign(-q/2.0 - std::sqrt(det));
    x = y - b2 / 3.0;
  }else{
    double tmp_a = std::sqrt(-p/3.0);
    double tmp_b = -q / std::pow(tmp_a, 2.0);

    y = 2.0 * tmp_a * std::cos(1.0/3.0 * std::acos(tmp_b / 2.0 / tmp_a));
    x = y - b2 / 3.0;
  }
  return (x);
}

std::vector<double> bt_app(double sm, double hm, double r, double p, double u, double v){
  double a = hm*sm + (sm - 3.0*hm*sm) * p + (2.0*hm*sm - sm) * p * p;
  double b = -hm*sm*p + (2*hm*sm - sm) * p * p;
  double c = 1.0 - p;
  double d = p;

  double tmp_a = 2*(a-c*r-u);
  double tmp_b = 2*(c*r+u);
  double tmp_c = 2*(d*r+v);
  double tmp_d = 2*(b-d*r-v);

  double ret1, ret0;
  if(tmp_a + tmp_d <= 0 && tmp_a*tmp_d - tmp_b*tmp_c >=0){
    ret1 = 0.0;
    ret0 = 0.0;
  }else if(tmp_b == 0){
    ret1 = tmp_a;
    ret0 = 0.0;
  }else{
    ret1 = cubic(1.0, -2*tmp_a, tmp_a*tmp_a-tmp_b*tmp_d, tmp_a*tmp_b*tmp_d - tmp_b*tmp_b*tmp_c);
    ret0 = (ret1 * ret1 - tmp_a * ret1) / tmp_b;
  }
  std::vector<double> ret = {ret1, ret0};
  return(ret);
}

int main(){
  double sm = 0.02;
  double hm = 0.5;

  double r = 0.02;
  double u = 0.000001;
  double v = 0.000001;

  int seps = 10000;

  std::ofstream ofs("const.txt");
  ofs << "p\tphi_B\tphi_b" << std::endl;

  for(int i = 0; i <= seps; i++){
    double p = i * 1.0 / seps;
    std::vector<double> phi;
    phi = bt_app(sm, hm, r, p, u, v);

    ofs << p << "\t" << phi.at(0) << "\t" << phi.at(1) << std::endl;
  }

  return(0);
}
