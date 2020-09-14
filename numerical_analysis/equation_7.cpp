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

double Mp(double sm, double hm, double sf, double hf, double p, double u, double v){
  double ret = ((hm * sm + hf * sf) * p +
                ((1.0 - 3.0 * hm) * sm + (1.0 - 3.0 * hf) * sf) * p * p +
                ((2.0 * hm - 1.0) * sm + (2.0 * hf - 1.0) * sf) * p * p * p) / 2.0 +
                 - u * p + v * (1.0 - p);
  return(ret);
}

double fpt1(double sm, double hm,
            double r, double p, double phi1, double phi0, double u){
  double a = hm*sm + (sm - 3.0*hm*sm) * p + (2.0*hm*sm - sm) * p * p;
  double c = 1.0 - p;

  double ret = (phi1 * phi1 / 2.0 - (a - c*r - u) * phi1 - (c*r + u) * phi0);
  return(ret);
}

double fpt2(double sm, double hm,
            double r, double p, double phi1, double phi0, double v){
  double b = -hm*sm*p + (2*hm*sm - sm) * p * p;
  double d = p;

  double ret = (phi0 * phi0 / 2.0 - (b - d*r - v) * phi0 - (d*r + v) * phi1);
  return(ret);
}

int main(){
  double sm = 0.02;
  double hm = 0.5;
  double sf = -0.01;
  double hf = 0.5;

  double r = 0.02;
  double u = 0.000001;
  double v = 0.000001;

  double eps = 1.0 / 10000;
  double h = 0.5;

  std::ofstream ofs("change.txt");
  ofs << "p\tphi_B\tphi_b" << std::endl;

  double equ_p = 0.99960016;
  std::vector<double> phi_equ;
  phi_equ = bt_app(sm, hm, r, equ_p, u, v);

  {
    // p < equ_p
    double phi1 = phi_equ.at(0);
    double phi0 = phi_equ.at(1);
    double p = equ_p - eps;
    if(p > 0.0){
      ofs << p << "\t" << phi1 << "\t" << phi0 << std::endl;
    }

    // Runge-Kutta method
    while(p > eps){
      double k11 = -h * fpt1(sm, hm, r, p, phi1, phi0, u);
      double k01 = -h * fpt2(sm, hm, r, p, phi1, phi0, v);
      double kp1 = -h * Mp(sm, hm, sf, hf, p, u, v);

      double k12 = -h * fpt1(sm, hm, r, p + kp1 / 2.0, phi1 + k11 / 2.0, phi0 + k01 / 2.0, u);
      double k02 = -h * fpt2(sm, hm, r, p + kp1 / 2.0, phi1 + k11 / 2.0, phi0 + k01 / 2.0, v);
      double kp2 = -h * Mp(sm, hm, sf, hf, p + kp1 / 2.0, u, v);

      double k13 = -h * fpt1(sm, hm, r, p + kp2 / 2.0, phi1 + k12 / 2.0, phi0 + k02 / 2.0, u);
      double k03 = -h * fpt2(sm, hm, r, p + kp2 / 2.0, phi1 + k12 / 2.0, phi0 + k02 / 2.0, v);
      double kp3 = -h * Mp(sm, hm, sf, hf, p + kp2 / 2.0, u, v);

      double k14 = -h * fpt1(sm, hm, r, p + kp3, phi1 + k13, phi0 + k03, u);
      double k04 = -h * fpt2(sm, hm, r, p + kp3, phi1 + k13, phi0 + k03, v);
      double kp4 = -h * Mp(sm, hm, sf, hf, p + kp3, u, v);

      phi1 = phi1 + k11 / 6.0 + k12 / 3.0 + k13 / 3.0 + k14 / 6.0;
      phi0 = phi0 + k01 / 6.0 + k02 / 3.0 + k03 / 3.0 + k04 / 6.0;
      p = p + kp1 / 6.0 + kp2 / 3.0 + kp3 / 3.0 + kp4 / 6.0;
      ofs << p << "\t" << phi1 << "\t" << phi0 << std::endl;
    }
  }

  {
    // p > equ_p
    double phi1 = phi_equ.at(0);
    double phi0 = phi_equ.at(1);
    double p = equ_p + eps;
    if(p < 1.0){
      ofs << p << "\t" << phi1 << "\t" << phi0 << std::endl;
    }

    // Runge-Kutta method
    while(p < 1.0 - eps){
      double k11 = -h * fpt1(sm, hm, r, p, phi1, phi0, u);
      double k01 = -h * fpt2(sm, hm, r, p, phi1, phi0, v);
      double kp1 = -h * Mp(sm, hm, sf, hf, p, u, v);

      double k12 = -h * fpt1(sm, hm, r, p + kp1 / 2.0, phi1 + k11 / 2.0, phi0 + k01 / 2.0, u);
      double k02 = -h * fpt2(sm, hm, r, p + kp1 / 2.0, phi1 + k11 / 2.0, phi0 + k01 / 2.0, v);
      double kp2 = -h * Mp(sm, hm, sf, hf, p + kp1 / 2.0, u, v);

      double k13 = -h * fpt1(sm, hm, r, p + kp2 / 2.0, phi1 + k12 / 2.0, phi0 + k02 / 2.0, u);
      double k03 = -h * fpt2(sm, hm, r, p + kp2 / 2.0, phi1 + k12 / 2.0, phi0 + k02 / 2.0, v);
      double kp3 = -h * Mp(sm, hm, sf, hf, p + kp2 / 2.0, u, v);

      double k14 = -h * fpt1(sm, hm, r, p + kp3, phi1 + k13, phi0 + k03, u);
      double k04 = -h * fpt2(sm, hm, r, p + kp3, phi1 + k13, phi0 + k03, v);
      double kp4 = -h * Mp(sm, hm, sf, hf, p + kp3, u, v);

      phi1 = phi1 + k11 / 6.0 + k12 / 3.0 + k13 / 3.0 + k14 / 6.0;
      phi0 = phi0 + k01 / 6.0 + k02 / 3.0 + k03 / 3.0 + k04 / 6.0;
      p = p + kp1 / 6.0 + kp2 / 3.0 + kp3 / 3.0 + kp4 / 6.0;
      ofs << p << "\t" << phi1 << "\t" << phi0 << std::endl;
    }
  }

  return(0);
}
