#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

double Mp(double sm, double hm, double sf, double hf, double p){
  double ret = ((hm * sm + hf * sf) * p +
                ((1.0 - 3.0 * hm) * sm + (1.0 - 3.0 * hf) * sf) * p * p +
                ((2.0 * hm - 1.0) * sm + (2.0 * hf - 1.0) * sf) * p * p * p) / 2.0;
  return(ret);
}

std::vector<double> renew_mat(double sm, double hm, double p, double r, const std::vector<double>& mat){
  double a = hm*sm + (sm - 3.0*hm*sm) * p + (2.0*hm*sm - sm) * p * p;
  double b = -hm*sm*p + (2*hm*sm - sm) * p * p;
  double c = 1.0 - p;
  double d = p;

  std::vector<double> ret(4);
  ret.at(0) = mat.at(0) * (1.0 + a - c * r) + mat.at(1) * c * r;
  ret.at(1) = mat.at(0) * d * r + mat.at(1) * (1.0 + b - d * r);
  ret.at(2) = mat.at(2) * (1.0 + a - c * r) + mat.at(3) * c * r;
  ret.at(3) = mat.at(2) * d * r + mat.at(3) * (1.0 + b - d * r);

  return(ret);
}

int main(){
  double sm = 0.02;
  double hm = 0.5;
  double sf = -0.01;
  double hf = 0.5;

  double r = 0.02;

  std::ofstream ofs("eq_e4_positive.txt");
  ofs << "p\tEBB\tEbB\tEBb\tEbb" << std::endl;

  double p = 1.0 - 0.00001;

  std::vector<double> mat = {1.0, r/(r+(1-hm)*sm), 0.0, 0.0};

  while(p > 0.01){
    p -= Mp(sm, hm, sf, hf, p);
    mat = renew_mat(sm, hm, p, r, mat);
    ofs << p << "\t" << mat.at(0) << "\t" << mat.at(1) << "\t" << mat.at(2) << "\t" << mat.at(3) << std::endl;
  }
  return(0);
}
