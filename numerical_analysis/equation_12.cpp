#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<iomanip>

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
                ((2.0 * hm - 1.0) * sm + (2.0 * hf - 1.0) * sf) * p * p * p) / 2.0
                - u*p + v*(1-p);
  return(ret);
}

std::vector<double> set_for_coef(int seps){
  std::vector<double> ret = {1.0};
  for(int i = 1; i < seps; i++){
    ret.push_back(1.0);
  }
  ret.push_back(1.0);
  return(ret);
}

std::vector<double> set_for_moments(int seps, int degree){
  std::vector<double> ret = {0.0};
  for(int i = 1; i < seps; i++){
    double p = i * 1.0 / seps;
    ret.push_back(std::pow(p, degree));
  }
  ret.push_back(1.0);
  return(ret);
}

std::vector<double> set_for_constant_p(double s, double h, double r, int seps, double u, double v){
  std::vector<double> ret = {};
  for(int i = 0; i <= seps; i++){
    double p = i * 1.0 / seps;
    std::vector<double> phis = bt_app(s, h, r, p, u, v);
    ret.push_back(p * phis.at(0) + (1-p) * phis.at(1));
  }
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

std::vector<double> set_for_changing_p(double sm, double hm, double sf, double hf, double r, int seps, double u, double v, double equ_p){
  int equ_low_index = std::floor(equ_p * seps);
  double h = 0.1;
  std::vector<double> phis = bt_app(sm, hm, r, equ_p, u, v);
  double phi1 = phis.at(0);
  double phi0 = phis.at(1);
  std::vector<double> ret(seps + 1, 0.0);

  //register previous step
  double regi_p;
  int regi_floor;
  double regi_phi1;
  double regi_phi0;

  double p = 1.0 * equ_low_index / seps;
  ret.at(equ_low_index) = p * phi1 + (1.0 - p) * phi0;

  regi_p = p;
  regi_floor = equ_low_index - 1;
  regi_phi1 = phi1;
  regi_phi0 = phi0;

  if(equ_low_index > 0){
    p = 1.0 * (equ_low_index - 1) / seps;
    ret.at(equ_low_index - 1) = p * phi1 + (1.0 - p) * phi0;

    regi_p = p;
    regi_floor = equ_low_index - 2;
  }

  while(p > 1.0 / seps){
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

    if(std::floor(equ_p * seps) != regi_floor){
      for(int i = regi_floor; i > std::floor(p * seps); i--){
        double focal_p = 1.0 * i / seps;
        double focal_phi1 = (regi_phi1 - phi1) / (regi_p - p) * (focal_p - p) + phi1;
        double focal_phi0 = (regi_phi0 - phi0) / (regi_p - p) * (focal_p - p) + phi0;
        ret.at(i) = focal_p * focal_phi1 + (1.0 - focal_p) * focal_phi0;
      }
    }

    regi_p = p;
    regi_floor = std::floor(p * seps);
    regi_phi1 = phi1;
    regi_phi0 = phi0;
  }

  if(equ_low_index != 0){
    ret.at(0) = ret.at(1);
  }

  phi1 = phis.at(0);
  phi0 = phis.at(1);

  p = 1.0 * (equ_low_index + 1) / seps;
  ret.at(equ_low_index + 1) = p * phi1 + (1.0 - p) * phi0;

  regi_p = p;
  regi_floor = equ_low_index + 1;
  regi_phi1 = phi1;
  regi_phi0 = phi0;

  if(equ_low_index < seps - 1){
    p = 1.0 * (equ_low_index + 2) / seps;
    ret.at(equ_low_index + 2) = p * phi1 + (1.0 - p) * phi0;

    regi_p = p;
    regi_floor = equ_low_index + 2;
  }

  while(p < 1.0 - 1.0 / seps){
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

    if(std::floor(equ_p * seps) != regi_floor){
      for(int i = regi_floor; i < std::floor(p * seps); i++){
        double focal_p = 1.0 * (i + 1) / seps;
        double focal_phi1 = (regi_phi1 - phi1) / (regi_p - p) * (focal_p - p) + phi1;
        double focal_phi0 = (regi_phi0 - phi0) / (regi_p - p) * (focal_p - p) + phi0;
        ret.at(i + 1) = focal_p * focal_phi1 + (1.0 - focal_p) * focal_phi0;
      }
    }

    regi_p = p;
    regi_floor = std::floor(p * seps);
    regi_phi1 = phi1;
    regi_phi0 = phi0;
  }
  if(equ_low_index != seps){
    ret.at(seps) = ret.at(seps - 1);
  }

  return(ret);
}

double return_index(double x, int seps, const std::vector<double>& int_values){
  int low_index = std::floor(x * seps);
  if(low_index < seps){
    double reminder = x * seps - low_index;
    return(int_values.at(low_index) * (1-reminder) + int_values.at(low_index+1) * reminder);
  }else if(low_index == seps){
    return(int_values.at(seps));
  }else{
    std::cerr << "OUT OF RANGE!" << '\n';
    return(0);
  }
}

double log_function(double sm, double hm, double sf, double hf, double pop_size, double y,
                    int seps, const std::vector<double>& int_values){
  double x = (y + 1.0) / 2.0;
  double sum = 0.0;
  if(std::abs(sm + sf) < 0.0000001 && std::abs(hm - hf) < 0.0000001){
    double add = (-(1-2*hf)*(1-2*hf)*x*x*x*x - 4*hf*(1-2*hf)*x*x*x -4*hf*hf*x*x)*pop_size*sf*sf;
    sum += add;
  }else{
    double add = 2*pop_size*(hm*sm+hf*sf)*x + pop_size*((1-2*hm)*sm+(1-2*hf)*sf)*x*x;
    sum += add;
  }
  sum += std::log(return_index(x, seps, int_values));
  return(sum);
}

double integrate_D(double sm, double hm, double sf, double hf, double u, double v, double pop_size,
                    int seps, const std::vector<double>& int_values, int n, double h, double& base){
    double sum = 0.0;
    double pi = std::acos(-1);
    double log_const = std::log(pi * h) - (u + v) * std::log(2.0);

    double max_log = 0.0;

    for(int k = -n; k <= n; k++){
        double x = std::tanh(pi * std::sinh(k * h) / 2.0);
        double log_value = log_function(sm, hm, sf, hf, pop_size, x, seps, int_values) +
                            (v - u) * pi * std::sinh(k * h) / 2.0 +
                            std::log(std::cosh(k * h));

        double log_value2 = (u + v) * std::log(std::cosh(pi * std::sinh(k * h) / 2.0));

        // cosh(x) is approximated by exponential if |x| is large
        if(std::sinh(k * h) < -100.0){
            log_value2 = (u + v) * (-pi * std::sinh(k * h) / 2.0 - std::log(2.0));
        }
        if(std::sinh(k * h) > 100.0){
            log_value2 = (u + v) * (pi * std::sinh(k * h) / 2.0 - std::log(2.0));
        }

        if(k == -n){
          max_log = log_const + log_value - log_value2;
        }else if(log_const + log_value - log_value2 > max_log){
          max_log = log_const + log_value - log_value2;
        }
        sum += std::exp(log_const + log_value - log_value2 - base);
    }
    base = max_log;
    return(sum);
}

int main(){
  double sm = 0.02;
  double hm = 0.5;
  double sf = -0.02;
  double hf = 0.5;

  double u = 0.000001;
  double v = 0.000001;

  double r = 0.000;
  int pop_size = 100000;

  double equ_p = 0.5;

  int seps = 10000;
  int i_max = 1000;
  std::vector<double> int_values;

  std::ofstream ofs("unconditional.txt");
  ofs << "r\tphi\t1/C" << std::endl;

  for(int i = -i_max; i <= 0; i++){
    r = 0.5 * pow(10.0, 4.0 * i / i_max);

    int_values = set_for_coef(seps);

    /* calculate 1/C in Equation 10 (or Equation 11) */
    /* 1/C is represented by coef * exp(coef_base) */
    /* integration is done by using double exponential formulas */
    double coef_base = 0.0;
    integrate_D(sm, hm, sf, hf, 4.0*u*pop_size, 4.0*v*pop_size, pop_size, seps, int_values, 2000, 0.01, coef_base);
    double coef = integrate_D(sm, hm, sf, hf, 4.0*u*pop_size, 4.0*v*pop_size, pop_size, seps, int_values, 2000, 0.01, coef_base);

    /* calculate p * phi_B(p) + (1-p) * phi_b(p) */
    /* choose Equation 3 or 7 */
    int_values = set_for_constant_p(sm, hm, r, seps, u, v);
    //int_values = set_for_changing_p(sm, hm, sf, hf, r, seps, u, v, equ_p);

    /* calculate phi/C (Equation 12) */
    /* phi/C is represented by ret * exp(ret_base) */
    /* integration is done by using double exponential formulas */
    double ret_base = 0.0;
    integrate_D(sm, hm, sf, hf, 4.0*u*pop_size, 4.0*v*pop_size, pop_size, seps, int_values, 2000, 0.01, ret_base);
    double ret = integrate_D(sm, hm, sf, hf, 4.0*u*pop_size, 4.0*v*pop_size, pop_size, seps, int_values, 2000, 0.01, ret_base);

    ofs << r << "\t" << ret / coef * std::exp(ret_base - coef_base)
      << "\t" << coef * std::exp(coef_base) << std::endl;
  }
  return(0);
}
