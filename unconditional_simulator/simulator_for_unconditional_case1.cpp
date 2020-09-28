#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<cmath>
#include<string>

class Genotype_def{
public:
  std::vector<bool> my;
  std::vector<bool> my_neo;
  std::vector<bool> ma;
  std::vector<bool> fy;
  std::vector<bool> fy_neo;
  std::vector<bool> fa;
  std::vector<bool> female;

  Genotype_def(){
    for(int i = 0; i < 64; i++){
      my.push_back( (i / 32) % 2 );
      my_neo.push_back( (i / 16) % 2 );
      ma.push_back( (i /  8) % 2 );
      fy.push_back( (i /  4) % 2 );
      fy_neo.push_back( (i /  2) % 2 );
      fa.push_back( (i /  1) % 2 );

      if(my.at(i) + my_neo.at(i) + fy.at(i) + fy_neo.at(i) > 0){
        female.push_back(0);
      }else{
        female.push_back(1);
      }
    }
  }
};

std::vector<double> male_gamete(const std::vector<int>& freqs,
                                const Genotype_def geno, const double sm,
                                const double hm, const double rm, const double u, const double v){
  //genotype frequencies among male
  std::vector<double> sex_freqs(64);
  double sum = 0.0;
  for(int i = 0; i < 64; i++){
    if(geno.female.at(i) == 0){
      double fitness;
      if(geno.ma.at(i) + geno.fa.at(i) == 2){
        fitness = 1.0 + sm;
      }else if(geno.ma.at(i) + geno.fa.at(i) == 1){
        fitness = 1.0 + hm * sm;
      }else{
        fitness = 1.0;
      }

      sex_freqs.at(i) += freqs.at(i) * fitness;
      sum += freqs.at(i) * fitness;
    }
  }


  for(int i = 0; i < 64; i++){
    sex_freqs.at(i) /= sum;
  }

  std::vector<double> gametes(8);
  std::vector<double> gametes_mut(8);

  for(int i = 0; i < 64; i++){
    //fff
    gametes.at(geno.fy.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //ffm
    gametes.at(geno.fy.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //fmf
    gametes.at(geno.fy.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //fmm
    gametes.at(geno.fy.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //mff
    gametes.at(geno.my.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //mfm
    gametes.at(geno.my.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //mmf
    gametes.at(geno.my.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //mmm
    gametes.at(geno.my.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
  }

  gametes_mut.at(0) = (1-v) * gametes.at(0) + u * gametes.at(1);
  gametes_mut.at(1) = (1-u) * gametes.at(1) + v * gametes.at(0);
  gametes_mut.at(2) = (1-v) * gametes.at(2) + u * gametes.at(3);
  gametes_mut.at(3) = (1-u) * gametes.at(3) + v * gametes.at(2);
  gametes_mut.at(4) = (1-v) * gametes.at(4) + u * gametes.at(5);
  gametes_mut.at(5) = (1-u) * gametes.at(5) + v * gametes.at(4);
  gametes_mut.at(6) = (1-v) * gametes.at(6) + u * gametes.at(7);
  gametes_mut.at(7) = (1-u) * gametes.at(7) + v * gametes.at(6);

  return(gametes_mut);
}

std::vector<double> female_gamete(const std::vector<int>& freqs,
                                const Genotype_def geno, const double sf,
                                const double hf, const double rf, const double u, const double v){
  //genotype frequencies among female
  std::vector<double> sex_freqs(64);
  double sum = 0.0;
  for(int i = 0; i < 64; i++){
    if(geno.female.at(i) == 1){
      double fitness;
      if(geno.ma.at(i) + geno.fa.at(i) == 2){
        fitness = 1.0 + sf;
      }else if(geno.ma.at(i) + geno.fa.at(i) == 1){
        fitness = 1.0 + hf * sf;
      }else{
        fitness = 1.0;
      }

      sex_freqs.at(i) += freqs.at(i) * fitness;
      sum += freqs.at(i) * fitness;
    }
  }
  for(int i = 0; i < 64; i++){
    sex_freqs.at(i) /= sum;
  }

  std::vector<double> gametes(8);
  std::vector<double> gametes_mut(8);

  for(int i = 0; i < 64; i++){
    //fff
    gametes.at(geno.fy.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //ffm
    gametes.at(geno.fy.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //fmf
    gametes.at(geno.fy.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //fmm
    gametes.at(geno.fy.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //mff
    gametes.at(geno.my.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //mfm
    gametes.at(geno.my.at(i) * 4 + geno.fy_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //mmf
    gametes.at(geno.my.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //mmm
    gametes.at(geno.my.at(i) * 4 + geno.my_neo.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
  }

  gametes_mut.at(0) = (1-v) * gametes.at(0) + u * gametes.at(1);
  gametes_mut.at(1) = (1-u) * gametes.at(1) + v * gametes.at(0);
  gametes_mut.at(2) = (1-v) * gametes.at(2) + u * gametes.at(3);
  gametes_mut.at(3) = (1-u) * gametes.at(3) + v * gametes.at(2);
  gametes_mut.at(4) = (1-v) * gametes.at(4) + u * gametes.at(5);
  gametes_mut.at(5) = (1-u) * gametes.at(5) + v * gametes.at(4);
  gametes_mut.at(6) = (1-v) * gametes.at(6) + u * gametes.at(7);
  gametes_mut.at(7) = (1-u) * gametes.at(7) + v * gametes.at(6);

  return(gametes_mut);
}

std::vector<int> next_gen(const std::vector<double>& male_gamete,
                            const std::vector<double>& female_gamete,
                            const int n, std::mt19937& mt){
  std::vector<double> deterministic(64);
  for(int i = 0; i < 64; i++){
    int male = i / 8;
    int female = i % 8;

    deterministic.at(i) = male_gamete.at(male) * female_gamete.at(female);
  }
  std::vector<int> next_generarion(64);
  int n_remain = n;

  for(int i = 0; i < 63; i++){
    double det_remain = 0.0;
    for(int j = i; j < 64; j++){
      det_remain += deterministic.at(j);
    }

    std::binomial_distribution<> choose(n_remain, deterministic.at(i) / det_remain);
    next_generarion.at(i) = choose(mt);
    n_remain -= next_generarion.at(i);
  }
  next_generarion.at(63) = n_remain;
  return(next_generarion);
}

void initialize_freqs_xy(std::vector<int>& freqs, const int n, std::mt19937& mt,
                          const double m00, const double m01, const double m10, const double m11,
                          const double f00, const double f01, const double f10, const double f11){
  //XY system
  std::vector<double> m_gametes = {m00, m01, 0.0, 0.0, m10, m11, 0.0, 0.0};
  std::vector<double> f_gametes = {f00, f01, 0.0, 0.0, f10, f11, 0.0, 0.0};

  freqs = next_gen(m_gametes, f_gametes, n, mt);
}

void introduce_y(std::vector<int>& freqs, const int n, std::mt19937& mt){
  std::uniform_int_distribution<> mut(0, n - 1);
  std::bernoulli_distribution p(0.5);

  int mutant = mut(mt);
  int k = 0;
  while(mutant >= freqs.at(k)){
    mutant -= freqs.at(k);
    k++;
  }
  freqs.at(k)--;
  if(p(mt)){
    freqs.at(k + 16)++;
  }else{
    freqs.at(k + 2)++;
  }
}

std::vector<int> return_yy(const std::vector<int>& freqs, const Genotype_def geno){
  int num_y = 0;
  int num_y_neo = 0;
  for(int i = 0; i < 64; i++){
    num_y += freqs.at(i) * (geno.my.at(i) + geno.fy.at(i));
    num_y_neo += freqs.at(i) * (geno.my_neo.at(i) + geno.fy_neo.at(i));
  }
  std::vector<int> ret = {num_y, num_y_neo};
  return(ret);
}

int main(int argc, char *argv[]){
  //invasion of ZW system
  //parameters
  int n;
  double sf;
  double hf;
  double sm;
  double hm;

  double rf = 0.0;
  double rm;
  double u;
  double v;

  if(argc == 9){
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%lf", &sf);
    sscanf(argv[3], "%lf", &hf);
    sscanf(argv[4], "%lf", &sm);
    sscanf(argv[5], "%lf", &hm);
    sscanf(argv[6], "%lf", &rm);
    sscanf(argv[7], "%lf", &u);
    sscanf(argv[8], "%lf", &v);
  }else{
    std::cerr << "Number of variables not matched!" << std::endl;
    return(1);
  }
  std::cout << n << "\t" << sf << "\t" << hf << "\t" << sm << "\t" << hm << "\t" << rm << "\t" << u << "\t" << v << std::endl;

  int max_rep = 1000;

  std::random_device seed;
  std::mt19937 mt(seed());

  std::vector<int> freqs;
  Genotype_def geno;

  std::ofstream ofs("establishment_prob_seq.txt", std::ios::app);

  int reps = 0;
  int regi_y_neo_fix = 0;

  std::ifstream ifs("distribution1.txt");
  std::string str;

  if(ifs.fail()){
    std::cerr << "file input failed" << std::endl;
    return(1);
  }

  while(getline(ifs, str)){
    int index;
    double m00, m01, m10, m11;
    double f00, f01, f10, f11;

    sscanf(str.c_str(), "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &index, &m00, &m01, &m10, &m11, &f00, &f01, &f10, &f11);

    //initialize
    initialize_freqs_xy(freqs, n, mt, m00, m01, m10, m11, f00, f01, f10, f11);
    introduce_y(freqs, n, mt);

    for(int i = 0; 1; i++){
      std::vector<double> m_gamete = male_gamete(freqs, geno, sm, hm, rm, u, v);
      std::vector<double> f_gamete = female_gamete(freqs, geno, sf, hf, rf, u, v);
      freqs = next_gen(m_gamete, f_gamete, n, mt);

      std::vector<int> yy = return_yy(freqs, geno);

      if(yy.at(0) == 0 || yy.at(1) == 0){
        if(yy.at(0) == 0){
          regi_y_neo_fix++;
        }
        break;
      }
    }
    reps++;
    if((regi_y_neo_fix >= max_rep && reps % 10000 == 0) || reps == 2000000){
      break;
    }
  }
  ofs << n << "\t" << sf << "\t" << hf << "\t" << sm << "\t" << hm << "\t" << rm << "\t" << max_rep << "\t1\t";
  ofs << regi_y_neo_fix * 1.0 / reps << "\t" << regi_y_neo_fix << "\t" << reps << "\t" << u << "\t" << v << std::endl;

  return(0);
}
