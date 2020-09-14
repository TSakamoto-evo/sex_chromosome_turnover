#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<string>
#include<iomanip>

class Genotype_def{
public:
  std::vector<bool> my;
  std::vector<bool> ma;
  std::vector<bool> fy;
  std::vector<bool> fa;
  std::vector<bool> female;

  Genotype_def(){
    for(int i = 0; i < 16; i++){
      my.push_back( (i / 8) % 2 );
      ma.push_back( (i / 4) % 2 );
      fy.push_back( (i / 2) % 2 );
      fa.push_back( (i / 1) % 2 );

      if(my.at(i) + fy.at(i) > 0){
        female.push_back(0);
      }else{
        female.push_back(1);
      }
    }
  }
};

std::vector<double> male_gamate(const std::vector<int>& freqs,
                                const Genotype_def geno, const double sm,
                                const double hm, const double u, const double v){
  //genotype frequencies among male
  std::vector<double> sex_freqs(16);
  double sum = 0.0;
  for(int i = 0; i < 16; i++){
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


  for(int i = 0; i < 16; i++){
    sex_freqs.at(i) /= sum;
  }

  std::vector<double> gametes(4);
  std::vector<double> gametes_mut(4);

  for(int i = 0; i < 16; i++){
    //ff
    gametes.at(geno.fy.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //fm
    gametes.at(geno.fy.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //mf
    gametes.at(geno.my.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //mm
    gametes.at(geno.my.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
  }

  gametes_mut.at(0) = (1-v) * gametes.at(0) + u * gametes.at(1);
  gametes_mut.at(1) = (1-u) * gametes.at(1) + v * gametes.at(0);
  gametes_mut.at(2) = (1-v) * gametes.at(2) + u * gametes.at(3);
  gametes_mut.at(3) = (1-u) * gametes.at(3) + v * gametes.at(2);

  return(gametes_mut);
}

std::vector<double> female_gamate(const std::vector<int>& freqs,
                                const Genotype_def geno, const double sf,
                                const double hf, const double u, const double v){
  //genotype frequencies among female
  std::vector<double> sex_freqs(16);
  double sum = 0.0;
  for(int i = 0; i < 16; i++){
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
  for(int i = 0; i < 16; i++){
    sex_freqs.at(i) /= sum;
  }

  std::vector<double> gametes(4);
  std::vector<double> gametes_mut(4);

  for(int i = 0; i < 16; i++){
    //ff
    gametes.at(geno.fy.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //fm
    gametes.at(geno.fy.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //mf
    gametes.at(geno.my.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
    //mm
    gametes.at(geno.my.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * sex_freqs.at(i);
  }

  gametes_mut.at(0) = (1-v) * gametes.at(0) + u * gametes.at(1);
  gametes_mut.at(1) = (1-u) * gametes.at(1) + v * gametes.at(0);
  gametes_mut.at(2) = (1-v) * gametes.at(2) + u * gametes.at(3);
  gametes_mut.at(3) = (1-u) * gametes.at(3) + v * gametes.at(2);

  return(gametes_mut);
}

std::vector<int> next_gen(const std::vector<double>& male_gamate,
                            const std::vector<double>& female_gamate,
                            const int n, std::mt19937& mt){
  std::vector<double> deterministic(16);
  for(int i = 0; i < 16; i++){
    int male = i / 4;
    int female = i % 4;

    deterministic.at(i) = male_gamate.at(male) * female_gamate.at(female);
  }
  std::vector<int> next_generarion(16);
  int n_remain = n;

  for(int i = 0; i < 15; i++){
    double det_remain = 0.0;
    for(int j = i; j < 16; j++){
      det_remain += deterministic.at(j);
    }

    std::binomial_distribution<> choose(n_remain, deterministic.at(i) / det_remain);
    next_generarion.at(i) = choose(mt);
    n_remain -= next_generarion.at(i);
  }
  next_generarion.at(15) = n_remain;
  return(next_generarion);
}

void initialize_freqs_xy(std::vector<int>& freqs, const int n, std::mt19937& mt){
  //XY system
  std::vector<double> m_gamates = {0.25, 0.25, 0.25, 0.25};
  std::vector<double> f_gamates = {0.5, 0.5, 0.0, 0.0};

  freqs = next_gen(m_gamates, f_gamates, n, mt);
}

int main(int argc, char *argv[]){
  //invasion of ZW system

  //parameters
  int n;
  double sf;
  double hf;
  double sm;
  double hm;
  double u;
  double v;
  int num_parallels;

  if(argc == 10){
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%lf", &sf);
    sscanf(argv[3], "%lf", &hf);
    sscanf(argv[4], "%lf", &sm);
    sscanf(argv[5], "%lf", &hm);
    sscanf(argv[6], "%lf", &u);
    sscanf(argv[7], "%lf", &v);
    sscanf(argv[8], "%d", &num_parallels);
  }else{
    std::cerr << "Number of variables not matched!" << std::endl;
    return(1);
  }

  std::cout << n << "\t" << sf << "\t" << hf << "\t" << sm << "\t" << hm << "\t" << u << "\t" << v << std::endl;

  int max_rep = 1000 / num_parallels;
  int max_gen = 5000000 / max_rep / num_parallels;

  std::random_device seed;
  std::mt19937 mt(seed());

  std::vector<int> freqs;
  Genotype_def geno;

  std::string string0 = "distribution";
  std::string string1 = argv[9];
  std::string string2 = ".txt";

  std::ofstream ofs1(string0 + string1 + string2);

  for(int k = 0; k < max_rep; k++){
    initialize_freqs_xy(freqs, n, mt);
    for(int i = -10; i < max_gen; i++){
      for(int j = 0; j < n; j++){
        std::vector<double> m_gamate = male_gamate(freqs, geno, sm, hm, u, v);
        std::vector<double> f_gamate = female_gamate(freqs, geno, sf, hf, u, v);
        freqs = next_gen(m_gamate, f_gamate, n, mt);

        if(i >= 0 && j == n - 1){
          ofs1 << std::setprecision(10) << k*max_gen+i << "\t" << m_gamate.at(0) << "\t" << m_gamate.at(1)
          << "\t" << m_gamate.at(2) << "\t" << m_gamate.at(3)
          << "\t" << f_gamate.at(0) << "\t" << f_gamate.at(1)
          << "\t" << f_gamate.at(2) << "\t" << f_gamate.at(3) << std::endl;
        }
      }
    }
  }

  return(0);
}
