#include<iostream>
#include<fstream>
#include<vector>
#include<random>

class Genotype_def{
public:
  std::vector<bool> my;
  std::vector<bool> mw;
  std::vector<bool> ma;
  std::vector<bool> fy;
  std::vector<bool> fw;
  std::vector<bool> fa;
  std::vector<bool> female;

  /* relate indices to genotypes */
  Genotype_def(){
    for(int i = 0; i < 64; i++){
      my.push_back( (i / 32) % 2 );
      mw.push_back( (i / 16) % 2 );
      ma.push_back( (i /  8) % 2 );
      fy.push_back( (i /  4) % 2 );
      fw.push_back( (i /  2) % 2 );
      fa.push_back( (i /  1) % 2 );

      if(my.at(i) + fy.at(i) > 0 && mw.at(i) + fw.at(i) == 0){
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
    gametes.at(geno.fy.at(i) * 4 + geno.fw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //ffm
    gametes.at(geno.fy.at(i) * 4 + geno.fw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //fmf
    gametes.at(geno.fy.at(i) * 4 + geno.mw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //fmm
    gametes.at(geno.fy.at(i) * 4 + geno.mw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //mff
    gametes.at(geno.my.at(i) * 4 + geno.fw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
    //mfm
    gametes.at(geno.my.at(i) * 4 + geno.fw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //mmf
    gametes.at(geno.my.at(i) * 4 + geno.mw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rm * sex_freqs.at(i);
    //mmm
    gametes.at(geno.my.at(i) * 4 + geno.mw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rm) * sex_freqs.at(i);
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
    gametes.at(geno.fy.at(i) * 4 + geno.fw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //ffm
    gametes.at(geno.fy.at(i) * 4 + geno.fw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //fmf
    gametes.at(geno.fy.at(i) * 4 + geno.mw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //fmm
    gametes.at(geno.fy.at(i) * 4 + geno.mw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //mff
    gametes.at(geno.my.at(i) * 4 + geno.fw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
    //mfm
    gametes.at(geno.my.at(i) * 4 + geno.fw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //mmf
    gametes.at(geno.my.at(i) * 4 + geno.mw.at(i) * 2 + geno.fa.at(i)) += 0.5 * 0.5 * rf * sex_freqs.at(i);
    //mmm
    gametes.at(geno.my.at(i) * 4 + geno.mw.at(i) * 2 + geno.ma.at(i)) += 0.5 * 0.5 * (1.0 - rf) * sex_freqs.at(i);
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

void initialize_freqs_xy(std::vector<int>& freqs, const int n, std::mt19937& mt, const double p){
  //XY system
  std::vector<double> m_gametes = {(1-p)/2.0, p/2.0, 0.0, 0.0, (1-p)/2.0, p/2.0, 0.0, 0.0};
  std::vector<double> f_gametes = {1-p, p, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  freqs = next_gen(m_gametes, f_gametes, n, mt);
}

void introduce_w(std::vector<int>& freqs, const int n, std::mt19937& mt){
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

int introduce_w_beneficial(std::vector<int>& freqs, std::mt19937& mt,
                            const Genotype_def geno){
  std::vector<int> index_ma, index_fa;
  std::vector<double> nums_ma, nums_fa;

  int sum_ma = 0;
  int sum_fa = 0;

  for(int i = 0; i < 64; i++){
    /* pick up haplotypes on which allele B(1) or allele b(0) exist */
    if(geno.ma.at(i) == 1){
      index_ma.push_back(i);
      nums_ma.push_back(freqs.at(i));
      sum_ma += freqs.at(i);
    }
    if(geno.fa.at(i) == 1){
      index_fa.push_back(i);
      nums_fa.push_back(freqs.at(i));
      sum_fa += freqs.at(i);
    }
  }
  if(sum_ma + sum_fa != 0){
    std::bernoulli_distribution select(sum_ma * 1.0 / (sum_ma + sum_fa));
    if(select(mt)){
      std::discrete_distribution<> dist(nums_ma.begin(), nums_ma.end());
      int k = index_ma.at(dist(mt));
      freqs.at(k)--;
      freqs.at(k + 16)++;
      return(1);
    }else{
      std::discrete_distribution<> dist(nums_fa.begin(), nums_fa.end());
      int k = index_fa.at(dist(mt));
      freqs.at(k)--;
      freqs.at(k + 2)++;
      return(1);
    }
  }else{
    return(0);
  }
}

std::vector<int> return_xw(const std::vector<int>& freqs, const Genotype_def geno){
  int num_x = 0;
  int num_w = 0;
  for(int i = 0; i < 64; i++){
    num_x += freqs.at(i) * (2 - geno.my.at(i) - geno.fy.at(i));
    num_w += freqs.at(i) * (geno.mw.at(i) + geno.fw.at(i));
  }
  std::vector<int> ret = {num_x, num_w};
  return(ret);
}

int main(){
  //invasion of ZW system
  //parameters
  int n = 10000;
  double sf = 0.02;
  double hf = 1.0;
  double sm = -0.02;
  double hm = 0.0;
  double u = 0.000001;
  double v = 0.000001;

  double r = 0.005;

  std::vector<double> p_values = {0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999};

  int max_rep = 1000;

  std::random_device seed;
  std::mt19937 mt(seed());

  std::vector<int> freqs;
  Genotype_def geno;

  std::ofstream ofs1("genotype.txt");
  std::ofstream ofs2("fixation_prob.txt", std::ios::app);

  ofs1 << "index\tmy\tfy\tmw\tfw\tma\tfa\tfemale" << std::endl;
  for(int i = 0; i < 64; i++){
    ofs1 << i << "\t" << geno.my.at(i) << "\t" << geno.fy.at(i) << "\t" << geno.mw.at(i)
    << "\t" << geno.fw.at(i) << "\t" << geno.ma.at(i) << "\t" << geno.fa.at(i) << "\t" <<
    geno.female.at(i) << std::endl;
  }

  for(int p_index = 0; p_index < static_cast<int>(p_values.size()); p_index++){
    double p = p_values.at(p_index);
    int reps = 0;
    int regi_w_fix = 0;

    while(1){
      int check = 0;

      while(check == 0){
        initialize_freqs_xy(freqs, n, mt, p);
        check = introduce_w_beneficial(freqs, mt, geno);
      }

      for(int i = 0; 1; i++){
        std::vector<double> m_gamete = male_gamete(freqs, geno, sm, hm, r, u, v);
        std::vector<double> f_gamete = female_gamete(freqs, geno, sf, hf, r, u, v);
        freqs = next_gen(m_gamete, f_gamete, n, mt);

        std::vector<int> xw = return_xw(freqs, geno);

        if(xw.at(0) == 0 || xw.at(1) == 0){
          if(xw.at(0) == 0){
            regi_w_fix++;
          }
          break;
        }
      }
      reps++;
      if((regi_w_fix >= max_rep && reps % 10000 == 0 && reps >= 100000) || reps == 2000000){
        break;
      }
    }
    ofs2 << n << "\t" << sf << "\t" << hf << "\t" << sm << "\t" << hm << "\t" << r << "\t";
    ofs2 << u << "\t" << v << "\t" << p << "\t" << max_rep << "\t1\t";
    ofs2 << regi_w_fix * 1.0 / reps << "\t" << regi_w_fix << "\t" << reps << std::endl;
  }

  return(0);
}
