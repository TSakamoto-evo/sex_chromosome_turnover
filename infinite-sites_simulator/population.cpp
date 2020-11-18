#include "population.hpp"

Population::Population( Parameters paras, double b_freq ){
  para = paras;
  para.num_sites = 0;

  std::random_device seed;
  std::mt19937 mt_tmp(seed());
  mt = mt_tmp;

  std::vector<double> index_pos_tmp;
  index_pos = index_pos_tmp;

  std::vector<double> fit_tmp;
  female_fit = fit_tmp;
  male_fit = fit_tmp;

  std::bernoulli_distribution pb(b_freq);
  std::bernoulli_distribution half(0.5);
  std::vector<bool> tmp_neutral(0);

  std::vector<Individual> inds_tmp;
  for(int i = 0; i < para.pop_size; i++){
    int mother_xy = 0;
    int father_xy = half(mt);
    int mother_a = 0;
    int father_a = 0;
    int mother_b = pb(mt);
    int father_b = pb(mt);

    inds_tmp.emplace_back(mother_xy*2+father_xy, mother_a*2+father_a, mother_b*2+father_b, tmp_neutral, tmp_neutral);
  }
  inds = inds_tmp;
}

void Population::set_selection_vectors(){
  /* selection */
  std::vector<double> new_female_fit;
  std::vector<double> new_male_fit;

  for(int i = 0; i < para.pop_size; i++){
    int locus_b = inds.at(i).return_b();
    if(inds.at(i).return_sex() == 1){
      new_male_fit.push_back(0.0);
      if(locus_b == 0){
        new_female_fit.push_back(1.0);
      }else if(locus_b == 3){
        new_female_fit.push_back(1.0 + para.sf);
      }else{
        new_female_fit.push_back(1.0 + para.sf * para.hf);
      }
    }else{
      new_female_fit.push_back(0.0);
      if(locus_b == 0){
        new_male_fit.push_back(1.0);
      }else if(locus_b == 3){
        new_male_fit.push_back(1.0 + para.sm);
      }else{
        new_male_fit.push_back(1.0 + para.sm * para.hm);
      }
    }
  }

  female_fit = new_female_fit;
  male_fit = new_male_fit;
}

std::vector<bool> Population::make_haplotype_recomb( int parent, int& locus_xy, int& locus_a, int& locus_b, int reco_num ){
  std::bernoulli_distribution half(0.5);
  int left_edge = half(mt);

  if(half(mt)){
    locus_xy = inds.at(parent).return_xy() / 2;
  }else{
    locus_xy = inds.at(parent).return_xy() % 2;
  }

  if(reco_num == 0){
    if(left_edge == 1){
      locus_a = inds.at(parent).return_a() / 2;
      locus_b = inds.at(parent).return_b() / 2;
      return( inds.at(parent).return_neutral1() );
    }else{
      locus_a = inds.at(parent).return_a() % 2;
      locus_b = inds.at(parent).return_b() % 2;
      return( inds.at(parent).return_neutral0() );
    }
  }else{
    std::vector<double> r_pos;
    std::uniform_real_distribution<> recomb_pos(0.0, 1.0);

    for(int i = 0; i < reco_num; i++){
      r_pos.push_back(recomb_pos(mt));
    }

    std::vector<bool> neutral_father = inds.at(parent).return_neutral0();
    std::vector<bool> neutral_mother = inds.at(parent).return_neutral1();

    std::vector<bool> neutral_tmp(para.num_sites);

    if(left_edge == 1){
      {
        int count = 0;
        for(double x: r_pos){
          if(x < 0.25){
            count++;
          }
        }

        if(count % 2 == 0){
          locus_a = inds.at(parent).return_a() / 2;
        }else{
          locus_a = inds.at(parent).return_a() % 2;
        }
      }

      {
        int count = 0;
        for(double x: r_pos){
          if(x < 0.75){
            count++;
          }
        }

        if(count % 2 == 0){
          locus_b = inds.at(parent).return_b() / 2;
        }else{
          locus_b = inds.at(parent).return_b() % 2;
        }
      }

      for(int i = 0; i < para.num_sites; i++){
        int count = 0;
        for(double x: r_pos){
          if(x < index_pos.at(i)){
            count++;
          }
        }

        if(count % 2 == 0){
          neutral_tmp.at(i) = neutral_mother.at(i);
        }else{
          neutral_tmp.at(i) = neutral_father.at(i);
        }
      }
    }else{
      {
        int count = 0;
        for(double x: r_pos){
          if(x < 0.25){
            count++;
          }
        }

        if(count % 2 == 0){
          locus_a = inds.at(parent).return_a() % 2;
        }else{
          locus_a = inds.at(parent).return_a() / 2;
        }
      }

      {
        int count = 0;
        for(double x: r_pos){
          if(x < 0.75){
            count++;
          }
        }

        if(count % 2 == 0){
          locus_b = inds.at(parent).return_b() % 2;
        }else{
          locus_b = inds.at(parent).return_b() / 2;
        }
      }

      for(int i = 0; i < para.num_sites; i++){
        int count = 0;
        for(double x: r_pos){
          if(x < index_pos.at(i)){
            count++;
          }
        }

        if(count % 2 == 0){
          neutral_tmp.at(i) = neutral_father.at(i);
        }else{
          neutral_tmp.at(i) = neutral_mother.at(i);
        }
      }
    }

    return(neutral_tmp);
  }
}

void Population::mating(){
  std::discrete_distribution<> choose_female(female_fit.begin(), female_fit.end());
  std::discrete_distribution<> choose_male(male_fit.begin(), male_fit.end());

  std::poisson_distribution<> reco_gen(para.recomb_rate * 2 * para.pop_size);
  int num_reco = reco_gen(mt);

  std::vector<int> mother_reco(para.pop_size, 0);
  std::vector<int> father_reco(para.pop_size, 0);

  std::uniform_int_distribution<> choose_ind(0, para.pop_size - 1);
  std::bernoulli_distribution half(0.5);

  for(int i = 0; i < num_reco; i++){
    int ind = choose_ind(mt);
    if(half(mt)){
      mother_reco.at(ind)++;
    }else{
      father_reco.at(ind)++;
    }
  }

  std::vector<Individual> next_inds;
  for(int i = 0; i < para.pop_size; i++){
    int mother = choose_female(mt);
    int m_xy, m_a, m_b;
    std::vector<bool> m_neutral = make_haplotype_recomb(mother, m_xy, m_a, m_b, mother_reco.at(i));

    int father = choose_male(mt);
    int f_xy, f_a, f_b;
    std::vector<bool> f_neutral = make_haplotype_recomb(father, f_xy, f_a, f_b, father_reco.at(i));

    next_inds.emplace_back(2*m_xy+f_xy, 2*m_a+f_a, 2*m_b+f_b, f_neutral, m_neutral);
  }

  inds = next_inds;
}

void Population::mutation(){
  std::poisson_distribution<> mut_gen(para.mut_rate * 2 * para.pop_size);
  std::poisson_distribution<> mut_gen_b(para.mut_rate_at_b * 2 * para.pop_size);

  std::uniform_int_distribution<> choose_ind(0, para.pop_size - 1);
  std::uniform_real_distribution<> choose_pos(0.0, 1.0);
  std::bernoulli_distribution half(0.5);

  int mut_num = mut_gen(mt);
  int mut_num_b = mut_gen_b(mt);

  for(int i = 0; i < mut_num; i++){
    int ind = choose_ind(mt);
    if(half(mt)){
      inds.at(ind).add_neutral_site(1, 0);
    }else{
      inds.at(ind).add_neutral_site(0, 1);
    }

    for(int j = 0; j < para.pop_size; j++){
      if(j != ind){
        inds.at(j).add_neutral_site(0, 0);
      }
    }
    index_pos.push_back(choose_pos(mt));
    para.num_sites++;
  }

  for(int i = 0; i < mut_num_b; i++){
    int ind = choose_ind(mt);
    if(half(mt)){
      inds.at(ind).mutation_b(1);
    }else{
      inds.at(ind).mutation_b(0);
    }
  }
}

void Population::erase_unused_positions(){
  std::vector<int> not_used;
  for(int i = para.num_sites - 1; i >= 0; i--){
    int exist = 0;
    int ref = inds.at(0).return_pos(i);

    for(int j = 1; j < para.pop_size; j++){
      if(inds.at(j).return_pos(i) != ref){
        exist = 1;
        break;
      }
    }

    if(exist == 0){
      not_used.push_back(i);
    }
  }

  for(int i = 0; i < static_cast<int>(not_used.size()); i++){
    for(int j = 0; j < para.pop_size; j++){
      inds.at(j).erase_neutral_site(not_used.at(i));
    }
    index_pos.erase(index_pos.begin() + not_used.at(i));
    para.num_sites--;
  }
}

void Population::one_generation(){
  set_selection_vectors();
  mating();
  erase_unused_positions();
  mutation();
}

std::vector<double> Population::calculate_pi(int seps, int num_samples){
  std::vector<std::vector<bool>> samples;
  std::bernoulli_distribution half(0.5);

  std::vector<int> ind_indices;
  for(int i = 0; i < para.pop_size - 1; i++){
    ind_indices.push_back(i);
  }
  std::shuffle( ind_indices.begin(), ind_indices.end(), mt );

  std::vector<int> diffs(seps);

  for(int i = 0; i < num_samples; i++){
    int ind = ind_indices.at(i);
    if(half(mt)){
      samples.push_back(inds.at(ind).return_neutral0());
    }else{
      samples.push_back(inds.at(ind).return_neutral1());
    }
  }

  for(int i = 0; i < para.num_sites; i++){
    int derived = 0;
    for(int j = 0; j < num_samples; j++){
      if(samples.at(j).at(i) == 1){
        derived++;
      }
    }
    int window = std::floor(index_pos.at(i) * seps);
    diffs.at(window) += derived * (num_samples - derived);
  }

  std::vector<double> ret;
  for(int i = 0; i < seps; i++){
    ret.push_back( 2.0 * diffs.at(i) / num_samples / (num_samples - 1) );
  }
  return(ret);
}

std::vector<std::vector<double>> Population::calculate_pi_conditional_a(int seps, int num_samples){
  std::vector<std::vector<bool>> samples0;
  std::vector<std::vector<bool>> samples1;

  std::vector<int> type0;
  std::vector<int> type1;

  for(int i = 0; i < para.pop_size; i++){
    int locus_a = inds.at(i).return_a();
    int hap1 = locus_a / 2;
    int hap0 = locus_a % 2;

    if(hap0 == 1){
      type1.push_back(i);
    }else{
      type0.push_back(i);
    }

    if(hap1 == 1){
      type1.push_back(i + para.pop_size);
    }else{
      type0.push_back(i + para.pop_size);
    }
  }

  std::shuffle( type0.begin(), type0.end(), mt);
  std::shuffle( type1.begin(), type1.end(), mt);

  if( static_cast<int>(type0.size()) > num_samples && static_cast<int>(type1.size()) > num_samples ){
    for(int i = 0; i < num_samples; i++){
      int index0 = type0.at(i);
      int index1 = type1.at(i);

      if(index0 >= para.pop_size){
        int ind = index0 - para.pop_size;
        samples0.push_back(inds.at(ind).return_neutral1());
      }else{
        int ind = index0;
        samples0.push_back(inds.at(ind).return_neutral0());
      }

      if(index1 >= para.pop_size){
        int ind = index1 - para.pop_size;
        samples1.push_back(inds.at(ind).return_neutral1());
      }else{
        int ind = index1;
        samples1.push_back(inds.at(ind).return_neutral0());
      }
    }

    std::vector<int> diffs_w0(seps);
    std::vector<int> diffs_w1(seps);
    std::vector<int> diffs_b(seps);

    for(int i = 0; i < para.num_sites; i++){
      int derived0 = 0;
      int derived1 = 0;

      for(int j = 0; j < num_samples; j++){
        if(samples0.at(j).at(i) == 1){
          derived0++;
        }
        if(samples1.at(j).at(i) == 1){
          derived1++;
        }
      }
      int window = std::floor(index_pos.at(i) * seps);
      diffs_w0.at(window) += derived0 * (num_samples - derived0);
      diffs_w1.at(window) += derived1 * (num_samples - derived1);
      diffs_b.at(window) += derived0 * (num_samples - derived1) + derived1 * (num_samples - derived0);
    }

    std::vector<double> ret_w0;
    std::vector<double> ret_w1;
    std::vector<double> ret_b;

    for(int i = 0; i < seps; i++){
      ret_w0.push_back( 2.0 * diffs_w0.at(i) / num_samples / (num_samples - 1) );
      ret_w1.push_back( 2.0 * diffs_w1.at(i) / num_samples / (num_samples - 1) );
      ret_b.push_back( 1.0 * diffs_b.at(i) / num_samples / num_samples );
    }

    std::vector<std::vector<double>> ret = {ret_w0, ret_w1, ret_b};
    return(ret);
  }else{
    std::vector<double> ret_w0(seps, 0.0);
    std::vector<double> ret_w1(seps, 0.0);
    std::vector<double> ret_b(seps, 0.0);

    std::vector<std::vector<double>> ret = {ret_w0, ret_w1, ret_b};
    return(ret);
  }
}

void Population::introduce_a(){
  std::uniform_int_distribution<> choose_ind(0, para.pop_size - 1);
  std::bernoulli_distribution half(0.5);

  int ind = choose_ind(mt);
  inds.at(ind).mutation_a( half(mt) );
}

std::vector<int> Population::population_conf(){
  std::vector<int> ret(8);

  for(int i = 0; i < para.pop_size; i++){
    if(inds.at(i).return_sex() == 1){
      int a1 = inds.at(i).return_a() / 2;
      int a0 = inds.at(i).return_a() % 2;
      int xy1 = inds.at(i).return_xy() / 2;
      int xy0 = inds.at(i).return_xy() % 2;

      ret.at(2*a1+xy1)++;
      ret.at(2*a0+xy0)++;
    }else{
      int a1 = inds.at(i).return_a() / 2;
      int a0 = inds.at(i).return_a() % 2;
      int xy1 = inds.at(i).return_xy() / 2;
      int xy0 = inds.at(i).return_xy() % 2;

      ret.at(4+2*a1+xy1)++;
      ret.at(4+2*a0+xy0)++;
    }
  }
  return(ret);
}

std::vector<int> Population::population_conf_b(){
  std::vector<int> ret(8);

  for(int i = 0; i < para.pop_size; i++){
    if(inds.at(i).return_sex() == 1){
      int a1 = inds.at(i).return_a() / 2;
      int a0 = inds.at(i).return_a() % 2;
      int b1 = inds.at(i).return_b() / 2;
      int b0 = inds.at(i).return_b() % 2;

      ret.at(2*a1+b1)++;
      ret.at(2*a0+b0)++;
    }else{
      int a1 = inds.at(i).return_a() / 2;
      int a0 = inds.at(i).return_a() % 2;
      int b1 = inds.at(i).return_b() / 2;
      int b0 = inds.at(i).return_b() % 2;

      ret.at(4+2*a1+b1)++;
      ret.at(4+2*a0+b0)++;
    }
  }
  return(ret);
}

void Population::renew_random(){
  std::random_device seed;
  std::mt19937 mt_tmp(seed());
  mt = mt_tmp;
}
