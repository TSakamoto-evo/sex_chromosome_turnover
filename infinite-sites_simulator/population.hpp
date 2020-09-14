#ifndef POPULATION
#define POPULATION

#include "parameters.hpp"
#include "individual.hpp"
#include <random>

class Population{
private:
  Parameters para;
  std::vector<Individual> inds;
  std::vector<double> index_pos;
  std::mt19937 mt;

  std::vector<double> female_fit;
  std::vector<double> male_fit;

public:
  Population( Parameters paras, double b_freq );
  void set_selection_vectors();
  void mating();
  void mutation();
  void erase_unused_positions();
  void one_generation();
  std::vector<bool> make_haplotype_recomb( int parent, int& locus_xy, int& locus_a, int& locus_b , int reco_num);
  int return_mut_sites(){ return para.num_sites; };
  std::vector<double> calculate_pi(int seps, int num_samples);
  std::vector<std::vector<double>> calculate_pi_conditional_a(int seps, int num_samples);
  void introduce_a();
  std::vector<int> population_conf();
  std::vector<int> population_conf_b();
  void renew_random();
};

#endif
