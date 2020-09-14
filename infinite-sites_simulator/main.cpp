#include "parameters.hpp"
#include "population.hpp"
#include <fstream>
#include <iostream>

int main(){
  Parameters paras;

  paras.pop_size = 1000;
  paras.recomb_rate = 0.002;
  paras.mut_rate = 0.01;
  paras.mut_rate_at_b = 0.000001;

  paras.sm = 0.025;
  paras.hm = 0.5;
  paras.sf = -0.02;
  paras.hf = 0.5;

  int seps = 100;
  int samples = 20;

  std::ofstream ofs1("pi_dist.txt");
  std::ofstream ofs2("pi_w0.txt");
  std::ofstream ofs3("pi_w1.txt");
  std::ofstream ofs4("pi_b.txt");
  std::ofstream ofs5("freqs.txt");

  Population pop(paras, 0.5);

  // burn-in 200n generations
  for(int i = 0; i < 200 * paras.pop_size; i++){
    pop.one_generation();

    if( i % paras.pop_size == 0 ){
      std::cout << i << "\t" << pop.return_mut_sites() << std::endl;
    }
  }

  int polymorphic_b = 0;

  // continue burn-in period if locus B is monomorphic
  while(polymorphic_b == 0){
    std::vector<int> conf_ab = pop.population_conf_b();
    int allele0 = conf_ab.at(0) + conf_ab.at(2) + conf_ab.at(4) + conf_ab.at(6);
    int allele1 = conf_ab.at(1) + conf_ab.at(3) + conf_ab.at(5) + conf_ab.at(7);

    if(allele0 > 4 && allele1 > 4){
      polymorphic_b = 1;
    }else{
      for(int i = 0; i < 2 * paras.pop_size; i++){
        pop.one_generation();

        if( i % paras.pop_size == 0 ){
          std::cout << i << "\t" << pop.return_mut_sites() << std::endl;
        }
      }
    }
  }


  Population pop_regi = pop;

  // register nucleotide diversities just before allele A arises
  {
    std::vector<double> pi_dist = pop.calculate_pi(seps, samples);

    ofs1 << -1;
    for(int j = 0; j < seps; j++){
      ofs1 << "\t" << pi_dist.at(j);
    }
    ofs1 << std::endl;

    std::vector<int> conf_ab = pop.population_conf_b();
    ofs5 << -1;
    for(int i: conf_ab){
      ofs5 << "\t" << i;
    }
    ofs5 << std::endl;

    std::vector<std::vector<double>> pi_dists = pop.calculate_pi_conditional_a(seps, samples);

    ofs2 << -1; ofs3 << -1; ofs4 << -1;
    for(int j = 0; j < seps; j++){
      ofs2 << "\t" << pi_dists.at(0).at(j);
      ofs3 << "\t" << pi_dists.at(1).at(j);
      ofs4 << "\t" << pi_dists.at(2).at(j);
    }
    ofs2 << std::endl; ofs3 << std::endl; ofs4 << std::endl;
  }

  // repeat introduction of allele A until allele A establishes
  while(1){
    int est = 1;
    int gene = 0;
    pop.introduce_a();

    while(est == 1){
      pop.one_generation();
      gene++;
      std::vector<int> conf = pop.population_conf();

      if(conf.at(5) + conf.at(7) == 0){
        std::vector<int> conf_ab = pop.population_conf_b();
        int allele0 = conf_ab.at(0) + conf_ab.at(2) + conf_ab.at(4) + conf_ab.at(6);
        int allele1 = conf_ab.at(1) + conf_ab.at(3) + conf_ab.at(5) + conf_ab.at(7);

        if(allele0 > 0 && allele1 > 0){
          est = 2;
        }else{
          est = 3;
        }
      }
      if(conf.at(6) + conf.at(7) == 0){
        est = 0;
      }
    }

    std::cout << est;

    if(est == 2){
      std::cout << "\n" << gene << std::endl;
      break;
    }else{
      pop = pop_regi;
      pop.renew_random();
    }
  }

  // run simulation for 50N generations
  for(int i = 0; i < 50 * paras.pop_size; i++){
    pop.one_generation();
    if( i % paras.pop_size == 0 ){
      std::cout << i << "\t" << pop.return_mut_sites() << std::endl;
    }

    if( i % (paras.pop_size / 10) == 0){
      std::vector<double> pi_dist = pop.calculate_pi(seps, samples);

      ofs1 << i;
      for(int j = 0; j < seps; j++){
        ofs1 << "\t" << pi_dist.at(j);
      }
      ofs1 << std::endl;


      std::vector<std::vector<double>> pi_dists = pop.calculate_pi_conditional_a(seps, samples);

      ofs2 << i; ofs3 << i; ofs4 << i;
      for(int j = 0; j < seps; j++){
        ofs2 << "\t" << pi_dists.at(0).at(j);
        ofs3 << "\t" << pi_dists.at(1).at(j);
        ofs4 << "\t" << pi_dists.at(2).at(j);
      }
      ofs2 << std::endl; ofs3 << std::endl; ofs4 << std::endl;

      std::vector<int> conf_ab = pop.population_conf_b();
      ofs5 << i;
      for(int i: conf_ab){
        ofs5 << "\t" << i;
      }
      ofs5 << std::endl;

    }
  }
  return(0);
}
