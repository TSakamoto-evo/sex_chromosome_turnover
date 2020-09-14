#ifndef INDIVIDUAL
#define INDIVIDUAL

#include "parameters.hpp"
#include <vector>

class Individual{
private:
  int locus_xy;
  int locus_a;
  int locus_b;
  std::vector<bool> neutral0;
  std::vector<bool> neutral1;

public:
  Individual( int input_xy, int input_a, int input_b,
    std::vector<bool> input_neutral0, std::vector<bool> input_neutral1 );
  int return_xy(){ return locus_xy; };
  int return_a(){ return locus_a; };
  int return_b(){ return locus_b; };
  std::vector<bool> return_neutral0(){ return neutral0; };
  std::vector<bool> return_neutral1(){ return neutral1; };
  int return_pos( int pos ){ return (neutral0.at(pos) + neutral1.at(pos)); };
  int return_sex();

  void add_neutral_site( bool add0, bool add1 ){ neutral0.push_back(add0); neutral1.push_back(add1); };
  void erase_neutral_site( int position );
  void mutation_b( bool which );
  void mutation_a( bool which );
};

#endif
