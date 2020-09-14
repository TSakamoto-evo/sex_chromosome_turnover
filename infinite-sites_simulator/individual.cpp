#include "individual.hpp"

Individual::Individual( int input_xy, int input_a, int input_b,
    std::vector<bool> input_neutral0, std::vector<bool> input_neutral1  ){

  locus_xy = input_xy;
  locus_a = input_a;
  locus_b = input_b;
  neutral0 = input_neutral0;
  neutral1 = input_neutral1;
}

void Individual::erase_neutral_site( int position ){
  neutral0.erase(neutral0.begin() + position);
  neutral1.erase(neutral1.begin() + position);
}

int Individual::return_sex(){
  /* For Case 1 */
  if(locus_xy == 0 && locus_a == 0){
    return(1);
  }else{
    return(0);
  }

  /* For Case 2 */
  /*
  if(locus_a != 0){
    return(1);
  }else if(locus_xy == 0){
    return(1);
  }else{
    return(0)
  }
  */
}

void Individual::mutation_b( bool which ){
  int locus_b1 = locus_b / 2;
  int locus_b0 = locus_b % 2;

  if( which ){
    locus_b = 2 * (1 - locus_b1) + locus_b0;
  }else{
    locus_b = 2 * locus_b1 + (1 - locus_b0);
  }
}

void Individual::mutation_a( bool which ){
  int locus_a1 = locus_a / 2;
  int locus_a0 = locus_a % 2;

  if( which ){
    locus_a = 2 * (1 - locus_a1) + locus_a0;
  }else{
    locus_a = 2 * locus_a1 + (1 - locus_a0);
  }
}
