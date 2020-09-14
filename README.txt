###################################################
Codes for numerical analyses and simulations in
Sakamoto and Innan (2020)
###################################################

Below, usage for each program is provided.

1. numerical_analysis/equation_3.cpp
Equation 3 is calculated.
Set parameter values appropriately in the main function.
For Case 2, sm and hm should be interchanged with sf and hf.

compile: g++ equation_3.cpp -std=c++11 -O3 -o XXX.out
parameters: sm, hm, r, u and v


2. numerical_analysis/equation_7.cpp
Equation 7 is calculated.
Set parameter values appropriately in the main function.
Equilibrium value of p must also be provided such that Mp = 0 in the main function.
For Case 2, sm and hm should be interchanged with sf and hf.

compile: g++ equation_7.cpp -std=c++11 -O3 -o XXX.out
parameters: sm, hm, sf, hf, r, u, v and equ_p

3. numerical_analysis/equation_12.cpp
Equation 12 is calculated.
Set parameter values appropriately in the main function.
Choose Equation 3 or 7 for calculation for phis.
  (Equation 7 should be used except Equalizing selection case.)
Equilibrium value of p must also be provided such that Mp = 0 in the main function
  if Equation 7 is used for calculating phis.
For Case 2, sm and hm should be interchanged with sf and hf.

compile: g++ equation_12.cpp -std=c++11 -O3 -o XXX.out
parameters: sm, hm, sf, hf, u, v, pop_size and equ_p

4. numerical_analysis/equation_e4_negative.cpp
E(X^i_j) in Equation E4 is calculated.
It should be used for negative selection case.
Set parameter values appropriately in the main function.
For Case 2, sm and hm should be interchanged with sf and hf.

compile: g++ equation_e4_negative.cpp -std=c++11 -O3 -o XXX.out
parameters: sm, hm, sf, hf and r

5. numerical_analysis/equation_e4_positive.cpp
E(X^i_j) in Equation E4 is calculated.
It should be used for positive selection case.
Set parameter values appropriately in the main function.
For Case 2, sm and hm should be interchanged with sf and hf.

compile: g++ equation_e4_positive.cpp -std=c++11 -O3 -o XXX.out
parameters: sm, hm, sf, hf and r

6. conditional_simulator/simulator_for_conditional_case1.cpp
Establishment probability is obtained by stochastic simulation of Wright-Fisher model.
Set parameter values appropriately in the main function.

compile: g++ simulator_for_conditional_case1.cpp -std=c++11 -O3 -o XXX.out
parameters: n, sm, hm, sf, hf, u, v and r

7. conditional_simulator/simulator_for_conditional_case1.cpp
Establishment probability is obtained by stochastic simulation of Wright-Fisher model.
Set parameter values appropriately in the main function.

compile: g++ simulator_for_conditional_case1.cpp -std=c++11 -O3 -o XXX.out
parameters: n, sm, hm, sf, hf, u, v and r

8. unconditional_simulator
Establishment probability is obtained by stochastic simulation of Wright-Fisher model.
Set parameter values appropriately by command line arguments.
First, run generate_stationary.cpp to obtain stationary distribution of frequency of allele B.
Then, unconditional probability is obtained by simulator_for_unconditional_caseX.cpp.

generate_stationary.cpp
  compile: g++ generate_stationary.cpp -std=c++11 -O3 -o XXX.out
  execution: XXX.out n sf hf sm hm u v 1 1

simulator_for_unconditional_caseX.cpp
  compile: g++ simulator_for_unconditional_caseX.cpp -std=c++11 -O3 -o XXX.out
  execution: XXX.out n sf hf sm hm r u v

9. infinite-sites simulator
Nucleotide diversity within and between populations are simulated.
Set parameter values appropriately in main function.
In default, Case 1 is simulated, but Case 2 can also be simulated by modifying return_sex function in individual.cpp appropriately.
Parameters mut_rate and recomb_rate represents rates for entire region, while mut_rate_at_b represents rate at locus B/b.

compile: *.cpp -std=c++11 -O3 -o XXX.out
parameters: pop_size, recomb_rate, mut_rate, mut_rate_at_b, sm, hm, sf, hf
