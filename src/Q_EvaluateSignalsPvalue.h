/**
 * @file
 * @authors  Peter Hansen <peter.hansen@charite.de>
 * @version 0.01
 *
 * \brief Is used by Q_EvaluateSignals.h and contains everything
 * which has to do with the calculation of p-values.
 * 
 */
 
#ifndef Q_EVALUATE_SIGNALS_PVALUE_H
#define Q_EVALUATE_SIGNALS_PVALUE_H

#include "Q_ReadInFiles.h"
#include <boost/math/distributions/binomial.hpp>
#include <cmath>


long double get_probability_of_success(int t_f, int t_r, int l, int q_min, int q_max);
int get_saturation_score(Chromosome &chromosome, int radius);
int get_saturation_pvalues(Chromosome &chromosome, int q_min, int q_max, int length, int t_f, int t_r, int c_f, int c_r, bool ctrl);

struct Pvalues
{
	// indicates whether a control is given or not
	bool has_control;

	// length of the considered regions (q_max)
	int n;
	
	// probability of success for treatment
	long double p_t;
	
	// probability of success for control
	long double p_c;
	
	// map containing p-values für differences -n,...,0,...n
	std::map<int, long double> p_values;
	std::map<int, long double> p_values_10log;
	
	//long double get_p_value(int d);
		
	// constructor only treatment
	Pvalues(int length, int q_min, int q_max, int t_f, int t_r);
	
	// constructor treatment and control
	Pvalues(int length, int q_min, int q_max, int t_f, int t_r, int c_f, int c_r);
	
	// function that returns the P-value for a given distance
	// and throws errors if a p-value is not defined for a given distance
	// check if d is in the allowed range -n,...,0,...n
	long double get_p_value(int k, bool lg10);
	
};


#endif
