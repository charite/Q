#include "Q_EvaluateSignalsPvalue.h"
#include "Q_ReadInFiles.h"

int get_saturation_score(Chromosome &chromosome, int radius)
{
	for(int i=0;i<chromosome.sum_num;i++)
	{
		chromosome.SUMMITS[i].saturation_score = (chromosome.SUMMITS[i].q_beg_chip+chromosome.SUMMITS[i].q_end_chip)-(chromosome.SUMMITS[i].q_beg_ctrl+chromosome.SUMMITS[i].q_end_ctrl);
		if(chromosome.SUMMITS[i].saturation_score<0){chromosome.SUMMITS[i].saturation_score=0;}
		chromosome.SUMMITS[i].saturation_score = (double) (chromosome.SUMMITS[i].saturation_score/(2*radius));
	}
	return 0;
}

/////////////////////////////////////////////////////////////////


// constructor for the case without control
Pvalues::Pvalues(int length, int q_min, int q_max, int t_f, int t_r)
{
	has_control=false;
	n=2*q_max;
		
	// calculate probility of success for treatment
	// --------------------------------------------
	long double lambda_t = (q_max-q_min)*(t_r/(long double)length);
	long double rate_t   = (t_f/(long double)length)*(1-expl(-lambda_t));
	p_t = 1-exp(-rate_t*2);

	// init binomial distribution
	// --------------------------
	boost::math::binomial_distribution<long double> BinomialTreatment(n,p_t);

	// init p-values from n to 0
	// -------------------------
	long double pval=0;
	for(int k=n;0<=k;k--)
	{
		pval=pval+pdf(BinomialTreatment,(k));
		p_values[k]=pval;
		p_values_10log[k]=-log10l(pval);
	}
}

// constructor for the case treatment and control
Pvalues::Pvalues(int length, int q_min, int q_max, int t_f, int t_r, int c_f, int c_r)
{
	n=2*q_max;
	has_control=true;
	
	// calculate probility of success for treatment and control
	// --------------------------------------------------------
	
	long double lambda_t = (q_max-q_min) * (t_r/(long double)length);
	long double rate_t = (t_f/(long double)length) * (1-expl(-lambda_t));
	p_t = 1 - exp(-rate_t*2);
	
	long double lambda_c = (q_max-q_min) * (c_r/(long double)length);
	long double rate_c = (c_f/(long double)length) * (1-expl(-lambda_c));
	p_c = 1 - exp(-rate_c*2);
	
	// init binomial distribution
	// --------------------------
	
	boost::math::binomial_distribution<long double> BinomialTreatment(n,p_t);
	boost::math::binomial_distribution<long double> BinomialControl(n,p_c);	
	
	// init p-values from -n to n
	// --------------------------
	
	std::map<int, long double> DiffBinom;
	for(int d=-n;d<=n;d++)
	{
		for(int k=0;k<=(n-abs(d));k++)
		{
			if(0<=d)
			{
				DiffBinom[d]=DiffBinom[d]+pdf(BinomialTreatment,(k+abs(d)))*pdf(BinomialControl,k);
			}
			else
			{
				DiffBinom[d]=DiffBinom[d]+pdf(BinomialControl,k+abs(d))*pdf(BinomialTreatment,k);
			}
		}
	}
	
	long double pval=0;
	for(int d=n;d>=-n;d--)
	{
		pval=pval+DiffBinom[d];
		p_values[d]=pval;
		p_values_10log[d]=-log10l(pval);
	}
}

long double Pvalues::get_p_value(int k, bool lg10)
{
	if(has_control==false)
	{
		if( k<0 || n<k )
		{
			std::cerr << "Error in Pvalues::get_p_value(int k): ";
			std::cerr << "k="<< k << " is outside the allowed range (0 to " << n << ")\n";
			return(-1);
		}
		else
		{
			if(lg10){return(p_values_10log[k]);}
			else{return(p_values[k]);}
		}
	}
	else
	{
		if( k<-n || n<k )
		{
			std::cerr << "Error in Pvalues::get_p_value(int k): ";
			std::cerr << "k="<< k << " is outside the allowed range (" << -n << " to " << n << ")\n";
			return(-1);
		}
		else
		{
			if(lg10){return(p_values_10log[k]);}
			else{return(p_values[k]);}
		}
	}
}


int get_saturation_pvalues(Chromosome &chromosome, int q_min, int q_max, int length, int t_f, int t_r, int c_f, int c_r, bool ctrl)
{
	if(ctrl)
	{
		Pvalues pvalues(length, q_min, q_max, t_f, t_r, c_f, c_r);
		for(int i=0;i<chromosome.sum_num;i++)
		{
			chromosome.SUMMITS[i].p_value=pvalues.get_p_value(
			chromosome.SUMMITS[i].q_beg_chip+chromosome.SUMMITS[i].q_end_chip
			-(chromosome.SUMMITS[i].q_beg_ctrl+chromosome.SUMMITS[i].q_end_ctrl),false);
		}
	}
	else
	{
		Pvalues pvalues(length, q_min, q_max, t_f,t_r);
		for(int i=0;i<chromosome.sum_num;i++)
		{
			chromosome.SUMMITS[i].p_value=pvalues.get_p_value(chromosome.SUMMITS[i].q_beg_chip+chromosome.SUMMITS[i].q_end_chip, false);
		}
	}
	return 0;
}

long double get_probability_of_success(int t_f, int t_r, int l, int q_min, int q_max)
{
	long double lambda_t=(q_max-q_min) * (long double)t_r/l;
//	std::cout << "lambda_t: " << lambda_t << "\n";
	double r_t = 2 * (long double)t_f/l * (1-expl(-lambda_t));
	long double p_t = 1-expl(-r_t);
	return p_t;
}
