#include "Q_ReadInFiles.h"
#include "Q_ParseCommandLine.h"
#include "Q_WriteResults.h"


int writeChromosomeInfo(seqan::CharString chromosome_info_file, std::vector<Chromosome> &chromosome, int chr_num, int q_min, int q_max, bool keep_dup)
{
	std::ofstream OUT;
	OUT.open (toCString(chromosome_info_file));	

	OUT	<< "chromosome"                       << "\t"
		<< "length"                           << "\t"
		<< "# summits"                        << "\t"		
		<< "# f-hits-treatment"               << "\t"
		<< "# r-hits-treatment"               << "\t"
		<< "# hits-treatment"                 << "\t"
		<< "# redundant-hits-removed-treatment"    << "\t"
		<< "duplication-rate-treatment"            << "\t"
		<< "# f-hits-control"                 << "\t"
		<< "# r-hits-control"                 << "\t"	
		<< "# hits-control"                   << "\t"
		<< "# redundant-hits-removed-control" << "\t"			
		<< "duplication-rate-control"         << "\t"
		<< "# qfrags-treatment-expected"      << "\t"
		<< "# qfrags-treatment-observed"      << "\t"		
		<< "q-enrichment-score"               << "\n";


	long int len_sum                          = 0;
	int sum_num_sum                           = 0;
	int f_hit_chip_num_sum                    = 0;
	int r_hit_chip_num_sum                    = 0;
	int hit_chip_num_sum                      = 0;
	int redundant_hits_removed_num_chip_sum   = 0;	
	int redundant_hits_num_chip_sum           = 0;	
	int f_hit_ctrl_num_sum                    = 0;
	int r_hit_ctrl_num_sum                    = 0;	
	int hit_ctrl_num_sum                      = 0;	
	int redundant_hits_removed_num_ctrl_sum   = 0;
	int redundant_hits_num_ctrl_sum           = 0;		
	int qfrag_num_sum                         = 0;

	// sum up infos for all chromosomes
	for(int i=0;i<chr_num;i++)
	{
		// calculate duplication rates
		float dup_rate_chip=0;
		if(0<chromosome[i].hit_num_chip)
		{
			dup_rate_chip=(float)chromosome[i].redundant_hits_num_chip/(chromosome[i].hit_num_chip+chromosome[i].redundant_hits_removed_num_chip);
		}

		float dup_rate_ctrl=0;
		if(0<chromosome[i].hit_num_ctrl)
		{
			dup_rate_ctrl=(float)chromosome[i].redundant_hits_num_ctrl/(chromosome[i].hit_num_ctrl+chromosome[i].redundant_hits_removed_num_ctrl);
		}
		
		// calculate expected number of qfrags
		double qfrags_num_exp=0;
		if(0<chromosome[i].f_hit_num_chip && 0<chromosome[i].r_hit_num_chip)
		{
			qfrags_num_exp=chromosome[i].f_hit_num_chip * (q_max-q_min) * (double)chromosome[i].r_hit_num_chip/chromosome[i].len;
		}
		double q_enrichment_score=0;
		if(0<qfrags_num_exp)
		{
			q_enrichment_score=chromosome[i].qfrag_num_chip/qfrags_num_exp;
		}
		if(keep_dup)
		{
			q_enrichment_score=-1;
		}
		
		OUT.precision(2);
		OUT.setf(std::ios_base::fixed);	
		OUT	<< chromosome[i].name                                  << "\t"
			<< chromosome[i].len                                   << "\t"
			<< chromosome[i].sum_num                               << "\t"			
			<< chromosome[i].f_hit_num_chip                        << "\t"
			<< chromosome[i].r_hit_num_chip                        << "\t"
			<< chromosome[i].hit_num_chip                          << "\t"
			<< chromosome[i].redundant_hits_removed_num_chip       << "\t"
			<< dup_rate_chip                                       << "\t"
			<< chromosome[i].f_hit_num_ctrl                        << "\t"	
			<< chromosome[i].r_hit_num_ctrl                        << "\t"			
			<< chromosome[i].hit_num_ctrl                          << "\t"			
			<< chromosome[i].redundant_hits_removed_num_ctrl       << "\t"
			<< dup_rate_ctrl                                       << "\t"
			<< std::fixed << (int)round(qfrags_num_exp)            << "\t"				
			<< chromosome[i].qfrag_num_chip                        << "\t"
			<< q_enrichment_score                                  << "\n";
			
		len_sum            = len_sum            + chromosome[i].len;
		sum_num_sum        = sum_num_sum        + chromosome[i].sum_num;		
		f_hit_chip_num_sum = f_hit_chip_num_sum + chromosome[i].f_hit_num_chip;
		r_hit_chip_num_sum = r_hit_chip_num_sum + chromosome[i].r_hit_num_chip;
		hit_chip_num_sum   = hit_chip_num_sum   + chromosome[i].hit_num_chip;
		redundant_hits_removed_num_chip_sum   = redundant_hits_removed_num_chip_sum + chromosome[i].redundant_hits_removed_num_chip;		
		redundant_hits_num_chip_sum   = redundant_hits_num_chip_sum + chromosome[i].redundant_hits_num_chip;		
		f_hit_ctrl_num_sum = f_hit_ctrl_num_sum + chromosome[i].f_hit_num_ctrl;
		r_hit_ctrl_num_sum = r_hit_ctrl_num_sum + chromosome[i].r_hit_num_ctrl;
		hit_ctrl_num_sum   = hit_ctrl_num_sum   + chromosome[i].hit_num_ctrl;
		redundant_hits_removed_num_ctrl_sum   = redundant_hits_removed_num_ctrl_sum + chromosome[i].redundant_hits_removed_num_ctrl;
		redundant_hits_num_ctrl_sum   = redundant_hits_num_ctrl_sum + chromosome[i].redundant_hits_num_ctrl;		
		qfrag_num_sum      = qfrag_num_sum      + chromosome[i].qfrag_num_chip;
	}
	
	// calculate duplication rates for all chromosomes
	float dup_rate_chip=0;
	if(0<hit_chip_num_sum)
	{
		dup_rate_chip=(float)redundant_hits_num_chip_sum/(redundant_hits_removed_num_chip_sum+hit_chip_num_sum);
	}
	
	float dup_rate_ctrl=0;
	if(0<hit_ctrl_num_sum)
	{
		dup_rate_ctrl=(float)redundant_hits_num_ctrl_sum/(redundant_hits_removed_num_ctrl_sum+hit_ctrl_num_sum);
	}
	
	// calculate expected number of qfrags
	double qfrags_num_exp=0;
	if(0<f_hit_chip_num_sum && 0<r_hit_chip_num_sum)
	{
		qfrags_num_exp=f_hit_chip_num_sum * (q_max-q_min) * (double)r_hit_chip_num_sum/len_sum;
	}
	double q_enrichment_score=0;
	if(0<qfrags_num_exp)
	{
		q_enrichment_score=qfrag_num_sum/qfrags_num_exp;
	}
	if(keep_dup)
	{
		q_enrichment_score=-1;
	}

	OUT.precision(2);
	OUT	<< "ALL"                                    << "\t"
		<< len_sum                                  << "\t"
		<< sum_num_sum                              << "\t"		
		<< f_hit_chip_num_sum                       << "\t"
		<< r_hit_chip_num_sum                       << "\t"
		<< hit_chip_num_sum                         << "\t"		
		<< redundant_hits_removed_num_chip_sum      << "\t"
		<< dup_rate_chip                            << "\t"		
		<< f_hit_ctrl_num_sum                       << "\t"	
		<< r_hit_ctrl_num_sum                       << "\t"		
		<< hit_ctrl_num_sum                         << "\t"
		<< redundant_hits_removed_num_ctrl_sum      << "\t"
		<< dup_rate_ctrl                            << "\t"
		<< std::fixed << (int)round(qfrags_num_exp) << "\t"
		<< qfrag_num_sum                            << "\t"					
		<< q_enrichment_score                       << "\n\n";

		
	OUT.close();
	return 0;
}

///////////////////////////////////////////////////////////////////////

int writeChromosomeInfoShort(seqan::CharString out_prefix, seqan::CharString chromosome_info_file, std::vector<Chromosome> &chromosome, int chr_num, int q_min, int q_max, bool keep_dup, int fragment_length_avg)
{
	std::ofstream OUT;
	OUT.open (toCString(chromosome_info_file));	

	OUT	<< "chromosome"                            << "\t"
		<< "# hits-treatment"                      << "\t"
		<< "# redundant-hits-removed-treatment"    << "\t"
		<< "duplication-rate-treatment"            << "\t"
		<< "# hits-control"                        << "\t"
		<< "# redundant-hits-removed-control"      << "\t"			
		<< "duplication-rate-control"              << "\t"
		<< "# q-hits-treatment-expected"           << "\t"
		<< "# q-hits-treatment-observed"           << "\t"		
		<< "q-enrichment-score (QES)"              << "\t"
		<< "q-sample-score (QSS)"                  << "\t"
		<< "estimated number signal hits (ENSH)"   << "\t"		
		<< "average fragment length"               << "\t"
		<< "relative strand correlation (RSC)"           << "\n";				
		

	long int len_sum                          = 0;
	int sum_num_sum                           = 0;
	int f_hit_chip_num_sum                    = 0;
	int r_hit_chip_num_sum                    = 0;
	int hit_chip_num_sum                      = 0;
	int redundant_hits_removed_num_chip_sum   = 0;	
	int redundant_hits_num_chip_sum           = 0;	
	long int f_hit_ctrl_num_sum                    = 0;
	long int r_hit_ctrl_num_sum                    = 0;	
	int hit_ctrl_num_sum                      = 0;	
	int redundant_hits_removed_num_ctrl_sum   = 0;
	int redundant_hits_num_ctrl_sum           = 0;		
	int qfrag_num_sum                         = 0;
	int obs_q_hit_num_sum                         = 0;
	
	// sum up infos for all chromosomes
	for(int i=0;i<chr_num;i++)
	{
		// calculate duplication rates
		float dup_rate_chip=0;
		if(0<chromosome[i].hit_num_chip)
		{
			dup_rate_chip=(float)chromosome[i].redundant_hits_num_chip/(chromosome[i].hit_num_chip+chromosome[i].redundant_hits_removed_num_chip);
		}

		float dup_rate_ctrl=0;
		if(0<chromosome[i].hit_num_ctrl)
		{
			dup_rate_ctrl=(float)chromosome[i].redundant_hits_num_ctrl/(chromosome[i].hit_num_ctrl+chromosome[i].redundant_hits_removed_num_ctrl);
		}
		
		double exp_q_hit_num =-1;
		int obs_q_hit_num =-1;
		double QES =-1;
		double QSS = -1;
		int ENSH = -1;
		if(!keep_dup)
		{
			// calculate expected number of qhits
			exp_q_hit_num=getExpNumQHits(
				chromosome[i].f_hit_num_chip,
				chromosome[i].r_hit_num_chip,
				chromosome[i].len,
				q_min,q_max);

			// get observed number of qhits
			obs_q_hit_num=getObsNumQHits(chromosome[i].CHIP_HITS);
		
			// get q-enrichment-score
			QES = getQEnrichmentScore(
				obs_q_hit_num,
				exp_q_hit_num,
				chromosome[i].f_hit_num_chip+chromosome[i].r_hit_num_chip);
		
			// get q-sample score
			QSS = QES * (1-dup_rate_chip);
		
			// get estimated number of signal hits
			ENSH = round(chromosome[i].hit_num_chip * QSS);
		}
			
		OUT.precision(2);
		OUT.setf(std::ios_base::fixed);	
		OUT	<< chromosome[i].name                                  << "\t"
			<< chromosome[i].hit_num_chip                          << "\t"
			<< chromosome[i].redundant_hits_removed_num_chip       << "\t"
			<< dup_rate_chip                                       << "\t"
			<< chromosome[i].hit_num_ctrl                          << "\t"			
			<< chromosome[i].redundant_hits_removed_num_ctrl       << "\t"
			<< dup_rate_ctrl                                       << "\t"
			<< (int)round(exp_q_hit_num)                           << "\t"
			<< obs_q_hit_num                                       << "\t"
			<< QES                                                 << "\t"
			<< QSS                                                 << "\t"
			<< ENSH                                                << "\t"							
			<< "NA"                                                << "\t"
			<< "NA"                                                << "\n";
			
		len_sum            = len_sum            + chromosome[i].len;
		sum_num_sum        = sum_num_sum        + chromosome[i].sum_num;		
		f_hit_chip_num_sum = f_hit_chip_num_sum + chromosome[i].f_hit_num_chip;
		r_hit_chip_num_sum = r_hit_chip_num_sum + chromosome[i].r_hit_num_chip;
		hit_chip_num_sum   = hit_chip_num_sum   + chromosome[i].hit_num_chip;
		redundant_hits_removed_num_chip_sum   = redundant_hits_removed_num_chip_sum + chromosome[i].redundant_hits_removed_num_chip;		
		redundant_hits_num_chip_sum   = redundant_hits_num_chip_sum + chromosome[i].redundant_hits_num_chip;		
		f_hit_ctrl_num_sum = f_hit_ctrl_num_sum + chromosome[i].f_hit_num_ctrl;
		r_hit_ctrl_num_sum = r_hit_ctrl_num_sum + chromosome[i].r_hit_num_ctrl;
		hit_ctrl_num_sum   = hit_ctrl_num_sum   + chromosome[i].hit_num_ctrl;
		redundant_hits_removed_num_ctrl_sum   = redundant_hits_removed_num_ctrl_sum + chromosome[i].redundant_hits_removed_num_ctrl;
		redundant_hits_num_ctrl_sum   = redundant_hits_num_ctrl_sum + chromosome[i].redundant_hits_num_ctrl;		
		qfrag_num_sum      = qfrag_num_sum      + chromosome[i].qfrag_num_chip;
		obs_q_hit_num_sum = obs_q_hit_num_sum + obs_q_hit_num;
	}
	
	// calculate duplication rates for all chromosomes
	float dup_rate_chip=0;
	if(0<hit_chip_num_sum)
	{
		dup_rate_chip=(float)redundant_hits_num_chip_sum/(redundant_hits_removed_num_chip_sum+hit_chip_num_sum);
	}
	
	float dup_rate_ctrl=0;
	if(0<hit_ctrl_num_sum)
	{
		dup_rate_ctrl=(float)redundant_hits_num_ctrl_sum/(redundant_hits_removed_num_ctrl_sum+hit_ctrl_num_sum);
	}
	
	double exp_q_hit_num =-1;
	double QES =-1;
	double QSS = -1;
	int ENSH = -1;
	if(!keep_dup)
	{
		// calculate expected number of qhits
		exp_q_hit_num=getExpNumQHits(
			f_hit_chip_num_sum,
			r_hit_chip_num_sum,
			len_sum,
			q_min,q_max);
				
		// get q-enrichment-score
		QES = getQEnrichmentScore(
			obs_q_hit_num_sum,
			exp_q_hit_num,
			f_hit_chip_num_sum+r_hit_chip_num_sum);

		// get q-sample score
		QSS = QES * (1-dup_rate_chip);

		// get estimated number of signal hits
		ENSH = round(hit_chip_num_sum * QES);
	}
	else
	{
		obs_q_hit_num_sum=-1;
	}			

	OUT.precision(2);
	OUT	<< out_prefix                               << "\t"
		<< hit_chip_num_sum                         << "\t"		
		<< redundant_hits_removed_num_chip_sum      << "\t"
		<< dup_rate_chip                            << "\t"		
		<< hit_ctrl_num_sum                         << "\t"
		<< redundant_hits_removed_num_ctrl_sum      << "\t"
		<< dup_rate_ctrl                            << "\t"
		<< (int)round(exp_q_hit_num)                << "\t"
		<< obs_q_hit_num_sum                        << "\t"
		<< QES                                      << "\t"
		<< QSS                                      << "\t"	
		<< ENSH                                     << "\t"				
		<< fragment_length_avg                      << "\t"
		<< chromosome[0].rsc			            << "\n";		
		
	OUT.close();
	return 0;
}

long double getExpNumQHits(long int TF, long int TR, long int L, int q_min, int q_max)
{
	long double lambda_t = (q_max-q_min) * (long double)TR/L;
	long double exp_q_hit_num = (TF+TR) * (1-expl(-lambda_t));
	return exp_q_hit_num;
}

double getQEnrichmentScore(int obs_q_hit_num, double exp_q_hit_num, int UMNR_hit_num)
{
	double QES = (obs_q_hit_num-exp_q_hit_num)/UMNR_hit_num;
	if(QES<0){return 0;}
	else{return QES;}
}

int getObsNumQHits(std::vector<Hit> HITS)
{
	int obs_q_hit_num=0;
	for(unsigned int i=0;i<HITS.size();i++)
	{
		if(HITS[i].is_q_hit)
		{
			obs_q_hit_num++;
		}
	}
	return obs_q_hit_num;
}

///////////////////////////////////////////////////////////////////////

int writeSummitInfo(seqan::CharString summit_info_file, std::vector<Summit_chr> &summit, Options &options)
{
	int sum_num=summit.size();
	
	std::ofstream OUT;
	OUT.open (toCString(summit_info_file));
    
	OUT	<< "\"Chromosome\""        << "\t"
		<< "\"pos\""               << "\t"
		<< "\"pos+1\""             << "\t"
		<< "\"q_cov_chip\""        << "\t"
		<< "\"q_cov_ctrl\""        << "\t"
		<< "\"q_beg_chip\""        << "\t"
		<< "\"q_end_chip\""        << "\t"
		<< "\"q_beg_ctrl\""        << "\t"					
		<< "\"q_end_ctrl\""        << "\t"	
		<< "\"saturation-score\""  << "\t"
		<< "\"p-value\""           << "\t"	
		<< "\"q-value\""           << "\n";				
    
	for(int i=0;i<sum_num;i++)
	{
		if(options.p_value_cutoff != -1)
		{
			if(-log10l(summit[i].summit.q_value)<=options.p_value_cutoff)
			{
				break;
			}
		}	
		else if(i==options.top_n)
		{
			break;
		}
		
		OUT	<< summit[i].chr_name                        << "\t"
			<< summit[i].summit.pos              << "\t"			
			<< summit[i].summit.pos+1            << "\t"
			<< summit[i].summit.q_cov_chip       << "\t"
			<< summit[i].summit.q_cov_ctrl       << "\t"
			<< summit[i].summit.q_beg_chip       << "\t"				
			<< summit[i].summit.q_end_chip       << "\t"
			<< summit[i].summit.q_beg_ctrl       << "\t"				
			<< summit[i].summit.q_end_ctrl       << "\t"				
			<< summit[i].summit.saturation_score << "\t"
			<< -log10l(summit[i].summit.p_value) << "\t"
			<< -log10l(summit[i].summit.q_value) << "\n";				
	}
	OUT.close();

	return 0;
}


int writeNarrowPeak(seqan::CharString narrowPeak_file, std::vector<Summit_chr> &summit, int radius, Options &options)
{
	int sum_num=summit.size();

	std::ofstream OUT;
	OUT.open (toCString(narrowPeak_file));

	for(int i=0;i<sum_num;i++)
	{
		if(options.p_value_cutoff != -1)
		{
			if(-log10l(summit[i].summit.q_value)<=options.p_value_cutoff)
			{
				break;
			}
		}	
		else if(i==options.top_n)
		{
			break;
		}
		
		int sta;
		if(summit[i].summit.pos-radius<0)
		{
			sta=0;
		}
		else
		{
			sta=summit[i].summit.pos-radius;
		}
		int end;
		if(summit[i].chr_length<summit[i].summit.pos+radius)
		{
			end=summit[i].chr_length;
		}
		else
		{
			end=summit[i].summit.pos+radius;
		}
		OUT.setf(std::ios_base::fixed);		
		OUT	<< summit[i].chr_name                                 << "\t"
			<< sta                                                << "\t"			
			<< end                                                << "\t"
			<< "."                                                << "\t"
			<< int(round(summit[i].summit.saturation_score*1000)) << "\t"
			<< "."                                                << "\t"
			<< summit[i].summit.saturation_score                  << "\t"
			<< -log10l(summit[i].summit.p_value)                  << "\t"
			<< -log10l(summit[i].summit.q_value)                  << "\t"
			<< radius                                             << "\n";
	}
	OUT.close();
		
	return 0;
}

int writeReadsToBED(Chromosome &chromosome, std::string out_prefix, int extension_length, int shift_size)
{
	// fiddle seqan::CharString to std::string
	char* foo = toCString(chromosome.name);
	std::string chr_name = foo;
	
	std::ofstream OUT;
	std::string chip_out = out_prefix + "-" + chr_name + "-chip.bed";
	OUT.open(chip_out.c_str());
	
	for(int j=0;j<chromosome.hit_num_chip;j++)
	{	
		int sta,end;
		if(chromosome.CHIP_HITS[j].strand == 0)
		{
			sta=chromosome.CHIP_HITS[j].pos+shift_size;
			end=sta+extension_length;
			
		}
		else
		{
			end=chromosome.CHIP_HITS[j].pos-shift_size+1;
			sta=end-extension_length;
		}
		if(chromosome.len<sta || end<0){continue;}
		if(sta<0){sta=0;}
		if(chromosome.len<end){end=chromosome.len;}
		OUT << chromosome.name << "\t";
		OUT << sta << "\t";
		OUT << end << "\n";
	}
	OUT.close();
	
			
			
	std::string control_out = out_prefix + "-" + chr_name + "-control.bed";
	OUT.open (control_out.c_str());
	for(int j=0;j<chromosome.hit_num_ctrl;j++)
	{
		int sta,end;
		if(chromosome.CTRL_HITS[j].strand == 0)
		{
			sta=chromosome.CTRL_HITS[j].pos+shift_size;
			end=sta+extension_length;
			
		}
		else
		{
			end=chromosome.CTRL_HITS[j].pos-shift_size+1;
			sta=end-extension_length;
		}
		if(chromosome.len<sta || end<0){continue;}
		if(sta<0){sta=0;}
		if(chromosome.len<end){end=chromosome.len;}
		OUT << chromosome.name << "\t";
		OUT << sta << "\t";
		OUT << end << "\n";
	}
	OUT.close();	
	return 0;
}

int writeQFragsToBED(Chromosome &chromosome, std::string out_prefix, int min, int max)
{
	// fiddle seqan::CharString to std::string
	char* foo = toCString(chromosome.name);
	std::string chr_name = foo;
	
	std::ofstream OUT;
	std::string chip_out = out_prefix + "-" + chr_name + "-chip.bed";
	OUT.open (chip_out.c_str());
	
	int c_hit=0;
	int pos=0;
	while(pos<chromosome.len)
	{
		// at position pos are one or more hits
		while(c_hit<chromosome.hit_num_chip && pos==chromosome.CHIP_HITS[c_hit].pos)
		{
			// current hit is on the forward strand
			if(chromosome.CHIP_HITS[c_hit].strand==0)
			{
				// look for qfrags
				int j=c_hit+1; // start examination with the next hit
				
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_chip && chromosome.CHIP_HITS[j].pos <= (chromosome.CHIP_HITS[c_hit].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CHIP_HITS[j].strand==1 && chromosome.CHIP_HITS[c_hit].pos+min<chromosome.CHIP_HITS[j].pos)
					{
						// write qfrags
						OUT << chromosome.name << "\t" << chromosome.CHIP_HITS[c_hit].pos << "\t" << chromosome.CHIP_HITS[j].pos+1 << "\n";  
					}
					j++;
				}
			}
			c_hit++;			
		}
		pos++;
	}
	
	OUT.close();
	
	std::string control_out = out_prefix + "-" + chr_name + "-control.bed";
	OUT.open (control_out.c_str());
	
	c_hit=0;
	pos=0;
	while(pos<chromosome.len)
	{
		// at position pos are one or more hits
		while(c_hit<chromosome.hit_num_ctrl && pos==chromosome.CTRL_HITS[c_hit].pos)
		{
			// current hit is on the forward strand
			if(chromosome.CTRL_HITS[c_hit].strand==0)
			{
				// look for qfrags
				int j=c_hit+1; // start examination with the next hit
				
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_ctrl && chromosome.CTRL_HITS[j].pos <= (chromosome.CTRL_HITS[c_hit].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CTRL_HITS[j].strand==1 && chromosome.CTRL_HITS[c_hit].pos+min<chromosome.CTRL_HITS[j].pos)
					{
						// write qfrags
						OUT << chromosome.name << "\t" << chromosome.CTRL_HITS[c_hit].pos << "\t" << chromosome.CTRL_HITS[j].pos << "\n";  
					}
					j++;
				}
			}
			c_hit++;			
		}
		pos++;
	}
	
	OUT.close();
	
	return 0;
}


int writePlainShiftedExtendedReadsAndQfragsToBED(std::vector<Chromosome> &chromosome, Options &options)
{
	// Find the array index for chromosome ID
	int chr_idx=-1;
	for(unsigned int i=0;i<chromosome.size();i++){if(chromosome[i].name==options.write_bed){chr_idx=i;}}
	if(chr_idx==-1)
	{
		std::cout << "No matches for " << options.write_bed << "\n";
	}
	
	// Calulate average fragment length from lower and upper bound
	int fl=options.lowerbound+((options.upperbound-options.lowerbound)/2);
	if(options.verbose)
	{
		std::cout << "\tAverage fragment length: " << fl << "\n";
	}
	
	// Get read length from first alignment in ChIP BAM file
	seqan::BamStream bamStreamInChIP(toCString(options.chip_sample));
	seqan::BamAlignmentRecord record;
	readRecord(record, bamStreamInChIP);
	int rl=length(record.seq);
	if(options.verbose)
	{
		std::cout << "\tRead length: " << rl << "\n";
	}
	
	// Create names for for output files from outprefix
	seqan::CharString out_matched_reads=options.out_prefix;
	append(out_matched_reads,"-matched-reads");
	seqan::CharString out_shifted_reads=options.out_prefix;
	append(out_shifted_reads,"-shifted-reads");
	seqan::CharString out_extended_reads=options.out_prefix;
	append(out_extended_reads,"-extended-reads");
	seqan::CharString out_qfrags=options.out_prefix;
	append(out_qfrags,"-qfrags");

	// Call functions for writing the BED files
	writeReadsToBED(chromosome[chr_idx],toCString(out_matched_reads), rl, 0);
	writeReadsToBED(chromosome[chr_idx],toCString(out_shifted_reads), rl, fl/2);
	writeReadsToBED(chromosome[chr_idx],toCString(out_extended_reads), fl, 0);
	writeQFragsToBED(chromosome[chr_idx],toCString(out_qfrags),options.lowerbound,options.upperbound);	
	
	return 0;
}

int writeHitDist(std::vector<Chromosome> &chromosome, int radius, std::string out_prefix)
{	
	// Prepare vectors for counting
	// ----------------------------
	
	std::vector<int> pos;
	pos.resize(2*radius+1);
	for(int i=0;i<(2*radius+1);i++){pos[i]=i-radius;}
	
	std::vector<int> fwd_chip;
	fwd_chip.resize(2*radius+1);
	
	std::vector<int> rev_chip;
	rev_chip.resize(2*radius+1);
	
	std::vector<int> fwd_ctrl;
	fwd_ctrl.resize(2*radius+1);
	
	std::vector<int> rev_ctrl;
	rev_ctrl.resize(2*radius+1);

	
	// Count
	// -----
		
	for(unsigned int i=0; i<chromosome.size();i++)
	{
		// for each summit count the hits relative to summit position
		for(int j=2;j<chromosome[i].sum_num;j++)
		{
			if(chromosome[i].len<=chromosome[i].SUMMITS[j].pos){continue;}

			// jump to the first hit after the summit
			int hit_idx=chromosome[i].SUMMITS[j].anchor_hit_chip;
			
			// go back to the first hit before sum_pos-radius
			while((chromosome[i].SUMMITS[j].pos-radius-1)<chromosome[i].CHIP_HITS[hit_idx].pos)
			{
				hit_idx--;
				if(hit_idx<0) break;
			}
			hit_idx++;
			while(chromosome[i].CHIP_HITS[hit_idx].pos<=(chromosome[i].SUMMITS[j].pos+radius))
			{
				int position=chromosome[i].CHIP_HITS[hit_idx].pos-chromosome[i].SUMMITS[j].pos+radius;
				if(chromosome[i].CHIP_HITS[hit_idx].strand==0)
				{
					fwd_chip[position]++;
				}
				else
				{
					rev_chip[position]++;
				}
				hit_idx++;
				if(chromosome[i].hit_num_chip<=hit_idx) break;
			}

			if(chromosome[i].hit_num_ctrl==0) continue;


			// do the same for the control
			// ---------------------------
			
			// jump to the first hit after the summit
			hit_idx=chromosome[i].SUMMITS[j].anchor_hit_ctrl;
			
			// go back to the first hit before sum_pos-radius
			while((chromosome[i].SUMMITS[j].pos-radius-1)<chromosome[i].CTRL_HITS[hit_idx].pos)
			{
				hit_idx--;
				if(hit_idx<0) break;
			}
			hit_idx++;
			while(chromosome[i].CTRL_HITS[hit_idx].pos<=(chromosome[i].SUMMITS[j].pos+radius))
			{
				int position=chromosome[i].CTRL_HITS[hit_idx].pos-chromosome[i].SUMMITS[j].pos+radius;
				if(chromosome[i].CTRL_HITS[hit_idx].strand==0)
				{
					fwd_ctrl[position]++;
				}
				else
				{
					rev_ctrl[position]++;
				}
				hit_idx++;
				if(chromosome[i].hit_num_ctrl<=hit_idx) break;
			}
		}
	}
	
	// Write distribution to file
	// --------------------------
	
	std::ofstream OUT;
	std::string chip_out = out_prefix + "-hit-dist.tab";
	OUT.open(chip_out.c_str());
	
	for(int i=0;i<2*radius+1;i++)
	{
		OUT << pos[i] << "\t";
		OUT << fwd_chip[i] << "\t" << rev_chip[i] << "\t";
		OUT << fwd_ctrl[i] << "\t" << rev_ctrl[i] << "\n";		
	}
	OUT.close();
	return 0;
}



int writeBedGraph(std::vector<Chromosome> &chromosome, seqan::CharString out_prefix, int extension_length, bool is_control)
{
	// This function takes a vector of chromosomes and writes a bedGraph file.
	// If 'is_control' is false this function will use the treatment hits,
	// otherwise it will use the control hits.
	// Hits will be extended by 'extension_length'.
	// If you want the bedGraph for read coverage use the read length as extension length.
	// If you want the bedGraph for fragment coverage use the average fragment length.
	// fiddle seqan::CharString to std::string
	

	// create file names and out stream
	char* foo = toCString(out_prefix);
	std::string op = foo;	
	std::ofstream OUT;
	std::string out;
	if(is_control)
	{
		out = op + "-Q-control.bedgraph";
	}
	else
	{
		out = op + "-Q-treatment.bedgraph";
	}
	OUT.open (out.c_str());
	
	// iterate over all chromosomes
	for(unsigned int i=0; i<chromosome.size(); i++)
	{
		std::vector<Hit> HITS;
		if(is_control)
		{
			HITS=chromosome[i].CTRL_HITS;
		}
		else
		{
			HITS=chromosome[i].CHIP_HITS;
		}
		
		// init coverage vector
		std::vector<int> COV(chromosome[i].len,0);

		for(unsigned j=0;j<HITS.size();j++)
		{
			if(HITS[j].strand == 0)
			{
				for(int k=0;k<extension_length;k++)
				{
					if((int)COV.size()<(HITS[j].pos+k+1)) continue;
					COV[HITS[j].pos+k]++;
				}
			}
			else
			{
				for(int k=0;k<extension_length;k++)
				{
					if(((int)HITS[j].pos-k)<0) continue;
					COV[HITS[j].pos-k]++;
				}				
			}
		}		
		
		int c_pos=0;
		int c_cov=COV[0];		
		for(unsigned j=1;j<COV.size();j++)
		{
			if(COV[j]!=COV[j-1])
			{
				OUT << chromosome[i].name << "\t" << c_pos << "\t" << j << "\t" << c_cov << "\n";
				c_pos=j;
				c_cov=COV[j];
			}
		}
	}
	OUT.close();
	return 0;
}
