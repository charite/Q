#include "Q_EvaluateSignalsCount.h"
#include "Q_EvaluateSignalsPvalue.h"
#include "Q_ReadInFiles.h"
#include <boost/unordered_map.hpp>


int getQFragsCoverage(Chromosome &chromosome, int min, int max)
{
			//std::cout << chromosome.name << "start \n";
			
	// for each summit
	for(int i=0;i<chromosome.sum_num;i++)
	{
		// start with ChIP sample
		// ----------------------
		
		// jump to the first hit after the summit
		int hit_idx=chromosome.SUMMITS[i].anchor_hit_chip;
		
		// go back to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-max)<chromosome.CHIP_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;

		// now derive qfrags and count ends
		while(hit_idx<chromosome.hit_num_chip && chromosome.CHIP_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)		
		{
			// current hit is on the forward strand
			if(chromosome.CHIP_HITS[hit_idx].strand==0)
			{
				// look for qfrags
				int j=hit_idx+1; // start examination with the next hit
				
			
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_chip && chromosome.CHIP_HITS[j].pos <= (chromosome.CHIP_HITS[hit_idx].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CHIP_HITS[j].strand==1 && chromosome.CHIP_HITS[hit_idx].pos+min<chromosome.CHIP_HITS[j].pos)
					{
						// second hit is beyond pos
						if(chromosome.SUMMITS[i].pos<=chromosome.CHIP_HITS[j].pos)
						{
							chromosome.SUMMITS[i].q_cov_chip++;
						}
					}
					j++;
				}
				hit_idx++;			
			}
			else{hit_idx++;}
			if(chromosome.hit_num_chip<hit_idx){break;}
		}
		
		// the same for control
		// --------------------
		
		// if there is no control continue
		if(chromosome.hit_num_ctrl==0){continue;}
		
		// jump to the first hit after the summit
		hit_idx=chromosome.SUMMITS[i].anchor_hit_ctrl;
		
		// go back to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-max)<chromosome.CTRL_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;

		// now derive qfrags and count ends
		while(hit_idx<chromosome.hit_num_ctrl && chromosome.CTRL_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)
		{
			// current hit is on the forward strand
			if(chromosome.CTRL_HITS[hit_idx].strand==0)
			{
				// look for qfrags
				int j=hit_idx+1; // start examination with the next hit
				
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_ctrl && chromosome.CTRL_HITS[j].pos <= (chromosome.CTRL_HITS[hit_idx].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CTRL_HITS[j].strand==1 && chromosome.CTRL_HITS[hit_idx].pos+min<chromosome.CTRL_HITS[j].pos)
					{
						// second hit is beyond pos
						if(chromosome.SUMMITS[i].pos<chromosome.CTRL_HITS[j].pos)
						{
							chromosome.SUMMITS[i].q_cov_ctrl++;
							// write qfrags
						}
					}
					j++;
				}
				hit_idx++;			
			}
			else{hit_idx++;}
			if(chromosome.hit_num_ctrl<hit_idx){break;}
		}
	}

	// remove all summits with q_cov<1
	std::vector<Summit> S;
	for(int i=0;i<chromosome.sum_num;i++)
	{
		if(1<=chromosome.SUMMITS[i].q_cov_chip)
		{
			S.push_back(chromosome.SUMMITS[i]);
		}
	}
	chromosome.SUMMITS=S;
	chromosome.sum_num=S.size();
	
	return 0;
}


int getQfragEnds(Chromosome &chromosome, int min, int max)
{
	for(int i=0;i<chromosome.sum_num;i++)
	{
		int hit_idx=chromosome.SUMMITS[i].anchor_hit_chip;
		
		// go to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-max)<chromosome.CHIP_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;

		// now derive qfrags and count ends
		boost::unordered_map <int, int> begin; // key and value is hit index
		boost::unordered_map <int, int> end;
		while(hit_idx<chromosome.hit_num_chip && chromosome.CHIP_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)
		{
			// current hit is on the forward strand
			if(chromosome.CHIP_HITS[hit_idx].strand==0)
			{
				// look for qfrags
				int j=hit_idx+1; // start examination with the next hit
				
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_chip && chromosome.CHIP_HITS[j].pos < (chromosome.CHIP_HITS[hit_idx].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CHIP_HITS[j].strand==1 && chromosome.CHIP_HITS[hit_idx].pos+min<chromosome.CHIP_HITS[j].pos)
					{
						// second hit is beyond pos
						if(chromosome.SUMMITS[i].pos<chromosome.CHIP_HITS[j].pos)
						{
							begin[chromosome.CHIP_HITS[hit_idx].pos]=hit_idx;
							end[chromosome.CHIP_HITS[j].pos]=j;
						}
					}
					j++;
				}
				hit_idx++;			
			}
			else{hit_idx++;}
			if(chromosome.hit_num_chip<hit_idx){break;}
		}
		chromosome.SUMMITS[i].q_beg_chip=begin.size();
		chromosome.SUMMITS[i].q_end_chip=end.size();
	}
	
	// remove all summits that have less 
	// saturated positions than expected by chance
	double p_t=get_probability_of_success(chromosome.f_hit_num_chip,chromosome.r_hit_num_chip,chromosome.len,min,max);
	int exp_sat_pos=round(2*max*p_t);
	std::vector<Summit> S;
	for(int i=0;i<chromosome.sum_num;i++)
	{
		if(chromosome.SUMMITS[i].q_beg_chip>0 && chromosome.SUMMITS[i].q_end_chip>0
			&& (exp_sat_pos<=chromosome.SUMMITS[i].q_beg_chip+chromosome.SUMMITS[i].q_end_chip))
		{
			S.push_back(chromosome.SUMMITS[i]);
		}
	}
	//chromosome.SUMMITS=S;
	//chromosome.sum_num=S.size();	
	
	// if there is no control, return
	if(chromosome.hit_num_ctrl==0){return 0;}
	
	// else do it again for the control
	for(int i=0;i<chromosome.sum_num;i++)
	{
		int hit_idx=chromosome.SUMMITS[i].anchor_hit_ctrl;
		
		// go to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-max)<chromosome.CTRL_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;

		// now derive qfrags and count ends
		boost::unordered_map <int, int> begin; // key and value is hit index
		boost::unordered_map <int, int> end;
		while(hit_idx<chromosome.hit_num_ctrl && chromosome.CTRL_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)
		{
			// current hit is on the forward strand
			if(chromosome.CTRL_HITS[hit_idx].strand==0)
			{
				// look for qfrags
				int j=hit_idx+1; // start examination with the next hit
				
				// second hit for qfrag is at most max bases apart
				while(j<chromosome.hit_num_ctrl && chromosome.CTRL_HITS[j].pos < (chromosome.CTRL_HITS[hit_idx].pos+max))
				{
					// second hit of the qfrag is at least min bases apart
					// and on the reverse strand
					if(chromosome.CTRL_HITS[j].strand==1 && chromosome.CTRL_HITS[hit_idx].pos+min<chromosome.CTRL_HITS[j].pos)
					{
						// second hit is beyond pos
						if(chromosome.SUMMITS[i].pos<chromosome.CTRL_HITS[j].pos)
						{
							begin[chromosome.CTRL_HITS[hit_idx].pos]=hit_idx;
							end[chromosome.CTRL_HITS[j].pos]=j;
						}
					}
					j++;
				}
				hit_idx++;			
			}
			else{hit_idx++;}
			if(chromosome.hit_num_ctrl<hit_idx){break;}
		}
		chromosome.SUMMITS[i].q_beg_ctrl=begin.size();
		chromosome.SUMMITS[i].q_end_ctrl=end.size();
	}
	
	return 0;
}

int getKs(Chromosome &chromosome, int radius)
{
	for(int i=0;i<chromosome.sum_num;i++)
	{
		// jump to the first hit after the summit
		int hit_idx=chromosome.SUMMITS[i].anchor_hit_chip;
		
		// go back to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-radius)<chromosome.CHIP_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;
		while(chromosome.CHIP_HITS[hit_idx].pos<(chromosome.SUMMITS[i].pos+radius))
		{
			if(chromosome.CHIP_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)
			{
				if(chromosome.CHIP_HITS[hit_idx].strand==0)
				{
					chromosome.SUMMITS[i].kfu_chip++;
				}
				else
				{
					chromosome.SUMMITS[i].kru_chip++;
				}
			}
			else
			{
				if(chromosome.CHIP_HITS[hit_idx].strand==0)
				{
					chromosome.SUMMITS[i].kfd_chip++;
				}
				else
				{
					chromosome.SUMMITS[i].krd_chip++;
				}
			}
			hit_idx++;
			if(chromosome.hit_num_chip<hit_idx) break;
		}
		chromosome.SUMMITS[i].q_cov_chip=
			chromosome.SUMMITS[i].kfu_chip+
			chromosome.SUMMITS[i].krd_chip+
			chromosome.SUMMITS[i].kfd_chip+
			chromosome.SUMMITS[i].kru_chip;
	}

	
	// if there is no control, return
	if(chromosome.hit_num_ctrl==0){return 0;}
	
	// else do it again for the control
	for(int i=0;i<chromosome.sum_num;i++)
	{
		// jump to the first hit after the summit
		int hit_idx=chromosome.SUMMITS[i].anchor_hit_ctrl;
		
		// go back to the first hit before sum_pos-max
		while((chromosome.SUMMITS[i].pos-radius)<chromosome.CTRL_HITS[hit_idx].pos)
		{
			hit_idx--;
			if(hit_idx<0) break;
		}
		hit_idx++;

		while(chromosome.CTRL_HITS[hit_idx].pos<(chromosome.SUMMITS[i].pos+radius))
		{
			if(chromosome.CTRL_HITS[hit_idx].pos<chromosome.SUMMITS[i].pos)
			{
				if(chromosome.CTRL_HITS[hit_idx].strand==0)
				{
					chromosome.SUMMITS[i].kfu_ctrl++;
				}
				else
				{
					chromosome.SUMMITS[i].kru_ctrl++;
				}
			}
			else
			{
				if(chromosome.CTRL_HITS[hit_idx].strand==0)
				{
					chromosome.SUMMITS[i].kfd_ctrl++;
				}
				else
				{
					chromosome.SUMMITS[i].krd_ctrl++;
				}
			}
			hit_idx++;
			if(chromosome.hit_num_ctrl<hit_idx) break;
		}
	}	
	return 0;
}
