#include "Q_ReadInFiles.h"
#include "Q_GetFragLength.h"
#include <boost/dynamic_bitset.hpp>

int getFragLength(std::vector<Chromosome> &chromosome, int step_num, int thread_num, std::string out_prefix)
{
	// get number of chromosomes
	int chr_num=chromosome.size();
	
	// get read length
	int read_len=chromosome[0].read_len_chip;
	
	SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
	for(int i=0;i<chr_num;i++)
	{
		// init bitvectors: a position is set if covered by at least one hit
		int L = chromosome[i].len;
		boost::dynamic_bitset<> fwd(L);
		boost::dynamic_bitset<> rev(L);
		for(int j=0;j<chromosome[i].hit_num_chip;j++)
		{
			if(chromosome[i].CHIP_HITS[j].pos>=chromosome[i].len) continue;

			if(chromosome[i].CHIP_HITS[j].strand==0)
			{
				fwd.set(L-1-chromosome[i].CHIP_HITS[j].pos);
			}
			else
			{
				rev.set(L-1-chromosome[i].CHIP_HITS[j].pos);
			}
		}

		// calculate hamming distances and shift
		for(int step=0;step<step_num;step++)
		{
			(fwd >>= 1);
			chromosome[i].HAMMING_DISTANCES[step]=((fwd ^ rev).count());
		}
	}
	
	// combine hamming distances for all chromosomes
	std::map<int, int> GlobalHammingDistances;
	for(int step=0;step<step_num;step++)
	{
		int sum=0;
		int len_sum=0;
		for(int j=0;j<chr_num;j++)
		{
			sum=sum+chromosome[j].HAMMING_DISTANCES[step];
			len_sum=len_sum+chromosome[j].len;
		}
		GlobalHammingDistances[step]=sum;
	}
/*	
	// get fragment length
	int min_hd=GlobalHammingDistances[0];//put here 70
	int max_hd=GlobalHammingDistances[0];//put here 70
	int fl=0;
	int rl=0;
	for(int step=1;step<step_num;step++)
	{
		if(GlobalHammingDistances[step]<min_hd && (step<read_len))
		{
			min_hd=GlobalHammingDistances[step];
			rl=step;
		}
		
		if(step==read_len)
		{
			min_hd=GlobalHammingDistances[0];
		}
				
		if(GlobalHammingDistances[step]<min_hd && (70<=step))
		{
			min_hd=GlobalHammingDistances[step];
			fl=step;
		}
		if(max_hd<GlobalHammingDistances[step] && (70<=step))
		{
			max_hd=GlobalHammingDistances[step];
		}
	}
	
	// get RSC
	double RSC=(double)(max_hd-GlobalHammingDistances[fl])/(max_hd-GlobalHammingDistances[rl]);
	chromosome[0].rsc=RSC;
*/

	// get fragment length FL and HD[FL] run through from (step_num-1),...,2*read_len+1
	int fl=0;
	int min_hd_fl=GlobalHammingDistances[step_num-1];
	int max_hd_fl=GlobalHammingDistances[step_num-1];
	for(int step=(step_num-1);(2*read_len+1)<=step;step--)
	{
		if(GlobalHammingDistances[step]<min_hd_fl)
		{
			min_hd_fl=GlobalHammingDistances[step];
			fl=step;
		}
		if(max_hd_fl<GlobalHammingDistances[step])
		{
			max_hd_fl=GlobalHammingDistances[step];
		}
	}

	// get phantom peak RL and HD[RL] run through from 1,...,2*read_len
	int rl=0;
	int min_hd_rl=GlobalHammingDistances[0];
	for(int step=1;step<=(2*read_len);step++)
	{
		if(GlobalHammingDistances[step]<min_hd_rl)
		{
			min_hd_rl=GlobalHammingDistances[step];
			rl=step;			
		}		
	}

	// get RSC
	double RSC=(double)(max_hd_fl-GlobalHammingDistances[fl])/(max_hd_fl-GlobalHammingDistances[rl]);
	chromosome[0].rsc=RSC;
	
	// write R-script that generates a plot
	std::string R_out = out_prefix + "-Q-binding-characteristics.R";
	std::string PDF_out = out_prefix + "-Q-binding-characteristics.pdf";
	std::ofstream OUT;
	OUT.open(R_out.c_str());
	
	OUT << "HD<-c(";
	for(int step=0;step<step_num;step++)
	{
		OUT << GlobalHammingDistances[step];
		if(step<step_num-1){OUT << ",";}
	}
	OUT << ")\n\n";
	
	OUT << "STRAND_SHIFT<-c(";
	for(int step=1;step<=step_num;step++)
	{
		OUT << step;
		if(step<step_num){OUT << ",";}
	}
	OUT << ")\n\n";
	
	OUT << "RL<-" << read_len << "\n";
	OUT << "RL<-STRAND_SHIFT[which(HD==min(HD[1:RL]))]" << "\n";
	OUT << "FL<-STRAND_SHIFT[which(HD==min(HD[(RL*2+1):(length(HD))]))]"      << "\n";
	OUT << "MAX_HD<-max(HD)"      << "\n";
	OUT << "HD_FL<-HD[which(STRAND_SHIFT==FL)]"      << "\n";
	OUT << "HD_RL<-HD[which(STRAND_SHIFT==RL)]"      << "\n";
		
	OUT << "RSC<-(MAX_HD-HD_FL)/(MAX_HD-HD_RL)"      << "\n";
	OUT << "RSC<-format(RSC,digits=3,nsmall=3)"      << "\n";
	
	OUT << "MAIN_PLOT<-paste(\"FL = \",FL,\"|\")"      << "\n";
	OUT << "MAIN_PLOT<-paste(MAIN_PLOT,\"RSC = \",RSC)"      << "\n";	
		
	OUT << "pdf(\"" << PDF_out << "\",height=7,width=7)" << "\n";
	OUT << "plot(STRAND_SHIFT,HD,xlab=\"strand shift\",ylab=\"hamming distance\",type=\"l\",main=MAIN_PLOT,cex.axis=1.3,cex.lab=1.5)"               << "\n";
	
	OUT << "abline(";
	OUT 	<< "v=FL" << ",";
	OUT 	<< "col=\"red\",lty=2";
	OUT << ")" << "\n";
	OUT << "abline(h=HD_FL,lty=2,col=\"red\")\n";
	OUT << "abline(v=RL,lty=2,col=\"grey\")\n";
	OUT << "abline(h=HD_RL,lty=2,col=\"grey\")\n";
	OUT << "abline(h=MAX_HD,lty=2)\n";
	OUT << "abline(v=0,lty=2)\n";	
	OUT << "dev.off()"              << "\n";
	
	OUT.close();

	return(fl+1);
}

int getQfragLengthDistribution(std::vector<Chromosome> &chromosome, int step_num, int thread_num, std::string out_prefix)
{
	// get number of chromosomes
	int chr_num=chromosome.size();
	
	// get read length
	int read_len=0;
	for(int i=0;i<chr_num;i++)
	{
		if(read_len<chromosome[i].read_len_chip)
		{
			read_len=chromosome[i].read_len_chip;
		}
	}
	
	//std::cout << "   number of chromosomes: " << chr_num << "\n";
	//std::cout << "   read length: " << read_len << "\n\n";

	int max=step_num;

	// init vector for distribution
	std::vector<int> QLD_CHIP(step_num+1,0);
	std::vector<int> QLD_CTRL(step_num+1,0);


	// get distribution of qfrag lengths for each chromosome
	
	for(int i=0;i<chr_num;i++)
	{
		int c_hit=0;
		int pos=0;
		
		// determine numbers of qfrags for treatment
		while(pos<chromosome[i].len)
		{
			// at position pos are one or more hits
			while(c_hit<chromosome[i].hit_num_chip && pos==chromosome[i].CHIP_HITS[c_hit].pos)
			{
				// current hit is on the forward strand
				if(chromosome[i].CHIP_HITS[c_hit].strand==0)
				{
					// look for qfrags
					int j=c_hit+1; // start examination with the next hit
				
					// second hit for qfrag is at most max bases apart
					while(j<chromosome[i].hit_num_chip && chromosome[i].CHIP_HITS[j].pos <= (chromosome[i].CHIP_HITS[c_hit].pos+max))
					{
						// second hit of the qfrag is at least min bases apart
						// and on the reverse strand
						if(chromosome[i].CHIP_HITS[j].strand==1)// && chromosome[i].CHIP_HITS[c_hit].pos+min<=chromosome[i].CHIP_HITS[j].pos)
						{
							QLD_CHIP[chromosome[i].CHIP_HITS[j].pos-chromosome[i].CHIP_HITS[c_hit].pos+1]++;
							//std::cout << "   qlen: " << chromosome[i].CHIP_HITS[j].pos-chromosome[i].CHIP_HITS[c_hit].pos+1 << "\n";
						}
						j++;
					}
				}
				c_hit++;
			}
			pos++;
		}

		// if no control, then pseudo control
		std::vector<Hit> CTRL_HITS;
		if(chromosome[i].hit_num_ctrl==0)
		{
			CTRL_HITS=getPseudoControlSwitchStrandAndFlip(chromosome[i].CHIP_HITS,read_len,chromosome[i].len);
		}
		else
		{
			CTRL_HITS=chromosome[i].CTRL_HITS;
		}
		
		c_hit=0;
		pos=0;
		
		// determine numbers of qfrags for control
		while(pos<chromosome[i].len)
		{
			// at position pos are one or more hits
			while(c_hit<chromosome[i].hit_num_chip && pos==CTRL_HITS[c_hit].pos)
			{//std::cout << CTRL_HITS[c_hit].pos << "*\n";
			
				// current hit is on the forward strand
				if(CTRL_HITS[c_hit].strand==0)
				{
					// look for qfrags
					int j=c_hit+1; // start examination with the next hit
				
					// second hit for qfrag is at most max bases apart
					while(j<chromosome[i].hit_num_chip && CTRL_HITS[j].pos <= (CTRL_HITS[c_hit].pos+max))
					{
						// second hit of the qfrag is at least min bases apart
						// and on the reverse strand
						if(CTRL_HITS[j].strand==1)// && chromosome[i].CHIP_HITS[c_hit].pos+min<=chromosome[i].CHIP_HITS[j].pos)
						{
							QLD_CTRL[CTRL_HITS[j].pos-CTRL_HITS[c_hit].pos+1]++;
						}
						j++;
					}
				}
				c_hit++;
			}
			pos++;
		}		
	}
	
	int estimated_len_1=0;
	int estimated_len_2=0; 
	for(int k=2;k<=read_len;k++)
	{
			if(QLD_CHIP[estimated_len_1]<QLD_CHIP[k])
			{
				estimated_len_1=k;
			}
			if(QLD_CHIP[estimated_len_2]-QLD_CTRL[estimated_len_2]<QLD_CHIP[k]-QLD_CTRL[k])
			{
				estimated_len_2=k;
			}
	}	
	//std::cout << "   Estimated fragment length (1): " << estimated_len_1 << "\n";
	//std::cout << "   Estimated fragment length (2): " << estimated_len_2 << "\n\n";
	
	// write R-script that generates a plot
	std::string R_out = out_prefix + "-Q-qfrag-binding-characteristics.R";
	std::string PDF_out = out_prefix + "-Q-qfrag-binding-characteristics.pdf";
	std::ofstream OUT;
	
	OUT.open(R_out.c_str());

	OUT << "QLD_CHIP<-c(";
	for(int step=2;step<step_num;step++)
	{
		OUT << QLD_CHIP[step];
		if(step<step_num-1){OUT << ",";}
	}
	OUT << ")\n\n";

	OUT << "QLD_CTRL<-c(";
	for(int step=2;step<step_num;step++)
	{
		OUT << QLD_CTRL[step];
		if(step<step_num-1){OUT << ",";}
	}
	OUT << ")\n\n";
	
	OUT << "Q_LENGTH<-c(";
	for(int step=2;step<step_num;step++)
	{
		OUT << step;
		if(step<step_num-1){OUT << ",";}
	}
	OUT << ")\n\n";
	
	OUT << "RL<-" << read_len << "\n";
	
	
	OUT << "MAIN<-" << "\"" <<  out_prefix << "\"" << "\n";
	
	OUT << "YMAX<-max(c(QLD_CHIP,QLD_CTRL))\n\n";
	OUT << "YMIN<-min(c(QLD_CHIP,QLD_CTRL))\n\n";


	OUT << "pdf(\"" << PDF_out << "\",height=3.4,width=6)" << "\n";
	OUT << "par(mfrow=c(1,2))" << "\n";

	OUT << "QL_MAX<-Q_LENGTH[which(max(QLD_CHIP)==QLD_CHIP)]" << "\n";
	
	OUT << "plot(Q_LENGTH,QLD_CTRL, type=\"l\", xlab=\"qfrag length\",ylab=\"# qfrags\",ylim=c(YMIN,YMAX),col=\"darkgrey\",main=MAIN)\n\n";
	OUT << "lines(Q_LENGTH,QLD_CHIP)\n\n";
	OUT << "abline(v=RL,lty=2,col=\"grey\")\n\n";
	OUT << "abline(v=QL_MAX,lty=2)\n\n";
	OUT << "legend(\"topright\",legend=c(paste(\"argmax(x)=\",QL_MAX,sep=\"\")),bg=\"white\",bty=\"n\",cex=0.7)" << "\n";
	

	OUT << "QL_MAX<-Q_LENGTH[which(max(QLD_CHIP-QLD_CTRL)==QLD_CHIP-QLD_CTRL)] #plot(Q_LENGTH,QLD_CHIP-QLD_CTRL, type=\"l\", xlab=\"qfrag length\",ylab=\"# qfrags\")\n\n";
	OUT << "QL_MIN<-Q_LENGTH[which(min(QLD_CHIP-QLD_CTRL)==QLD_CHIP-QLD_CTRL)] #plot(Q_LENGTH,QLD_CHIP-QLD_CTRL, type=\"l\", xlab=\"qfrag length\",ylab=\"# qfrags\")\n\n";
	
	OUT << "plot(Q_LENGTH,QLD_CHIP-QLD_CTRL, type=\"l\", xlab=\"qfrag length\",ylab=\"# chip_qfrags - # psc_qfrags\",main=MAIN,xlim=c(0,2*RL+10))\n\n";
	OUT << "abline(h=0)\n\n";
	OUT << "abline(v=RL,lty=2,col=\"grey\")\n\n";
	OUT << "abline(v=2*RL,lty=2,col=\"grey\")\n\n";
	OUT << "abline(v=QL_MAX,lty=2)\n\n";
	OUT << "abline(v=min(QLD_CHIP-QLD_CTRL),lty=2)\n\n";
	OUT << "legend(\"topright\",legend=c(paste(\"argmax(x)=\",QL_MAX,sep=\"\")),bg=\"white\",bty=\"n\",cex=0.7)" << "\n";

	OUT.close();
	return(estimated_len_2);
}

std::vector<Hit> getPseudoControlSwitchStrandAndFlip(std::vector<Hit> CHIP_HITS, int read_length, int chr_length)
{
	std::vector<Hit> CTRL_HITS=CHIP_HITS;
	
	for (unsigned int i = 0; i < CHIP_HITS.size(); i++)
	{
		if(CHIP_HITS[i].strand==0)
		{
			CTRL_HITS[i].strand=1;
			CTRL_HITS[i].pos=CHIP_HITS[i].pos+CHIP_HITS[i].read_length-1;
			if(chr_length<CTRL_HITS[i].pos)
			{
				CTRL_HITS[i].pos=chr_length;
			}

		}
		if(CHIP_HITS[i].strand==1)
		{
			CTRL_HITS[i].strand=0;
			CTRL_HITS[i].pos=CHIP_HITS[i].pos-CHIP_HITS[i].read_length+1;
			if(CTRL_HITS[i].pos<0)
			{
				CTRL_HITS[i].pos=0;
			}
		}
	}
	sort(CTRL_HITS.begin(),CTRL_HITS.end(),compareHitsByPos);
	return(CTRL_HITS);
} 

int cleanChipHits(Chromosome &chromosome)
{
	std::vector<Hit> CLEAN_HITS;
	
	int c_hit=0;
	int c_pos=0;
	int chr_f_hit_num=0;
	int chr_r_hit_num=0;
	while(c_pos<chromosome.len)
	{
		int f_hit_num=0;
		int r_hit_num=0;			
		while(c_hit<chromosome.hit_num_chip && c_pos==chromosome.CHIP_HITS[c_hit].pos)
		{
			if(chromosome.CHIP_HITS[c_hit].strand==0)
			{
								//std::cout  << "*\n";

				f_hit_num++;
			}
			if(chromosome.CHIP_HITS[c_hit].strand==1)
			{
				r_hit_num++;
			}
			c_hit++;
		}
		

		
		if((0==f_hit_num && 1<=r_hit_num) || (0==r_hit_num && 1<=f_hit_num))
		{
			std::cout << chromosome.name << "\t";
			std::cout << c_pos << "\t";
			std::cout << "fhn:" << f_hit_num << "\t";
			std::cout << "rhn:" << r_hit_num << "\n";

			int hit_idx=c_hit-1;
			//push all hits from current position to clean hit set
			while(0<hit_idx && c_pos==chromosome.CHIP_HITS[hit_idx].pos)
			{
				// push
				CLEAN_HITS.push_back(chromosome.CHIP_HITS[hit_idx]);
				if(chromosome.CHIP_HITS[hit_idx].strand==0)
				{
					chr_f_hit_num++;
				}
				else
				{
					chr_r_hit_num++;
				}
				hit_idx--;
			}
		}	
		c_pos++;
	}
	chromosome.CHIP_HITS=CLEAN_HITS;
	chromosome.hit_num_chip=chr_f_hit_num+chr_r_hit_num;
	return(0);
}
