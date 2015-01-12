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
