#include "Q_ReadInFiles.h"
#include <boost/unordered_map.hpp>
#include "Q_GetFragLength.h"


bool compareHitsByPos(const Hit& a, const Hit& b)
{
    return a.pos < b.pos;
}


/**
 * \brief Function which reads the input alignment into an array of 
 * \ref Chromosome structs.
 * 
 * @param[in] &chromosome A vector of chromosomes
 * @param[in] &chr_num Total number of chromosomes
 * @param[in] input_file Path to the input file
 * 
 * This fuction uses the 
 * <a href="http://trac.seqan.de/wiki/Tutorial/BasicSamBamIO">sam/bam reader of seqan</a>.
 * 
 * This function takes an empty array of \ref Chromosome objects,
 * an empty variable \ref chr_num (declared in \ref Q.cpp) for the number
 * of chromosomes and an input file in sam or bam format.
 *  
 * Using the information in the header of the sam/bam file a 
 * \ref Chromosome object for each chromosome sequence is created and
 * added to an array of \ref Chromosome objects (declared in \ref Q.cpp).
 * Furthermore the name and the length of the chromosome are assigned 
 * to the chromosome object. In addition the variable chr_num is set to 
 * the number of chromosomes found in the sam/bam header. 
 * 
 * After the \ref Chromosome objects were created the function iterates 
 * over the sam/bam file line by line. If a line has a bitwise FLAG of 
 * 0 or 16 a \ref Hit object is created and added to the 
 * \ref Chromosome.HITS corresponding chromosome. That is, at first the
 * hits are sorted according to chromosome.
 * 
 * After iterating over the whole sam/bam file the \ref Hit vectors of
 * the Chromosome objects are sorted according position. This is done 
 * in parallel using SEQAN_OMP_PRAGMA using a number of threads that
 * is defined in \ref Options.thread_num.
 * 
 * 
 */
int ReadAlignmentFile(std::vector<Chromosome> &chromosome, int &chr_num, seqan::CharString chip_sample, seqan::CharString control_sample, bool keep_dup, int thread_num, bool use_pseudo_control)
{
	// Open input stream, BamStream can read SAM and BAM files
    seqan::BamStream bamStreamInChIP(toCString(chip_sample));
    if(!isGood(bamStreamInChIP))
    {
        std::cerr << "ERROR: Could not open " << chip_sample << " !\n";
        return 1;
    }

	// create chromosome objects
	chr_num=length(bamStreamInChIP.header.sequenceInfos);

	// maps chip-chromosome name to chip-rID
	boost::unordered_map <std::string, int> chr_map_aux;
	// maps control-chromosome rID to chip-chromosome rID
	boost::unordered_map <int, int> chr_map;

	for(int i=0;i<chr_num;i++)
	{
		Chromosome c;
		c.name=bamStreamInChIP.header.sequenceInfos[i].i1;
		c.len=bamStreamInChIP.header.sequenceInfos[i].i2;
		c.hit_num_chip=0;
		chromosome.push_back(c);
		chr_map_aux[toCString(c.name)]=i;
	}

	// read chip hits to chromosome objects
	seqan::BamAlignmentRecord record;
	while(!atEnd(bamStreamInChIP))
	{
		// get the read length from the first alignment record
		if(chromosome[record.rID].read_len_chip==-1)
		{
			chromosome[record.rID].read_len_chip=length(record.seq);
		}
	
	
		if (readRecord(record, bamStreamInChIP) != 0)
        {
            std::cerr << "ERROR: Could not read record!\n";
            return 1;
        }

		Hit hit;
		if(record.flag==16||record.flag==147||record.flag==83||record.flag==121||record.flag==153||record.flag==185||record.flag==115||record.flag==179||record.flag==81||record.flag==145||record.flag==113||record.flag==177)
		{
			hit.pos=record.beginPos+length(record.seq)-1;
			hit.strand=1;
			hit.read_length=length(record.seq);
		    chromosome[record.rID].CHIP_HITS.push_back(hit);
		    chromosome[record.rID].hit_num_chip++;
		    chromosome[record.rID].r_hit_num_chip++;
		}
		if(record.flag==0||record.flag==99||record.flag==163||record.flag==73||record.flag==89||record.flag==137||record.flag==131||record.flag==97||record.flag==161||record.flag==67||record.flag==65||record.flag==129)
		{
			hit.pos=record.beginPos;
			hit.strand=0;
			hit.read_length=length(record.seq);
			chromosome[record.rID].CHIP_HITS.push_back(hit);
		    chromosome[record.rID].hit_num_chip++;
		    chromosome[record.rID].f_hit_num_chip++;		    
		}
		
		if(chromosome[record.rID].read_len_chip<(int)length(record.seq))
		{
			chromosome[record.rID].read_len_chip=length(record.seq);
		}
	}


	// If there is a control,
    if(control_sample != "None")
    {
		// open an input stream for the control
		seqan::BamStream bamStreamInControl(toCString(control_sample));
		if(!isGood(bamStreamInControl))
		{
			std::cerr << "ERROR: Could not open " << control_sample << " !\n";
			return 1;
		}
		
		// iterate over control chromosomes
		int chr_num_ctrl=length(bamStreamInControl.header.sequenceInfos);
		for(int i=0;i<chr_num_ctrl;i++)
		{
			seqan::CharString name = bamStreamInControl.header.sequenceInfos[i].i1;
			// the chromosome 'name' has also hits for chip
			if(chr_map_aux.find(toCString(name)) != chr_map_aux.end())
			{
				chr_map[i] = chr_map_aux[toCString(name)];
			}
			else
			// the chromosome 'name' has no hits for chip
			{
				chr_map[i] = -1;
			}
		}
		
		seqan::BamAlignmentRecord record;
		while(!atEnd(bamStreamInControl))
		{
			if (readRecord(record, bamStreamInControl) != 0)
			{
				std::cerr << "ERROR: Could not read record!\n";
				return 1;
			}
			
			// if for the corresponding chromosome there is no hit for the ChIP sample
			if(chr_map[record.rID] == -1){continue;};
			//std::cout << record.rID << "\t" << chr_map[record.rID] << "\n";

			Hit hit;
			if(record.flag==16||record.flag==147||record.flag==83||record.flag==121||record.flag==153||record.flag==185||record.flag==115||record.flag==179||record.flag==81||record.flag==145||record.flag==113||record.flag==177)
			{
				hit.pos=record.beginPos+length(record.seq)-1;
				hit.strand=1;
				hit.read_length=length(record.seq);
				chromosome[chr_map[record.rID]].CTRL_HITS.push_back(hit);
				chromosome[chr_map[record.rID]].hit_num_ctrl++;
				chromosome[chr_map[record.rID]].r_hit_num_ctrl++;
			}
			if(record.flag==0||record.flag==99||record.flag==163||record.flag==73||record.flag==89||record.flag==137||record.flag==131||record.flag==97||record.flag==161||record.flag==67||record.flag==65||record.flag==129)
			{
				hit.pos=record.beginPos;
				hit.strand=0;
				hit.read_length=length(record.seq);
				chromosome[chr_map[record.rID]].CTRL_HITS.push_back(hit);
				chromosome[chr_map[record.rID]].hit_num_ctrl++;
				chromosome[chr_map[record.rID]].f_hit_num_ctrl++;		    
			}
		}
    }
  
	// pseudo control
	if(use_pseudo_control)
	{
		for(int i=0;i<chr_num;i++)
		{
			chromosome[i].CTRL_HITS=getPseudoControlSwitchStrandAndFlip(chromosome[i].CHIP_HITS,chromosome[i].read_len_chip,chromosome[i].len);
			chromosome[i].hit_num_ctrl=chromosome[i].hit_num_chip;
			chromosome[i].f_hit_num_ctrl=chromosome[i].f_hit_num_chip;
			chromosome[i].r_hit_num_ctrl=chromosome[i].r_hit_num_chip;
		}
	}

	
	// sort hits according to starting position
	SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
	for(int i=0;i<chr_num;i++)
	{
		sort(chromosome[i].CHIP_HITS.begin(),chromosome[i].CHIP_HITS.end(),compareHitsByPos);
		sort(chromosome[i].CTRL_HITS.begin(),chromosome[i].CTRL_HITS.end(),compareHitsByPos);
	}

	if(keep_dup)
	{
		SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
		for(int i=0;i<chr_num;i++)
		{
			chromosome[i].redundant_hits_num_chip=RemoveRedundantHits(chromosome[i].CHIP_HITS,true);
			
			if(control_sample != "None")
			{
				if(chromosome[i].CTRL_HITS.size()==0){continue;}
				chromosome[i].redundant_hits_num_ctrl=RemoveRedundantHits(chromosome[i].CTRL_HITS,true);
			}
		}
		return 0;
	}
	
	// if not keep-dup: remove redundant hits
	// --------------------------------------
	
	SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
	for(int i=0;i<chr_num;i++)
	{
		chromosome[i].redundant_hits_removed_num_chip=RemoveRedundantHits(chromosome[i].CHIP_HITS,false);
		chromosome[i].redundant_hits_num_chip=chromosome[i].redundant_hits_removed_num_chip;
		
		// re-determine f_hit_num_chip
		chromosome[i].f_hit_num_chip=0;
		chromosome[i].r_hit_num_chip=0;
		for(unsigned int j=0;j<chromosome[i].CHIP_HITS.size();j++)
		{
			if(chromosome[i].CHIP_HITS[j].strand)
			{
				chromosome[i].f_hit_num_chip++;
			}
			else
			{
				chromosome[i].r_hit_num_chip++;
			}
		}
		chromosome[i].hit_num_chip=chromosome[i].f_hit_num_chip+chromosome[i].r_hit_num_chip;
	}

    if(control_sample != "None")
    {
		SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
		for(int i=0;i<chr_num;i++)
		{
			// if there are hits for chip but not for control
			if(chromosome[i].CTRL_HITS.size()==0){continue;}

			chromosome[i].redundant_hits_removed_num_ctrl=RemoveRedundantHits(chromosome[i].CTRL_HITS,false);
			chromosome[i].redundant_hits_num_ctrl=chromosome[i].redundant_hits_removed_num_ctrl;
			
			// re-determine f_hit_num_chip
			chromosome[i].f_hit_num_ctrl=0;
			chromosome[i].r_hit_num_ctrl=0;
			for(unsigned int j=0;j<chromosome[i].CTRL_HITS.size();j++)
			{
				if(chromosome[i].CTRL_HITS[j].strand)
				{
					chromosome[i].f_hit_num_ctrl++;
				}
				else
				{
					chromosome[i].r_hit_num_ctrl++;
				}
			}
			chromosome[i].hit_num_ctrl=chromosome[i].f_hit_num_ctrl+chromosome[i].r_hit_num_ctrl;
		}	
	}
	
    return 0;
}

/**
 * \brief Function which removes redundant hits.
 * 
 * Function which removes redundant hits from a vector of hits sorted by position.
 * 
 * If keep is true the hits will not be removed,
 * but the number of redundant hits will be returned.
 * 
 */
int RemoveRedundantHits(std::vector<Hit> &HITS, bool keep)
{
	std::vector<Hit> HITS_RMDUP;
	
	// if the chromosome object has no hits
	if(HITS.size()==0)
	{
		HITS=HITS_RMDUP;
		return 0;
	}
	
	int seen_strands;
	int cur_pos=HITS[0].pos;
	HITS_RMDUP.push_back(HITS[0]);
	if(HITS[0].strand==true){seen_strands=0;}else{seen_strands=1;}
	
	for(unsigned int j=1;j<HITS.size();j++)
	{
		if(cur_pos!=HITS[j].pos)
		{
			cur_pos=HITS[j].pos;
			if(HITS[j].strand==true){seen_strands=0;}else{seen_strands=1;}
			HITS_RMDUP.push_back(HITS[j]);
		}
		else
		{
			// only f strand has been seen AND hit is on r strand
			if(seen_strands==0 && !HITS[j].strand)
			{
				HITS_RMDUP.push_back(HITS[j]);
				seen_strands=2;
			}
			// only r strand has been seen AND hit is on f strand
			else if(seen_strands==1 && HITS[j].strand)
			{
				HITS_RMDUP.push_back(HITS[j]);
				seen_strands=2;
			}		
		}	
	}
	int redundant_hits=HITS.size()-HITS_RMDUP.size();
	if(!keep)
	{
		HITS=HITS_RMDUP;
	}
	return redundant_hits;
}


/**
 * \brief Function which reads summits from a BED file into the summit 
 * vectors of the individual chromosomes.
 * 
 * This function takes an array of \ref Chromosome objects,
 * which are already filled with hits and a path to a BED file.
 * 
 * Instead of detecting the summits itself, the summits will be read
 * from the BED file into the summit vectors of the individual
 * chromosomes. The summits will be sorted by position and an
 * anchor hit will be determined for each summit for ChIP and control
 * sample.
 * 
 * For a given summit only position and anchor hit will be set,
 * all other attributes will have initial values.
 * 
 * Chromosome IDs in BED file must be consistent with IDs in SAM/BAM file.
 * 
 * Only the first three columns of the BED file will be considered
 * and summit position is always the center of a given region.
 * 
 */
int ReadBedFile(std::vector<Chromosome> &chromosome, seqan::CharString bed_hit_dist)
{
	std::map<std::string, std::vector<int> > Summits;
	std::ifstream fin;
	char* file_name = toCString(bed_hit_dist);
	fin.open(file_name);
	std::string line,name;
	int sta, end;
	while(!fin.eof())
	{
		getline(fin, line);
		if(line=="") break;
		std::istringstream tokens(line);
		tokens >> name;
		tokens >> sta;
		tokens >> end;
		Summits[name].push_back(sta+((end-sta)/2));
	}
	fin.close();
	
	// sort summit vectors by position
	typedef std::map<std::string, std::vector<int> >::iterator it_type;
	for(it_type iterator = Summits.begin(); iterator != Summits.end(); iterator++)
	{
		sort(Summits[iterator->first].begin(),Summits[iterator->first].end());
	}
	
	// Add summits to chromosomes
	for(unsigned int i=0;i<chromosome.size();i++)
	{
		std::vector<Summit> NewSummits;
		for(unsigned int j=0;j<Summits[toCString(chromosome[i].name)].size();j++)
		{
			Summit s;
			s.pos=Summits[toCString(chromosome[i].name)][j];
			NewSummits.push_back(s);
			chromosome[i].sum_num++;
		}
		chromosome[i].SUMMITS=NewSummits;
	}
	
	
	
	// get anchor hits for summits
	for(unsigned int i=0;i<chromosome.size();i++)
	{
		int c_hit=0; // ChIP sample
		for(int j=0;j<chromosome[i].sum_num;j++)
		{
			while(chromosome[i].CHIP_HITS[c_hit].pos<chromosome[i].SUMMITS[j].pos)
			{
				c_hit++;
			}
			chromosome[i].SUMMITS[j].anchor_hit_chip=c_hit;
		}
		
		if(chromosome[i].hit_num_ctrl==0){continue;}
		
		c_hit=0; // Control sample
		for(int j=0;j<chromosome[i].sum_num;j++)
		{
			while(chromosome[i].CTRL_HITS[c_hit].pos<chromosome[i].SUMMITS[j].pos)
			{
				c_hit++;
			}
			chromosome[i].SUMMITS[j].anchor_hit_ctrl=c_hit;
		}
	}
	return 0;	
}
