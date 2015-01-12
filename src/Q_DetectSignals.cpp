#include "Q_ReadInFiles.h"
#include "Q_DetectSignals.h"

/**
 * \brief Function which detects all summits.
 * 
 * @param[in] &chromosome A \ref Chromosome object
 * @param[in] min Minimal distance between two hits forming a qfrag.
 * @param[in] max Maximal distance between two hits forming a qfrag.
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * This function proceeds in four steps:
 *	- 1. Find summits in qfrags coverage profile.
 *	- 2. Refine summits
 *	- 3. Get hit indices for each summut as anchor for evaluation.
 * 
 * \section step1 Detection
 * In a first step the summits are detected. The qfrag coverage for a 
 * given position in the chromosome is the number of qfrags that cover
 * this position. A summit is a set of consecutive position in the 
 * chromosome with the same qfrags coverage, where the qfrags coverage
 * at the position before and the position behind the summit is lower.
 * The summit detection is carried out locally, that is the qfrag
 * coverage profile for the whole chromosome does not have to be stored
 * in main memory. Due to this fact hromosomes can be processed in parallel. 
 *  
 * \section step2 Refinement
 * A \ref Summit which has a \ref Summit with a higher qfrag coverage in 
 * a distance of mindist=max are discarded. The set of refinded summits 
 * is assigned to the chromosome.
 *  
 * \section step3 Anchorage
 * For each summit the index of the first hit behind the summit is 
 * determined and stored as a property of the summit.
 * 
 */
int getSummits(Chromosome &chromosome, int min, int max, int cutoff)
{

	int qfrag_num=0;               // number of qfrags -> store this in chromosome.qfrag_num
	int c_cov=0;                   // current coverage
	int c_hit=0;                   // index of current hit
	const int circle_size=1000;    // size of the vector used for local qfrags coverage
	
	int sum_num=0;                 // number of summits -> store this in chromosome.sum_num
	int sta=0;
	int state=0;
	
	if(chromosome.hit_num_chip==0)
	{
		chromosome.qfrag_num_chip=0;
		chromosome.sum_num=0;
		return 0;
	}
	
	// coverage profile for qfrags
	// ---------------------------
	
	std::vector<int> Circle(circle_size); // filled with 0 by default
	std::vector<Summit> RawSummits;

	// loop in circle while end of the chromosome is reached
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
						c_cov++;
						Circle[(chromosome.CHIP_HITS[j].pos+1)%circle_size]--;
						qfrag_num++;
						chromosome.CHIP_HITS[c_hit].is_q_hit=true;						
						chromosome.CHIP_HITS[j].is_q_hit=true;
					}
					j++;
				}
			}
			c_hit++;			
		}
		c_cov=Circle[pos%circle_size]+c_cov;
		Circle[pos%circle_size]=c_cov;
		
		chromosome.qfrag_num_chip=qfrag_num;

		
		// get summits
		// -----------
	
		if(1<pos)
		{
			// coverage increases for pos-2 to pos-1
			if(Circle[(pos-2)%circle_size]<Circle[(pos-1)%circle_size])
			{
				sta=pos-1;
				state=1;
			}
			// coverage decreases for pos-1 to pos
			if((state==1) && (Circle[(pos-1)%circle_size]>Circle[(pos)%circle_size])&& (cutoff<Circle[(pos-1)%circle_size]))// 
			{
				sum_num++;
				Summit s;
				s.pos=sta+((pos-sta)/2);
				s.q_cov=Circle[((pos-1)%circle_size)];
				RawSummits.push_back(s);
				state=0;
			}
			Circle[((pos-2)%circle_size)]=0;
		}
		pos++;
	}



	// refine summits
	// --------------
	
	std::vector<Summit> RefinedSummits;
	
	int mindist=max;
	
	for(int i=0;i<sum_num;i++)
	{
		int j=0;
		int kickout=1;
		// check summits to the left
		while((i-j>0) && ((RawSummits[i].pos-RawSummits[i-j].pos)<mindist))
		{
			if(RawSummits[i-j].q_cov>RawSummits[i].q_cov)
			{
				kickout=0;
				break;
			}
			j++;
		}
		// check summits to the right
		j=0;
		while((i+j<sum_num) && ((RawSummits[i+j].pos-RawSummits[i].pos)<mindist))
		{
			if(RawSummits[i+j].q_cov>RawSummits[i].q_cov)
			{
				kickout=0;
				break;
			}
			j++;
		}
		if(kickout)
		{
			Summit s;
			s.pos=RawSummits[i].pos;
			s.q_cov=RawSummits[i].q_cov;
			RefinedSummits.push_back(s);
		}

	}

	// merge summits that are less then mindist apart

	Summit current_summit;
	if(RefinedSummits.size()>0){current_summit=RefinedSummits[0];}
	sum_num = RefinedSummits.size();
	int i=0;
	int shift=0;
	while(i < sum_num-1)
	{
		// determine the distance between current summit and summit downstream
		int dist = RefinedSummits[i+1].pos - current_summit.pos;
		
		// If the distance is smaller than mindist,
		// shift current summit by dist/2 in downstream direction
		if(dist<mindist)
		{
			current_summit.pos = current_summit.pos + dist/2;
			shift++;
		}	
		// If the distance is grater or equal than mindist,
		// take current summit as final summit and proceed with the next summit
		else
		{
			shift=0;
			chromosome.SUMMITS.push_back(current_summit);
			current_summit=RefinedSummits[i+1];
		}
		i++;
	}
	
	chromosome.sum_num = chromosome.SUMMITS.size();


	
	// get anchors (first hit after summit)
	// ------------------------------------
	
	c_hit=0;
	for(int i=0;i<chromosome.sum_num;i++)
	{
		while(chromosome.CHIP_HITS[c_hit].pos<chromosome.SUMMITS[i].pos)
		{
			c_hit++;
		}
		chromosome.SUMMITS[i].anchor_hit_chip=c_hit;
	}
	
	
	if(chromosome.hit_num_ctrl==0){return 0;}
	
	// if there is a control, get the anchor hits for control
	c_hit=0;
	for(int i=0;i<chromosome.sum_num;i++)
	{
		while(chromosome.CTRL_HITS[c_hit].pos<chromosome.SUMMITS[i].pos)
		{
			c_hit++;
		}
		chromosome.SUMMITS[i].anchor_hit_ctrl=c_hit;
	}
	
	return 0;
}
