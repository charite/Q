/**
 * @file
 * @authors  Peter Hansen <peter.hansen@charite.de>
 * @version 0.01
 *
 * \brief Read in the alignment data
 *  
 * For each chromosome i.e. reference sequence a \ref Chromosome struct
 * is created. A \ref Chromosome object contains information about the 
 * corresponding chromosome as name, length or number of hits.
 * 
 * Furthermore a \ref Chromosome structure contains a vector of Hits \ref Chromosome.HITS,
 * where each \ref Hit has a position, which corresponds to the fragment
 * end that was sequenced, and a strand variable, which is true for
 * the forward strand and false for the reverse strand. In addition
 * a chromosome structure contains a vector of Summits \ref Chromosome.SUMMITS, where a summit 
 * is a detected signal.
 * 
 * \section future_plans Future plans
 * 
 * Actually this should be the class Chromosome and the single function
 * \ref ReadAlignmentFile should be a constructor for this class.
 * There could also be a class Feature which is a more general 
 * description of hits and summits. The array of chromosomes could be
 * replaced by a class Genome. This would be a reusable structure.  
 * 
 */

#ifndef READ_IN_FILES_H
#define READ_IN_FILES_H

#include <seqan/bam_io.h>


/**
 * \brief Struct to store a hit i.e. position and strand.
 * 
 * A hit corresponds to a line in a input alignment file. A hit has a 
 * position which corresponds to the 5' end of a fragment and a strand.
 * Position and strand are determined by means of the <a href="http://samtools.sourceforge.net/SAM1.pdf">SAM format 
 * specification</a>.
 * 
 * A hit is on the forward strand, if the bitwise FLAG in column two
 * of the sam file equals 0. A hit is on the reverse strand,
 * if the FLAG equals 16. Lines with a FLAG different from 0 or 16
 * are skipped.
 * 
 * The position POS in column four is the 1-based leftmost position
 * of the mapping position. For forward strand hits this corresponds
 * to the 5' end of the read and the position specified in the sam file
 * is also the position of the hit. For reverse strand hits POS 
 * corresponds to the 3' end of the read and the position of the hit is
 * POS + READ_LENGTH -1.
 * 
 */
struct Hit
{
	/** \brief The position of a hit corresponds to an 5' end of a fragment. */
	int pos;
	
	/** \brief A hit is either on the forward strand (T) or the reverse strand (F). */	
	bool strand;
	
	/** \brief True if hit belongs to a qfrag. */	
	bool is_q_hit;

	Hit() :
		is_q_hit(false)
	{}	
};


/**
 * \brief Auxiliary function for sorting hits.
 */
bool compareHitsByPos(const Hit& a, const Hit& b);


/**
 * \brief Struct to store information about a single detected position.
 */
struct Summit
{
	/** \brief Position of the summit .
	 * Is calculated in \ref getSummits */
	int pos;
	
	/** \brief Index of the first hit after the summit.
	* Is calculated in \ref getSummits */	
	int anchor_hit_chip;
	
	/** \brief Index of the first hit after the summit.
	* Is calculated in \ref getSummits */	
	int anchor_hit_ctrl;

	/** \brief Number of qfrags that cover the summit position.
	 * Is calculated in \ref getSummits */	
	int q_cov;
	
	/** \brief Number of qfrags that cover the summit position for the
	 * ChIP sample. Is calculated in \ref getQFragsCoverage */	
	int q_cov_chip;
	
	/** \brief Number of qfrags that cover the summit position for the
	 * control sample. Is calculated in \ref getQFragsCoverage */	
	int q_cov_ctrl;
	
	/** \brief Number of positions within a given interval (ChIP sample) before the 
	 * summit that are covered by at least one hit on the forward strand
	 * that is part of a qfrag that covers the summit. 
	 * Is calculated by the function \ref getQfragEnds*/
	int q_beg_chip;
	
	/** \brief Number of positions within a given interval (control sample) before the 
	 * summit that are covered by at least one hit on the forward strand
	 * that is part of a qfrag that covers the summit. 
	 * Is calculated by the function \ref getQfragEnds*/
	int q_beg_ctrl;
	
	/** \brief Number of positions within a given interval (ChIP sample) behind the 
	 * summit that are covered by at least one hit on the reverse strand
	 * that is part of a qfrag that covers the summit. 
	 * Is calculated by the function \ref getQfragEnds*/
	int q_end_chip;
	
	/** \brief Number of positions within a given interval (control sample) behind the 
	 * summit that are covered by at least one hit on the reverse strand
	 * that is part of a qfrag that covers the summit. 
	 * Is calculated by the function \ref getQfragEnds*/
	int q_end_ctrl;	

	/** \brief Number of forward hits summit upstream within distance of radius 
	 * Is calculated by the function \ref getKs (ChIP sample)*/
	int kfu_chip;
	
	/** \brief Number of forward hits summit downstream within distance of radius
	 * Is calculated by the function \ref getKs (ChIP sample)*/		
	int kfd_chip;
	
	/** \brief Number of reverse hits summit upstream within distance of radius
	 * Is calculated by the function \ref getKs (ChIP sample)*/
	int kru_chip;

	/** \brief Number of reverse hits summit downstream within distance of radius
	 * Is calculated by the function \ref getKs (ChIP sample)*/	
	int krd_chip;
	
	/** \brief Number of forward hits summit upstream within distance of radius 
	 * Is calculated by the function \ref getKs (control sample)*/
	int kfu_ctrl;
	
	/** \brief Number of forward hits summit downstream within distance of radius
	 * Is calculated by the function \ref getKs (control sample)*/		
	int kfd_ctrl;
	
	/** \brief Number of reverse hits summit upstream within distance of radius
	 * Is calculated by the function \ref getKs (control sample)*/
	int kru_ctrl;

	/** \brief Number of reverse hits summit downstream within distance of radius
	 * Is calculated by the function \ref getKs (control sample)*/	
	int krd_ctrl;
	
	long double saturation_score;	

	long double p_value;
	
	long double q_value;
	
	Summit() :
		pos(0), anchor_hit_chip(-1), anchor_hit_ctrl(-1), q_cov_chip(0), q_cov_ctrl(0), q_beg_chip(0), q_beg_ctrl(0), q_end_chip(0), q_end_ctrl(0), kfu_chip(0), kfd_chip(0), kru_chip(0), krd_chip(0), kfu_ctrl(0), kfd_ctrl(0), kru_ctrl(0), krd_ctrl(0), saturation_score(-1), p_value(-1), q_value(-1)
	{}
};



// create new struct Summit_chr
struct Summit_chr{
	Summit summit;
	seqan::CharString chr_name;
	int chr_length;
};


/**
 * \brief Struct to store information associated with a single chromosome.
 */
struct Chromosome
{
	/** \brief The name of the chromosome e.g. chr1 */
	seqan::CharString name;
	
	/** \brief The length of the chromosome */
	int len;

	/** \brief Length of the first read in treatment input data */
	int read_len_chip;

	/** \brief The total number of hits on the chromosome for ChIP sample*/	
	int hit_num_chip;

	/** \brief The total number of hits on the chromosome for control sample */	
	int hit_num_ctrl;

	/** \brief The total number of redundant ChIP hits that were removed*/	
	int redundant_hits_removed_num_chip;

	/** \brief The total number of redundant control hits that were removed*/	
	int redundant_hits_removed_num_ctrl;	

	/** \brief The total number of redundant ChIP hits that were removed*/	
	int redundant_hits_num_chip;

	/** \brief The total number of redundant control hits that were removed*/	
	int redundant_hits_num_ctrl;
	
	/** \brief The number of forward strand hits on the chromosome ChIP sample */
	int f_hit_num_chip;
	
	/** \brief The number of forward strand hits on the chromosome control sample */
	int f_hit_num_ctrl;
	
	/** \brief The number of reverse strand hits on the chromosome ChIP sample */	
	int r_hit_num_chip;
	
	/** \brief The number of reverse strand hits on the chromosome control sample */	
	int r_hit_num_ctrl;	
	
	/** \brief A vector containing all ChIP hits on the chromosome */		
	std::vector<Hit> CHIP_HITS;
	
	/** \brief A vector containing all control hits on the chromosome */		
	std::vector<Hit> CTRL_HITS;

	/** \brief The total observed number of qfrags on the chromosome for ChIP sample */	
	int qfrag_num_chip;
	
	/** \brief A vector containing all summits on the chromosome */			
	std::vector<Summit> SUMMITS;
	
	/** \brief A map containing the hamming distances for all shift sizes */			
	std::map<int, int> HAMMING_DISTANCES;

	/** \brief The total number of summits on the chromosome */	
	int sum_num;
	
	/** \brief Relative strand correlation RSC */	
	double rsc;
	
	
	Chromosome() :
		read_len_chip(-1), hit_num_chip(0), hit_num_ctrl(0), redundant_hits_removed_num_chip(0), redundant_hits_removed_num_ctrl(0), f_hit_num_chip(0), f_hit_num_ctrl(0), r_hit_num_chip(0), r_hit_num_ctrl(0), qfrag_num_chip(0), sum_num(0), rsc(0)
	{}
};

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
int ReadAlignmentFile(std::vector<Chromosome> &chromosome, int &chr_num, seqan::CharString chip_sample, seqan::CharString control_sample, bool keep_dup, int thread_num);


/**
 * \brief Function which removes redundant hits.
 * 
 * Function which removes redundant hits from a vector of hits sorted by position.
 * 
 */
int RemoveRedundantHits(std::vector<Hit> &HITS, bool keep);


/**
 * \brief Function which reads summits from a BED file. 
 * 
 * This function takes an array of \ref Chromosome objects,
 * filled with hits and an input file in bed format.
 * Only the first three columns of the BED file are used.
 * 
 */
int ReadBedFile(std::vector<Chromosome> &chromosome, seqan::CharString bed_hit_dist);

#endif
