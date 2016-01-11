 /**
 * @file
 * @version 0.02
 *
 * \brief Parse the command line and store the results into
 * a \ref Options struct.
 * 
 */

#ifndef Q_PARSE_COMMAND_LINE_H
#define Q_PARSE_COMMAND_LINE_H

#include <seqan/arg_parse.h>

/**
 * \brief Struct to store the command line options and arguments.
 */
struct Options
{
	/** File name of ChIP-seq sample (sam or bam) */
	seqan::CharString chip_sample;
	
	/** File name of ChIP-seq control sample (sam or bam) */
	seqan::CharString control_sample;

	/** If set, pseudo control will be generated from the treatment data */     
	bool use_pseudo_control;
	
	/** If TRUE, duplicate reads will be kept. */     
	bool keep_dup;
	
	/** Prefix for all output files that will be generated */
	seqan::CharString out_prefix;
    
    /** Name of the output file that contains general information
     * about the execution of qfrags */
	seqan::CharString out_info_file;

    /** Name of the output file that contains general information
     * about the execution of qfrags */
	seqan::CharString binding_characteristics_file;
	
    /** Name of the output file that contains information
     * about chromosomes */
	seqan::CharString out_chromosome_info_file;

    /** Name of the output file that contains detailed information
     * for each summit */
	seqan::CharString out_summit_info_file;

    /** Name of the output file that contains all summits in ENCODE
     * narrowPeak format */ 
	seqan::CharString out_narrowPeak_file;
	
	/** Number of steps for the determination of the average fragment length */ 
	unsigned step_num;
	
	/** If TRUE only the binding characteristics will be determined and peak calling will be skipped */
	bool binding_characteristics_only;

	/** If TRUE additional output will be printed to the screen */     
	bool make_qfrag_length_distribution;
	
	/** Average fragment length */ 
	int fragment_length_avg;

	/** Deviation of fragment lengths */ 
	unsigned fragment_length_dev;

	/** Minimal length of a qfrag */ 
	unsigned lowerbound;
	    
    /** Maximal length of a qfrag */     
	unsigned upperbound;
	
	/** Maximal number of peaks to be written to file */
	int top_n;
	
	/** Cutoff for adjusted P-value */
	double p_value_cutoff;

	/** If TRUE, duplicate reads will be kept. */     
	bool nexus_mode;
   
    /** Maximal number of chromosomes to be processed in parallel */         
	unsigned thread_num;
	
	/** If TRUE additional output will be printed to the screen */     
	bool verbose;
	
	/** Chromosome ID consistent with ID in SAM/BAM file.
	 * If set, BED files for reads, shifted reads, extended reads
	 * qfrags will be written for the given chromosome ID.
	 * Peak calling will be skipped. */
	seqan::CharString write_bed;
	
	/** For the summits of the given BED file,
	 * the distribution of hits around summits on forward and reverse
	 * strand will be written to a text file.
	 * 
	 * Summit position is always the center of a given region.
	 * Chromosome IDs in BED file must be consistent with 
	 * IDs in SAM/BAM file. Default radius around summits is 1000.
	 * The radius can be changed via the '-r' option.
	 * 
	 * Output is a tab separated table containing three columns and
	 * 2 times radius rows. The first column contains the relative
	 * positions to the summits. The second and third column
	 * contain the accumulated counts of hits for all summits
	 * in the BED file. Second column for forward and third column for
	 * reverse strand. 
	 */
	seqan::CharString bed_hit_dist;
	unsigned bed_radius;
	
	bool write_bedgraph_treatment;
	seqan::CharString out_bedgraph_treatment_file;
	bool write_bedgraph_control;
	seqan::CharString out_bedgraph_control_file;	
	
	Options() :
		control_sample("None"), keep_dup(true), step_num(1000), binding_characteristics_only(false), make_qfrag_length_distribution(false), fragment_length_avg(-1), fragment_length_dev(50), p_value_cutoff(-1.0), nexus_mode(false), thread_num(1), verbose(false), write_bed("None"), bed_hit_dist("None"), bed_radius(1000), write_bedgraph_treatment(false), write_bedgraph_control(false)
	{}
};


seqan::ArgumentParser::ParseResult parseCommandLine(Options &options, int argc, char const ** argv);

#endif
