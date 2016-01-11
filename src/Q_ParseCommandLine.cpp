#include "Q_ParseCommandLine.h"

/**
 * \brief Function which parses the command line arguments and stores
 * the results into an \ref Options struct.
 * 
 * This fuction uses the 
 * <a href="http://trac.seqan.de/wiki/Tutorial/ParsingCommandLineArguments">command line parser of seqan</a>.
 * 
 * At first the ArgumentParser is set up, that is all possible arguments
 * and option are defined. Furthermore the help screen is configured.
 * 
 * After parsing the command line and catching some exceptions the
 * content of the \ref Options object is written to a file called
 * <PREFIX>-Q-info.txt.
 */
seqan::ArgumentParser::ParseResult parseCommandLine(Options &options, int argc, char const ** argv)
{
	
	// Setup ArgumentParser
	// --------------------
	
	seqan::ArgumentParser parser("Q");
	setShortDescription(parser, "Saturation based ChIP-seq peak caller");
	setVersion(parser, "1.3.0");
	setDate(parser, "January 2016");
	addUsageLine(parser,"[\\fIOPTIONS\\fP] --treatment-sample [\\fIinput-file\\fP] --out-prefix [\\fISTRING\\fP]");
	addDescription(parser,"Q is a fast saturation-based ChIP-seq peak caller. Q works well in conjunction with the irreproducible discovery rate (IDR) procedure. Q was extensively tested on publicly available datasets from ENCODE and shown to perform well with respect to reproducibility of the called peak set, consistency of the peak sets with respect to predicted transcription factor binding motifs contained in them, as well as overall run time. Q is implemented in C++ making use of the Boost and SeqAn library. There are a number of useful features for the primary analysis of ChIP-seq data. Q can be run with or without data from a  control experiment. Duplicate reads are removed on the fly without altering the original BAM file, and the number of duplicated reads is then shown in Q's output. The average fragment length of the sequencing library, which is an essential parameter for peak calling and for downstream analysis, is estimated automatically from the data. This is done by examing the vector of read start positions along individual chromosomes and calcuting the shift that is associated with the smallest Hamming distance. This procedure yields an equivalent estimation of the average fragment length as the cross-correlation plot of SPP but is approximately three times faster. As a part of this procedure, Q also calculates the relative strand cross-correlation coefficient (RSC), which allows a global quality assessment of the enrichment. In addition Q offers its own quality metrics, which can be used for trouble-shooting and quality control of the results. If desired, Q also generates fragment coverage profiles which can be uploaded to UCSC's genome browser, where they can be displayed in the context of other related data such as for example ChIP-seq data for histone modifications and cofactors or expression data.");
	

	// Define options
	// --------------

	addSection(parser, "Peak Calling Options");
	
	addOption(parser, seqan::ArgParseOption(
		"t", "treatment-sample", "Input file (REQUIRED).",
		seqan::ArgParseArgument::INPUTFILE, "IN"));
	setRequired(parser, "t");
	//setValidValues(parser, "treatment-sample", "sam bam");

	addOption(parser, seqan::ArgParseOption(
		"c", "control-sample", "Input file.",
		seqan::ArgParseArgument::INPUTFILE, "IN"));
	//setValidValues(parser, "control-sample", "sam bam");

	addOption(parser, seqan::ArgParseOption(
		"l", "fragment-length-average", "The average length of fragments in the sequencing library of the treatment sample. If not given this value will be determined from the treatment data via binding characteristics.",
	seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "l", -1);

	addOption(parser, seqan::ArgParseOption(
		"x", "fragment-length-deviation", "The deviation of lengths of fragments in the sequencing library.",
	seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "x", 50);

	addOption(parser, seqan::ArgParseOption(
		"k", "keep-dup", "If set, duplicate reads will be kept."));
	
	addOption(parser, seqan::ArgParseOption(
		"n", "top-n", "Maximum number of top peaks to be written to file.",
	seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "n", 100000);
	
	addOption(parser, seqan::ArgParseOption(
		"pc", "p-value-cutoff", "Cutoff for the negative decadic logarithm of the adjusted p-values (i.e. 6 means 1e-06). Overrides --top-n.",
		seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
	setDefaultValue(parser, "pc", -1);

	addOption(parser, seqan::ArgParseOption(
		"nm", "nexus-mode", "If set, appropriate settings for ChIP-nexus will be used. Duplicate reads will be kept. If not set, the fragment length l will be estimated using the qfrag-length-distribution method and x will be set to 5."));
		
	addSection(parser, "General Options");

	addOption(parser, seqan::ArgParseOption(
		"o", "out-prefix", "Prefix for all output files.",
		seqan::ArgParseArgument::STRING, "STRING"));
	setDefaultValue(parser, "o", "OUT");		
	
	addOption(parser, seqan::ArgParseOption(
		"p", "thread-num", "Number of threads to use.",
		seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "p", 1);
			
	addOption(parser, seqan::ArgParseOption(
		"v", "verbose", "If set, information will be written to the screen."));
	
	addSection(parser, "Binding Characteristics Options For Fragment Length Estimation");
		
	addOption(parser, seqan::ArgParseOption(
		"s", "step-num", "Number of strand shifts for the determination of the average fragment length.",
		seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "s", 1000);
	
	addOption(parser, seqan::ArgParseOption(
		"bco", "binding-characteristics-only", "If set, only the binding characteristics will be determined and peak calling will be skipped."));

	addOption(parser, seqan::ArgParseOption(
		"qld", "qfrag-length-distribution", "If set, the distribution of qfrag lengths will be determined and peak calling will be skipped."));

	addSection(parser, "UCSC Track Options");
	
	addOption(parser, seqan::ArgParseOption(
		"wbt", "write-bedgraph-treatment", "If set, bedGraph file for the treatment sample will be genrated. Reads will be extended to the given or estimated average fragment length."));

	addOption(parser, seqan::ArgParseOption(
		"wbc", "write-bedgraph-control", "If set, bedGraph file for the control sample will be genrated. Reads will be extended to the given or estimated average fragment length of the treatment data."));

	addSection(parser, "Advanced Options");
				
	addOption(parser, seqan::ArgParseOption(
		"w", "write-bed", "Chromosome ID (e.g. chr1) must be consistent with ID in SAM/BAM file. If set, BED files for reads, shifted reads, extended reads and qfrags will be written for the given chromosome ID. Peak calling will be skipped.",
		seqan::ArgParseArgument::STRING, "STRING"));

	addOption(parser, seqan::ArgParseOption(
		"b", "bed-hit-dist", "Input BED file containing summits. Summit position is always the center of a given region. Chromosome IDs in BED file must be consistent with IDs in SAM/BAM file. Default radius around summits is 1000. The radius can be changed via the -r option. Distribution of hits around summits on forward and reverse strand will be written to a text file. Output is a tab separated table containing three columns and 2 times radius rows. The first column contains the relative positions to the summits. The second and third column contain the accumulated counts of hits for all summits in the BED file. Second column for forward and third column for reverse strand. Peak calling will be skipped.",
		seqan::ArgParseArgument::INPUTFILE, "IN"));

	addOption(parser, seqan::ArgParseOption(
		"r", "bed-radius", "Radius around summits for counting hits (-b).",
		seqan::ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "r", 1000);		

	addOption(parser, seqan::ArgParseOption(
		"psc", "use-pseudo-control", "If set, a pseudo control will be generated from the treatment data, by switching the strand of each read and shifting the 5' end towards 3' direction by one read length."));
	

	// Add a section for the description of the outout to help screen
	// --------------------------------------------------------------

	addTextSection(parser, "Output");
	addListItem(parser,
				"\\fBout-prefix-Q-runinfo.txt\\fP",
				"This file contains all chosen parameters.");
	addListItem(parser,
				"\\fBout-prefix-Q-binding-characteristics.R\\fP",
				"This file contains R code for generating the binding characteristics plot in pdf format. It will only be generated, if no average fragment length is specified.");					
	addListItem(parser,
				"\\fBout-prefix-Q-chromosome-info.tab\\fP",
				"This file contains information for each chromsome.");
	addListItem(parser,
				"\\fBout-prefix-Q-narrowPeak.bed\\fP",
				"This file contains information for each summit in ENCODE narrowPeak.");	
	addListItem(parser,
				"\\fBout-prefix-Q-summit-info.tab\\fP",
				"This file contains more detailed information for each summit.");	
		
	// Add examples command lines to help screen
	// -----------------------------------------

	addTextSection(parser, "Examples");
	addListItem(parser,
				"\\fBQ\\fP \\fB--treatment-sample\\fP \\fIinput-file\\fP \\fB--out-prefix\\fP \\fISTRING\\fP",
				"Minimal call using default parameters and no control sample. An average fragment length will be estimated from the ChIP data and used for peak calling.");
	addListItem(parser,
				"\\fBQ\\fP \\fB--treatment-sample\\fP \\fIinput-file\\fP \\fB--control-sample\\fP \\fIinput-file\\fP \\fB--out-prefix\\fP \\fISTRING\\fP",
				"Minimal call using a control sample and default parameters. An average fragment length will be estimated from the ChIP data and used for peak calling.");
	addListItem(parser,
				"\\fBQ\\fP \\fB-t\\fP \\fIinput-file\\fP \\fB-c\\fP \\fIinput-file\\fP \\fB-o\\fP \\fISTRING\\fP \\fB--fragment-length-average\\fP \\fIINT\\fP",
				"Minimal call using a control sample and default parameters. The specified average fragment length will be used for peak calling.");
	addListItem(parser,
				"\\fBQ\\fP \\fB--treatment-sample\\fP \\fIinput-file\\fP \\fB--out-prefix\\fP \\fISTRING\\fP \\fB--binding-characteristics-only\\fP",
				"Only the binding characteristics will be determined and peak calling will be skipped.");

		
	// Parse command line
	// ------------------    

	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);


	// Extract options if there were no exceptions during parsing
	// ----------------------------------------------------------

	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res;
	}

	getOptionValue(options.chip_sample, parser, "treatment-sample");

	getOptionValue(options.control_sample, parser, "control-sample");

	if (!isSet(parser,"control-sample"))
		options.control_sample="None";

	options.use_pseudo_control = isSet(parser, "use-pseudo-control");

		
	options.keep_dup = isSet(parser, "keep-dup");		
		
	getOptionValue(options.out_prefix, parser, "out-prefix");
	
	getOptionValue(options.fragment_length_avg, parser, "fragment-length-average");
	
	getOptionValue(options.fragment_length_dev, parser, "fragment-length-deviation");
	
	getOptionValue(options.top_n, parser, "top-n");
	
	getOptionValue(options.p_value_cutoff, parser, "p-value-cutoff");

	options.nexus_mode = isSet(parser, "nexus-mode");		
	
	getOptionValue(options.step_num, parser, "step-num");
	
	options.binding_characteristics_only = isSet(parser, "binding-characteristics-only");
	
	options.make_qfrag_length_distribution = isSet(parser, "qfrag-length-distribution");


	getOptionValue(options.thread_num, parser, "thread-num");

	options.verbose = isSet(parser, "verbose");
	
	options.write_bedgraph_treatment = isSet(parser, "write-bedgraph-treatment");
	
	options.write_bedgraph_control = isSet(parser, "write-bedgraph-control");	
		
	getOptionValue(options.write_bed, parser, "write-bed");
	
	getOptionValue(options.bed_hit_dist, parser, "bed-hit-dist");
	getOptionValue(options.bed_radius, parser, "bed-radius");

	options.out_info_file=options.out_prefix;
	append(options.out_info_file,"-Q-runinfo.txt");
	
	options.binding_characteristics_file=options.out_prefix;
	append(options.binding_characteristics_file,"-Q-binding-characteristics.R");

	options.out_chromosome_info_file=options.out_prefix;
	append(options.out_chromosome_info_file,"-Q-quality-statistics.tab");

	options.out_summit_info_file=options.out_prefix;
	append(options.out_summit_info_file,"-Q-summit-info.tab");

	options.out_narrowPeak_file=options.out_prefix;
	append(options.out_narrowPeak_file,"-Q-narrowPeak.bed");
	
	if(options.write_bedgraph_treatment)
	{
		options.out_bedgraph_treatment_file=options.out_prefix;
		append(options.out_bedgraph_treatment_file,"-Q-treatment.bedgraph");
	}
	
	if(options.write_bedgraph_control)
	{
		options.out_bedgraph_control_file=options.out_prefix;
		append(options.out_bedgraph_control_file,"-Q-control.bedgraph");
	}
	
	// If everything is ok, write options to file and return
	// -----------------------------------------------------
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime ); 
	timeinfo = localtime ( &rawtime );
	
	std::ofstream INFO; 
	INFO.open (toCString(options.out_info_file));
	INFO      << "Q - RUN INFO" << '\n' << "############\n" << '\n'

			  << "General information" << '\n' << "===================\n"
			  << "Q Version:              " << getVersion(parser) << '\n'
			  << "Authors:                " << "Peter Hansen <peter.hansen@charite.de>" << '\n'
			  << "                        " << "Peter Nick Robinson <peter.robinson@charite.de>" << '\n'                               
			  << "Execution time:         " << asctime(timeinfo) << '\n'

			  << "Input" << '\n' << "=====\n"
			  << "treatment-sample:       " << options.chip_sample << '\n'
			  << "control-sample:         " << options.control_sample << '\n'			  
			  << "bed-hit-dist:           " << options.bed_hit_dist << '\n' << '\n'
			  
			  << "Peak Calling Arguments" << '\n' << "======================\n";

			  if(options.nexus_mode)
			  {              
				INFO << "nexus-mode:                        " << options.nexus_mode << '\n';
				options.keep_dup=1;
				options.fragment_length_dev=5;
				if(options.fragment_length_avg!=-1)
				{
					INFO << "fragment-length-average:           " << options.fragment_length_avg << '\n';
				}
				else
				{
					INFO << "fragment-length-average:           " << "not specified, will be determined by qfrag-length-distribution plot" << '\n';
				}
			  }
			  
			  if(options.fragment_length_avg!=-1 && !options.nexus_mode)
			  {              
				INFO << "fragment-length-average:           " << options.fragment_length_avg << '\n';
			  }
			  else if(!options.nexus_mode)
			  {
				INFO << "fragment-length-average:           " << "not specified, will be determined by Hamming distance plot" << '\n';
			  }
	INFO      << "fragment-length-deviation:         " << options.fragment_length_dev << '\n'
			  << "keep-dup:                          " << options.keep_dup << '\n'			  
			  << "top-n:                             " << options.top_n << '\n';
			  if(options.p_value_cutoff!=-1)
			  {
				  INFO << "p-value-cutoff:                    " << options.p_value_cutoff << '\n' << '\n';
		      }
		      else
		      {
				  INFO << "p-value-cutoff:                    " << "not specified, will use top-n" << '\n' << '\n';
			  }
			  
INFO		  << "General Arguments" << '\n' << "=================\n"
			  << "out-prefix:                        " << options.out_prefix << '\n'                         			  
			  << "thread-num:                        " << options.thread_num << '\n'                         
			  << "verbose:                           " << options.verbose << '\n' << '\n'
			  
			  << "Binding Characteristics Arguments For Fragment Length Estimation" << '\n' << "================================================================\n" 			  
			  << "step-num:                          " << options.step_num << '\n'
			  << "binding-characteristics-only:      " << options.binding_characteristics_only << '\n'
			  << "qfrag-length-distribution:      " << options.make_qfrag_length_distribution << '\n' << '\n'
			  
			  << "UCSC Track Options" << '\n' << "==================\n" 			  
			  << "write-bedgraph-treatment:          " << options.write_bedgraph_treatment << '\n'
			  << "write-bedgraph-control:            " << options.write_bedgraph_control   << '\n' << '\n'
			  			  
			  << "Advanced Options" << '\n' << "================\n"
			  << "write-bed:                         " << options.write_bed << '\n'                         			  
    		  << "bed-radius:                        " << options.bed_radius << '\n' << '\n'
			  
			  << "Output" << '\n' << "======\n"
			  << "runinfo:                           " << options.out_info_file << '\n';
			  
	if(options.fragment_length_avg==-1)
	{
		INFO  << "binding-characteristics:           " << options.binding_characteristics_file << '\n';
	}		  
		
	if(!options.binding_characteristics_only)
	{  
		INFO  << "quality-statistics:                " << options.out_chromosome_info_file << '\n'
			  << "summit-info:                       " << options.out_summit_info_file << '\n'
			  << "narrowPeak:                        " << options.out_narrowPeak_file << '\n';
		if(options.write_bedgraph_treatment)
		{
			INFO << "bedGraph-treatment:                " << options.out_bedgraph_treatment_file << '\n';
		}
		if(options.write_bedgraph_control)
		{
			INFO << "bedGraph-control:                  " << options.out_bedgraph_control_file << '\n';
		}	
	}
	else
	{
		INFO  << '\n';
	}

	INFO  << '\n';	  
	INFO.close();

	return seqan::ArgumentParser::PARSE_OK;
}
