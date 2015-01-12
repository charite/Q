/**
 * @file
 * @authors  Peter Hansen <peter.hansen@charite.de>
 * @version 0.01
 *
 * \brief Contains everything what has to do with writing the output. 
 * 
 */
 

 /**
 * \brief Function writes global information for each chromosome to a
 * tab separated text file.
 * 
 * @param[in] chromosome_info_file Filename
 * @param[in] vector<Chromosome>&chromosome A vector with all \ref Chromosome objects
 * @param[in] chr_num Total number of chromosomes
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * The output file has one row for each chromosome and 7 columns
 * 	-# Name of the chromosome
 * 	-# Length of the chromosome
 * 	-# Number of forward hits for the chromosome
 * 	-# Number of reverse hits for the chromosome
 * 	-# Total number of hits for the chromosome
 * 	-# Total number of qfrags for the chromosome (given min and max)
 * 	-# Total number of summits for the chromosome
 * 
 * There is an additional first row containing the names of the columns.
 * 
 */
int writeChromosomeInfo(seqan::CharString chromosome_info_file, std::vector<Chromosome> &chromosome, int chr_num, int q_min, int q_max, bool keep_dup);
int writeChromosomeInfoShort(seqan::CharString out_prefix, seqan::CharString chromosome_info_file, std::vector<Chromosome> &chromosome, int chr_num, int q_min, int q_max, bool keep_dup,int fragment_length_avg);
int getObsNumQHits(std::vector<Hit> HITS);
long double getExpNumQHits(long int TF, long int TR, long int L, int q_min, int q_max);
double getQEnrichmentScore(int obs_q_hit_num, double exp_q_hit_num, int UMNR_hit_num);
int writeBedGraph(std::vector<Chromosome> &chromosome, seqan::CharString out_prefix, int extension_length, bool is_control);

 /**
 * \brief Function writes all information for all summits to a tab
 * separated text file.
 * 
 * @param[in] summit_info_file Filename
 * @param[in] vector<Chromosome>&chromosome A vector with all \ref Chromosome objects
 * @param[in] chr_num Total number of chromosomes
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * The output file has one row for each summit and 11 columns
 * 	-# Name of the chromosome
 * 	-# Position
 * 	-# Position + 1
 * 	-# qfrag coverage at summit position
 * 	-# log10 p-value
 * 	-# kfu
 * 	-# kfd
 * 	-# kru
 * 	-# krd
 * 	-# qbeg
 * 	-# qend
 * 
 */
int writeSummitInfo(seqan::CharString summit_info_file, std::vector<Summit_chr> &summit, Options &options);

 /**
 * \brief Function writes all summits to a tab separated text file 
 * <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format12">in
 * narrowPeak format</a>
 * 
 * @param[in] narrowPeak_file Filename
 * @param[in] vector<Chromosome>&chromosome A vector with all \ref Chromosome objects
 * @param[in] chr_num Total number of chromosomes
 * @param[in] radius Radius of the summit
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * The output file has one row for each summit and 10 columns
 * 	-# chrom
 * 	-# chromStart (summit position)
 * 	-# chromEnd (summit position + 1)
 * 	-# name (.)
 * 	-# score (qfrag coverage at summit position)
 * 	-# strand (.)
 * 	-# signalValue (??? see source file)
 * 	-# pValue (log10_pvalue)
 * 	-# qValue (-1)
 * 	-# peak (radius)
 * 
 */
int writeNarrowPeak(seqan::CharString narrowPeak_file, std::vector<Summit_chr> &summit, int radius, Options &options);

 /**
 * \brief Function for testing. Given a fragment length and chromosome ID,
 * this function writes BED files for plain reads, shifted reads,
 * extended reads _AND_ qfrags. 
 * 
 */
int writePlainShiftedExtendedReadsAndQfragsToBED(std::vector<Chromosome> &chromosome, Options &options);


 /**
 * \brief Function for testing writes BED files for mapped reads,
 * shifted and extended reads, for a given chromosome and an
 * average fragment length.
 * 
 * Input is
 * 	-# a chromosme object
 * 	-# a prefix for the output files, which have the suffixes
 * '_chip_hits.bed' and '_control_hits.bed', respectively.
 * 
 */
int writeReadsToBED(Chromosome &chromosome, std::string out_prefix, int extension_length, int shift_size);



 /**
 * \brief Function for testing writes all qfrags for given chromosome to
 * a BED file. Two BED files are created, one for the ChIP and one for
 * the control sample.
 * 
 * Input is
 * 	-# a chromosme object
 * 	-# a prefix for the output files, which have the suffixes
 * '_chip.bed' and '_control.bed', respectively.
 * 
 */
int writeQFragsToBED(Chromosome &chromosome, std::string out_prefix, int min, int max);

 /**
 * \brief Function for testing. Determines the cumulative absolute
 * frequencies for hits in a radius around all summits.
 * 
 * Output is a table with '2*radius+1' columns and three rows.
 * Each column stands for one position. The first row contains
 * the positions realtive to the summits, that is
 * 
 * radius,...,0,...,radius
 * 
 * For a given summit this corresponds to the counts at positions
 * 
 * summit_position-radius,...,summit_position,...,summit_position+radius
 * 
 * The second row contains the corresponding cumulated forward strand 
 * hit counts for each position. Whereas the second row contains
 * the cumulated reverse strand hit counts.
 * 
 * Input is
 * 	-# a chromosme object
 * 	-# a radius
 * 	-# a prefix for the output files, which have the suffixes
 * '_chip_qfrags.bed' and '_control_qfrags.bed', respectively.
 * 
 */
int writeHitDist(std::vector<Chromosome> &chromosome, int radius, std::string out_prefix);

