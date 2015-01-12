/**
 * @file
 * @version 0.01
 *
 * \brief Detect all summits for a given \ref Chromosome.
 *  
 * The detection of summits proceeds in three steps:
 * 	- 1. Find summits in qfrags coverage profile
 * 	- 2. Refine summits
 * 	- 3. Get hit indices for each summut as anchor for evaluation.
 * 
 * All these steps are implemented in a single function called 
 * \ref getSummits
 */
 
#ifndef Q_DETECT_SIGNALS_H
#define Q_DETECT_SIGNALS_H
 
 
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
int getSummits(Chromosome &chromosome, int min, int max, int cutoff);

#endif

