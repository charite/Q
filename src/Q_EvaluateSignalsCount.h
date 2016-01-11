/**
 * @file
 * @authors  Peter Hansen <peter.hansen@charite.de>
 * @version 0.01
 *
 * \brief Contains everthing what has to do with the evaluation of 
 * summits.
 * 
 * After the detection of summits each \ref \Chromosome object has a set
 * of summits \ref Chromosome::SUMMITS. In this unit each summit is 
 * evaluated, that is it is assigned a score.
 * 
 * The evaluation proceeds in two steps, which are still under
 * development.
 * 
 * The first step is to count hits and things similar to hits around
 * each summit. There are several different functions for this which
 * are implemented in \ref Q_EvaluateSignals.h. These functions are:
 * 	-# \ref getQfragEnds
 * 	-# \ref getKs
 * 	-# \ref getFeasibleQFrags
 * 
 * In a next step it is asked how likely it is to observe these counts
 * by chance. Here again there are several functions for calculating
 * p-values or scores which are implemented in 
 * \ref Q_EvaluateSignalsPvalue.h. These functions are:
 * 	-# \ref getLnPvalues1
 * 	-# \ref getLnPvalues2
 * 	-# \ref getLnPvalues3
 * 
 * @todo This unit is the most creepiest part of the program. It has a 
 * lot of functions. Document and check each function. Find a way to
 * compare the differnt scores produced by the functions. Decide which
 * score is the best.
 * 
 */
 
#ifndef Q_EVALUATE_SIGNALS_COUNTS_H
#define Q_EVALUATE_SIGNALS_COUNTS_H

#include "Q_ReadInFiles.h"


/**
 * \brief Given a chromosome for which all summits are already detected.
 * Find for each summit the numbers of qfrags that overlap the summit.
 * 
 * @param[in] &chromosome A \ref Chromosome object
 * @param[in] min Minimal fragment length
 * @param[in] max Maximal fragment length
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * 
 * The function iterates over all summits of a given chromosome. For each
 * summit \f$s_i\f$ at position \f$pos_i\f$ it determines the numbers of
 * qfrags that overlap \f$s_i\f$ for ChIP and control sample:
 * 	-# kfu - Number of forward strand hits in the interval \f$[pos_i-r,pos_i]\f$
 * 	-# kfd - Number of forward strand hits in the interval \f$[pos_i,pos_i+r]\f$
 *  
 * These numbers are stored in \ref Summit::kfu, \ref Summit::kfd, 
 * \ref Summit::kru and \ref Summit::krd.
 * 
 */
int getQFragsCoverage(Chromosome &chromosome, int min, int max);

 
 /**
 * \brief Given a chromosome for which all summits are already detected.
 * Find for each summit the numbers of positions within a given distance
 * before and behind the summit that are covered by at least one qfrag
 * end.
 * 
 * @param[in] &chromosome A \ref Chromosome object
 * @param[in] min Minimal distance \f$min\f$ between two hits forming a qfrag.
 * @param[in] max Maximal distance \f$max\f$ between two hits forming a qfrag.
 * 
 * \returns 0 for sucess and 1 otherwise
 * 
 * The function iterates over all summits of a given chromosome. For each
 * summit \f$s_i\f$ at position \f$pos_i\f$ it determines how many 
 * positions within the interval \f$[pos_i-max,pos_i]\f$ are covered by at 
 * least one hit on the forward strand which is part of a qfrag that
 * covers the summit \f$s_i\f$. Analogous, the number of positions 
 * within the interval \f$[pos_{i+1},pos_{i+1}+max]\f$ that are covered
 * by at least one hit on the reverse strand which is part of a qfrag 
 * that covers the summit \f$s_i\f$ is determined.
 * 
 * These numbers are stored in \ref Summit::q_beg and \ref Summit::q_end,
 * respectively. 
 * 
 * @todo This function has a misleading name. Not the qfrag ends are
 * determined but the numbers of saturated positions.
 * 
 * 
 */
int getQfragEnds(Chromosome &chromosome, int min, int max);

int getKs(Chromosome &chromosome, int radius);

#endif
