/**
 * @file
 * @authors  Peter Hansen <peter.hansen@charite.de>
 * @version 0.01
 *
 * \brief Estimate fragment length via cross-correlation
 *  
 * XXX
 *  
 */

#ifndef GET_FRAGMENT_LENGTH_H
#define GET_FRAGMENT_LENGTH_H

int getFragLength(std::vector<Chromosome> &chromosome, int step_num, int thread_num, std::string out_prefix);

int getQfragLengthDistribution(std::vector<Chromosome> &chromosome, int step_num, int thread_num, std::string out_prefix);

std::vector<Hit> getPseudoControlSwitchStrandAndFlip(std::vector<Hit> CHIP_HITS, int read_length, int chr_length); 

int cleanChipHits(Chromosome &chromosome); 


#endif
