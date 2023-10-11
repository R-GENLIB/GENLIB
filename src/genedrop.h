#ifndef GENEDROP
#define GENEDROP

#include <RcppCommon.h>

void gene_drop(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
				double* probRecomb, double* Morgan_Len, int BP_len, int model,  int nSimul,
				int convert, double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, double* R_matrix, 
				const std::string& out, const std::string& map_filepath, const std::string& ped_filepath, int seed);

// void pIBD_matrix(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
// 				double* probRecomb, double* Morgan_Len, int BP_len, int model,
// 				int convert, double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, double* R_matrix, 
// 				const std::string& out, const std::string& map_filepath, const std::string& ped_filepath, int seed);

// void pIBD_matrix2(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
// 				double* probRecomb, double* Morgan_Len, int BP_len, int model,
// 				int convert, double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, double* R_matrix, int seed);

#endif