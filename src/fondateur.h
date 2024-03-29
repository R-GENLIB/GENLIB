/*! \file fondateur.h
\brief Interface des fonctions de simulation calcul de probabilite

Interface de toutes les fonctions en rapport avec le gene fondateur

\author Sbastien Leclerc 
\contributor Jean-Franois Lefebvre

*/

#ifndef GENFOND
#define GENFOND

#include <RcppCommon.h>
#include <unordered_map>

int simul(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		int lSimul, double* pdRetConj,double* pdRetSimul,double* pdRetProp,double* probRecomb,double probSurvieHomo,int printprogress);

void simulhaplo(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
						int lSimul, double* probRecomb, double* Morgan_Len, int BP_len, int model, int convert,  
						double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, 
						std::unordered_map<int,haplotype*> *hapRef, std::string WD, int write_all_node, int seed);

int simulhaplo_traceback(std::string& path_ANH, std::string& path_PH, int& myPro, int& myAnc, 
				std::vector<int>& indVec, std::vector<int>& mererVec, std::vector<int>& pereVec,
				std::vector<int>& resultvec1, std::vector<int>& resultvec2, std::vector<int>& resultvec3);

void simulhaplo_compare_IBD (const int& pro1_ID, const int& pro2_ID, const int& BP_len, std::string& file_path, std::vector<int>& rvec1, std::vector<int>& rvec2, std::vector<double>& rvec3, std::vector<int>& rvec4);
// int getNumberRec(double* probRecomb, int sex);
// double getRandomNumber(int exponential);
//int descendreHaplotypes(CIndSimul* Ordre_tmp, double probHap); //, /**/std::unordered_map<std::string, haplotype*>/* const std::unordered_map<std::string, haplotype*> &*/*hapRef);
//void makeRecomb( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, double posRecomb, int& cle );
void no_convert( int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array);
void convert1(   int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array);

void makeRecombM( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int* posRecomb, int& cle );
void makeRecombF( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int* posRecomb, int& cle );

void recombine( haplotype* hapBegin, haplotype* hapEnd, haplotype* hapChild, int nbRecomb, int *posRecomb );

bool reconstruct(const std::string &WD);
bool ancestralseq(const std::string &fileName, std::unordered_map<float, std::string> &haploseqs);
std::vector<int> readSNPpos(const std::string &fileName);

int simulsingle(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
			 int lSimul, double* pdRetour,int printprogress);

int simulsingleFreq(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
				int lSimul, double* pdRetour,int printprogress);

SEXP simulsingleFct(int* SGenealogie, int * proposant, int lproposant, int* SplAncetre, int* SplAncEtatAll1, int* SplAncEtatAll2, int SlNAncetre,
				int SlSimul, SEXP SfctSousGrp, int Sprintprogress);

SEXP simulsingleProb(int* SGenealogie, int* SplProposant, int SlNProposant, int* SplAncetre,int SlNAncetres, int* SplAncEtat,SEXP mtProb,
				 int SlSimul, int Sprintprogress);

SEXP prob(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
	    double* pdRetConj,double* pdRetSimul,int printprogress,int onlyConj);

int CoefApparentement(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre,	double* pdRetour,int DuppDetection, int printprogress);

#endif



