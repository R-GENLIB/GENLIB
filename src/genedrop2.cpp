// #include <vector>
// #include <unordered_map>
// #include <fstream>
// #include <random>
// #include <algorithm>
// #include <Rcpp.h>

// #include "asa239.h"
// #include "base.h"
// #include "outils.h" 
// #include "cbignum.h"
// #include "hashtable.h"
// #include "userInterface.h"


// // ******************************************************************** 
// //
// //			Local func & class
// //
// // ********************************************************************

// typedef std::unordered_map<int, std::vector<char>> genotype_map;
// static void makeRecombF( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle, const int& BP_len );
// static void makeRecombM( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle, const int& BP_len );
// static void recombine(haplotype* hapBegin, haplotype* hapEnd, haplotype* hapChild, int nbRecomb, int* posRecomb, const int& BP_len);
// static void no_convert(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array);
// static void convert1(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array);
// static void read_ped_file(const std::string& ped_filepath, genotype_map& geno_map, const int& n_markers);
// static void read_map_file(const std::string& map_filepath, std::vector<int>& vec_pos);
// static void convert_hap(char* marker_vec, haplotype* hap, const genotype_map& founder_geno_map, const std::vector<int>& map);


// //should make static, and then store the crossover vector as a data member
// // instead of using a fx pointer to the member fx then make the gen_crossovers fx
// //     a member fx pointer to the appropriate model at time of initialization
// class Crossovers
// {
// 	public:		
// 		void Poisson_CO(const int&, double*, double*, int&, std::mt19937&, double*);
// 		void Poisson_ZT(const int&, double*, double*, int&, std::mt19937&, double*);
// 		void init_gamma(double&, double&, double&, double& );
// 		void Gamma_CO  (const int&, double*, double*, int&, std::mt19937&, double*);

//  	private:
// 		double first_arrival[2][10000]; //Numbins = 10,000. Need to calculate first arrival times for gamma process. They are distributed differently than the interrarival times
// };

// void Crossovers::Poisson_CO(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
// 	double pos;			
// 	static std::uniform_real_distribution<> u_dist(0, 1);
// 	static std::poisson_distribution<int> p1_dist(param[0]);
// 	static std::poisson_distribution<int> p2_dist(param[1]);

// 	if(sex==1){
// 		NumRecomb = p1_dist(mtgen); //number of recombination events 
// 		for( int h=0; h<NumRecomb; h++ ) {
// 			pos = u_dist(mtgen);
// 			CO_array[h] = pos; 
// 		}
		
// 		std::sort(CO_array, CO_array + NumRecomb);
// 	}
// 	else{
// 		NumRecomb = p2_dist(mtgen); //number of recombination events 
// 		for( int h=0; h<NumRecomb; h++ ) {
// 			pos = u_dist(mtgen);
// 			CO_array[h] = pos; 
// 		}
// 		std::sort(CO_array, CO_array + NumRecomb);		
// 	}
// }

// void Crossovers::Poisson_ZT(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
// 	double chiasma_pos[40];
// 	double pos;
// 	int	   num_chiasma = 0;		
// 	static std::uniform_real_distribution<> u_dist(0, 1);
// 	static std::poisson_distribution<int> p1_dist(param[0]);
// 	static std::poisson_distribution<int> p2_dist(param[1]);

// 	if(sex==1){
// 		num_chiasma = p1_dist(mtgen); //number of recombination events 
// 		while (num_chiasma == 0){
// 			num_chiasma = p1_dist(mtgen);
// 		}
		
// 		for( int h=0; h<num_chiasma; h++ ) {
// 			pos = u_dist(mtgen);
// 			chiasma_pos[h] = pos; 
// 		}
		
// 	}
// 	else{
// 		num_chiasma = p2_dist(mtgen); //number of recombination events 
// 		while (num_chiasma == 0){
// 			num_chiasma = p2_dist(mtgen);
// 		}
		
// 		for( int h=0; h<num_chiasma; h++ ) {
// 			pos = u_dist(mtgen);
// 			chiasma_pos[h] = pos; 
// 		}
// 	}
	
// 	int counter=0;
// 	double u_rand;
// 	NumRecomb = 0;
// 		for (int i=0; i<num_chiasma; i++){
// 			u_rand = u_dist(mtgen);
// 			if (u_rand < 0.50){
// 				CO_array[counter] = chiasma_pos[i];
// 				counter   += 1;
// 				NumRecomb += 1;
// 			}
// 		}
	
// 	std::sort(CO_array, CO_array + NumRecomb);
// }

// void Crossovers::init_gamma(double& paramF,  double& paramM,  double& Morgan_LenF,  double& Morgan_LenM){
// 	double x;
// 	int wasteman = 0;

// 	//Number of bins is 10,000, these are the bins for reimann sum
// 	//Morgan_Len/10000 is the delta x of the integral
// 	//gamma_q is 1-CDF. We multiply by 2 because we are dividing by mean of regular distribution (which has mean 1/2)
	
// 	x= Morgan_LenF/10000;
// 	first_arrival[0][0] = 2 * (1-gammad(2*paramF*x, paramF, &wasteman)) * Morgan_LenF/10000;

// 	for ( int i=1; i<10000; i++){
// 		x = Morgan_LenF*(i+1)/10000;
// 		first_arrival[0][i] = 2 * (1-gammad(2*paramF*x, paramF, &wasteman)) * Morgan_LenF/10000 + first_arrival[0][i-1];
// 	}	

// 	x= Morgan_LenM/10000;
// 	first_arrival[1][0] = 2 * (1-gammad(2*paramM*x, paramM, &wasteman)) * Morgan_LenM/10000;

// 	for ( int i=1; i<10000; i++){
// 		x = Morgan_LenM*(i+1)/10000;
// 		first_arrival[1][i] = 2 *  (1-gammad(2*paramM*x, paramM, &wasteman)) * Morgan_LenM/10000 + first_arrival[1][i-1];
// 	}
// }

// void Crossovers::Gamma_CO(const int &sex, double *param, double *Morgan_len, int& NumRecomb, std::mt19937& mtgen, double *CO_array){
// 	double u_rand;
// 	double length;
// 	double step;
// 	double interrarival;
// 	double current_pos;
// 	double chiasma_pos[20];
// 	int Num_Chiasmata = 0;
	
// 	static std::uniform_real_distribution<> u_dist(0, 1);
// 	static std::gamma_distribution<> g1_dist(param[0],1/(2*param[0]));
// 	static std::gamma_distribution<> g2_dist(param[1],1/(2*param[1]));

// 	if (sex==1){
// 		length = Morgan_len[0];
// 	} else{
// 		length = Morgan_len[1];
// 	}

// 	step = length/10000;
//     u_rand = u_dist(mtgen); 

// 	if ( u_rand > first_arrival[sex-1][9999]) NumRecomb = 0; //IF the first chiasmata is beyond the length of the chromsome then we have no recombination
// 	else {
// 		if (u_rand <= first_arrival[sex-1][0]){ 
// 			chiasma_pos[0] = 0.5 * step;
// 			Num_Chiasmata = 1;
// 		}
// 		else { //binary search to find the first position (sample between 0-1, then lookup the inverse CDF values)
// 			int mid, low = 0, high = 10000;
// 			while (high - low > 1) {
// 				mid = (high - low) / 2 + low;
// 				if (u_rand <= first_arrival[sex-1][mid])
// 					high = mid;
// 				else if (u_rand > first_arrival[sex-1][mid])
// 					low = mid;
// 			}
// 			Num_Chiasmata=1;
// 			//mid will end up being the index of the smallest value (step of reimann partial sums) we are less than or equal to
// 			chiasma_pos[0] = mid*step + 0.5*step;
// 		}

// 		current_pos = chiasma_pos[0];
// 		//After finding the first arrival time we can  use std::gamma_distribution for the rest
// 		if (sex==1){
// 			interrarival = g1_dist(mtgen);
// 			int counter = 1;
// 			while (current_pos + interrarival < Morgan_len[0]){
// 				Num_Chiasmata = Num_Chiasmata + 1;

// 				chiasma_pos[counter] = current_pos + interrarival;
// 				current_pos 		 = current_pos + interrarival;

// 				counter 	 = counter + 1;
// 				interrarival = g1_dist(mtgen);
// 			}

// 		}
// 		else{
// 			interrarival = g2_dist(mtgen);
// 			int counter = 1;
// 			while (current_pos + interrarival < Morgan_len[1]){
// 				Num_Chiasmata = Num_Chiasmata + 1;

// 				chiasma_pos[counter] = current_pos + interrarival;
// 				current_pos 		 = current_pos + interrarival;

// 				counter 	 = counter + 1;
// 				interrarival = g2_dist(mtgen);
// 			}
// 		}

// 		//After determining position of all chiasmata, we select each one with probability 0.5 of resolving as a crossover 
// 		int counter=0;
// 		NumRecomb = 0;
// 		for (int i=0; i<Num_Chiasmata; i++){
// 			u_rand = u_dist(mtgen);
// 			if (u_rand < 0.50){
// 				CO_array[counter] = chiasma_pos[i]/length;
// 				counter   += 1;
// 				NumRecomb += 1;
// 			}
// 		}
// 	}
// }

// static void read_map_file(const std::string& map_filepath, std::vector<int>& vec_pos){
//     std::ifstream read_map(map_filepath);
//     int pos;
//     std::string x,y,z;
//     while(read_map >> x >> y >> z >> pos) vec_pos.push_back(pos);    
// };

// static void read_ped_file(const std::string& ped_filepath, genotype_map& geno_map, const int& n_markers){
//     std::ifstream read_ped(ped_filepath);
//     std::string line;
    
//     int char_pos[7], IID;
//     char_pos[0] = 0;
// 	char *hap1, *hap2;
// 	int x = 0;
// 	geno_map.emplace(x, std::vector<char>('-', n_markers));
//     while(std::getline(read_ped, line)){
//         for (int i=1; i<7; i++){
//             char_pos[i] = line.find('\t', char_pos[i-1] + 1);
//         }        
//         IID = std::stoi(line.substr(char_pos[1]+1, char_pos[2]));
//         hap1 = (geno_map.emplace(IID, std::vector<char>(n_markers))).first->second.data();
//         hap2 = (geno_map.emplace(-1 * IID, std::vector<char>(n_markers))).first->second.data();
//         for (int i=0; i<n_markers; i++){
//             hap1[i] = line[char_pos[6] + 1 + 4*i];
//             hap2[i] = line[char_pos[6] + 1 + 4*i + 2];
//         }
//     }
// };

// int check_overlaps(const int& Lpos1, const int& Rpos1, const int& Lpos2, const int& Rpos2){
// 	// then its non-overlapping
// 	if ((Lpos1 >= Rpos2) || (Lpos2 >= Rpos1)) return 0;
// 	else return ((Rpos1 > Rpos2 ? Rpos2 : Rpos1) - (Lpos1 > Lpos2 ? Lpos1 : Lpos2) );
// }

// //returns the BP length of IBD shared segments for the pair
// int check_p_IBD(haplotype* hap1,  haplotype* hap2, const int& BP_len){
// 	haplotype* tmp1 = hap1;
// 	haplotype* tmp2 = hap2;
// 	int Lpos1=0, Lpos2=0, Rpos1, Rpos2;

// 	//total amount of IBD sharing (in BP) for the pair of haplotypes
// 	int total_IBD = 0; 
// 	while ((tmp1->pos != BP_len) | (tmp2->pos != BP_len)){
// 		Rpos1 = tmp1->pos;
// 		Rpos2 = tmp2->pos;

// 		if (tmp1->hap == tmp2->hap) {total_IBD = total_IBD + check_overlaps(Lpos1, Rpos1, Lpos2, Rpos2);}
// 		if (Lpos2 > Rpos1) {Lpos1 = Rpos1; tmp1 = tmp1->next_segment;}
// 		else if (Lpos1 > Rpos2) {Lpos2 = Rpos2; tmp2 = tmp2->next_segment;}
// 		else {
// 			if (Rpos2 > Rpos1) {Lpos1 = Rpos1; tmp1 = tmp1->next_segment;}
// 			else {Lpos2 = Rpos2; tmp2 = tmp2->next_segment;}
// 		}
// 	} 
// 	return total_IBD;
// };

// static void convert_hap(char* marker_vec, haplotype* hap, const genotype_map& founder_geno_map, const std::vector<int>& map ){
// 	int n_markers = map.size(); //****	
// 	haplotype* tmp = hap;
// 	const char* founder_hap = founder_geno_map.find(tmp->hap)->second.data();

// 	int marker_pos;
// 	for (int i=0; i < n_markers; i++){
// 		marker_pos = map[i];
// 		if (marker_pos > tmp->pos){
// 			tmp = tmp->next_segment;
// 			founder_hap = founder_geno_map.find(tmp->hap)->second.data();
// 		}
// 		marker_vec[i] = founder_hap[i];
// 	}
// };

// static void no_convert(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array){
// 	for(int k=0; k<nbrecomb; k++){
// 		BP_array[k] = (int) (CO_array[k]*bp_len);
// 		if (k>0){
// 			if(BP_array[k]==BP_array[k-1]){
// 				BP_array[k]=BP_array[k]+1;
// 			}
// 		}
// 	}
// }

// static void convert1(int& nbrecomb, double* CO_array, const double& Morgan_len, const int& bp_len, int* bp_map, double* cm_map, int* BP_array){
// 	for (int k=0; k<nbrecomb; k++){
// 		double physical_distance;
// 		double genetic_distance = CO_array[k]*Morgan_len*100;

// 		int map_index = 0;
// 		while(genetic_distance > cm_map[map_index]) map_index++;

// 		physical_distance = bp_map[map_index-1]  + (bp_map[map_index] - bp_map[map_index-1])*(genetic_distance - cm_map[map_index-1])/(cm_map[map_index] - cm_map[map_index-1]);
// 		BP_array[k] = (int)(physical_distance);

// 		if (k>0){
// 			if(BP_array[k]==BP_array[k-1]){
// 				BP_array[k]=BP_array[k]+1;
// 			}
// 		}
// 	}
// }

// //F for father, not female
// static void makeRecombF( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle, const int& BP_len)
// {
//     haplotype *perehap1, *perehap2;
// 	if(Ordre_tmp->pere != NULL){
// 		if (nbRecomb > 0){
// 			if (probHap < 0.5){
// 				perehap1=(*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
// 				perehap2=(*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
// 			} 
// 			else{
// 				perehap1=(*hapRef).find(Ordre_tmp->pere->clesHaplo_2)->second;
// 				perehap2=(*hapRef).find(Ordre_tmp->pere->clesHaplo_1)->second;
// 			}

// 			haplotype *hapChild_1 = new haplotype();
// 			haplotype *hapChild_deb1 = hapChild_1;
// 			recombine(perehap1, perehap2, hapChild_deb1, nbRecomb, posRecomb, BP_len);
// 			Ordre_tmp->clesHaplo_1 = cle;
// 			(*hapRef)[cle++] = hapChild_1;
// 		}
// 		else{ //If no recombination just pass one of father's chromosomes down 
// 			if(probHap<0.50){
// 				Ordre_tmp->clesHaplo_1=Ordre_tmp->pere->clesHaplo_1;
// 			}
// 			else{
// 				Ordre_tmp->clesHaplo_1=Ordre_tmp->pere->clesHaplo_2;
// 			}
// 		}
// 	}
// 	else Ordre_tmp -> clesHaplo_1 = 0;
// }

// //M for mother, not male
// static void makeRecombM( CIndSimul *Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, int *posRecomb, int &cle, const int& BP_len)
// {
//     haplotype *merehap1, *merehap2;
//  	if(Ordre_tmp->mere != NULL){
// 		if (nbRecomb > 0){
// 			if (probHap < 0.5){
// 				merehap1=(*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
// 				merehap2=(*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
// 			} 
// 			else{
// 				merehap1=(*hapRef).find(Ordre_tmp->mere->clesHaplo_2)->second;
// 				merehap2=(*hapRef).find(Ordre_tmp->mere->clesHaplo_1)->second;
// 			}

// 			haplotype *hapChild_2 = new haplotype();
// 			haplotype *hapChild_deb2 = hapChild_2;
// 			recombine(merehap1, merehap2, hapChild_deb2, nbRecomb, posRecomb, BP_len);
// 			Ordre_tmp->clesHaplo_2 = cle;
// 			(*hapRef)[cle++] = hapChild_2;
// 		}
// 		else{ //If no recombination just pass one of mother's chromosomes down 
// 			if(probHap<0.50){
// 				Ordre_tmp->clesHaplo_2=Ordre_tmp->mere->clesHaplo_1;
// 			}
// 			else{
// 				Ordre_tmp->clesHaplo_2=Ordre_tmp->mere->clesHaplo_2;
// 			}
// 		}
// 	}
// 	else Ordre_tmp -> clesHaplo_2 = 0;
// }

// static void recombine(haplotype* hapBegin, haplotype* hapEnd, haplotype* hapChild, int nbRecomb, int* posRecomb, const int& BP_len)
// {
// 	haplotype* hap_active = hapBegin;

// 	for (int i=0; i < nbRecomb; i++){
// 		int position = posRecomb[i];
// 		// de 0 a posRecomb on prend hapBegin
// 		while(position > hap_active->pos && hap_active->pos != BP_len) {
// 			(*hapChild).hap          = hap_active->hap;
// 			(*hapChild).pos          = hap_active->pos;
// 			(*hapChild).fixe         = 0;
// 			(*hapChild).next_segment = new haplotype();//[1];
// 			hapChild                 = hapChild->next_segment;
// 			hap_active               = hap_active->next_segment;
// 		}

// 		// on ajoute la recomb pour hapChild
// 		(*hapChild).hap          = hap_active->hap;
// 		(*hapChild).pos          = position;
// 		(*hapChild).fixe         = 0;

// 		if(i%2 == 0){
// 			hap_active = hapEnd;
// 		}
// 		else hap_active = hapBegin;

// 		// on met le pointeur de hapEnd a la bonne place.
// 		while(position >= hap_active->pos && hap_active->pos != BP_len) hap_active = hap_active->next_segment;

// 		// on verifie que l'haplotype qui suit n'est pas le meme. Si oui, on met la nouvelle position.
// 		if(hap_active->hap == hapChild->hap){
// 			(*hapChild).pos          = hap_active->pos;
// 		}
// 		else{
// 			(*hapChild).next_segment = new haplotype();//[1];
// 			hapChild                 = hapChild->next_segment;
// 			(*hapChild).hap          = hap_active->hap;
// 			(*hapChild).pos          = hap_active->pos;
// 			(*hapChild).fixe         = 0;
// 		}
// 	}

// 	while(hap_active->pos != BP_len) {
// 		hap_active               = hap_active->next_segment;
// 		(*hapChild).next_segment = new haplotype();//[1]; 
// 		hapChild                 = hapChild->next_segment;
// 		(*hapChild).hap          = hap_active->hap;
// 		(*hapChild).pos          = hap_active->pos;
// 		(*hapChild).fixe         = 0;
// 	}
// };

// void read_binary(genotype_map& founder_genotypes, std::ofstream& gen_bfile)
// {
// 	int proID, buffer, read_int;
// 	std::vector<int> segpos_vec, segnom_vec;
	
// 	gen_bfile.read(&proID, 4);
// 	gen_bfile.read(&buffer, 4);

// 	for (int i=0; i < buffer; i++){
// 		gen_bfile.read(&read_int, 4);
// 		segpos_vec.push_back(read_int);
// 		gen_bfile.read(&read_int, 4);
// 		segnom_vec.push_back(read_int);
// 	}
// 	gen_bfile.read(&read_int, 4);
		
// }

// void pIBD_matrix(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
// 				double* probRecomb, double* Morgan_Len, int BP_len, int model,  
// 				int convert, double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, double* R_matrix, 
// 				const std::string& out, const std::string& map_filepath, const std::string& ped_filepath, int seed) 
// {
// 	try{

// 	//CREATION DE TABLEAU D'INDIVIDU
// 	int lNIndividu;
// 	CIndSimul *Noeud=NULL;
// 	LoadGenealogie(Genealogie, GTRUE, &lNIndividu, &Noeud);

// 	//CREATION D'UN VECTEUR DE PROPOSANT
// 	CIndSimul **NoeudPro=NULL;
// 	LoadProposant(plProposant,lNProposant,&NoeudPro);

// 	//CREATION OF AN ANCESTOR VECTOR *** with the reference haplotypes for the chosen ancestors. ***
// 	CIndSimul **NoeudAnc=NULL;
// 	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

// 	//Creation des tableau
// 	INITGESTIONMEMOIRE;
// 	CIndSimul** Ordre = (CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));

// 	//Pour le sort spcial		
// 	int*	OrdreSaut	= (int*) memalloc(lNIndividu,sizeof(int*));				
// 	int NOrdre;
	
// 	int i;

// 	//Initialize all the nodes
// 	for(i=0;i<lNIndividu;i++)
// 	{
// 		Noeud[i].allele = 0;
// 		Noeud[i].etat=GENNONEXPLORER;
// 		Noeud[i].bFlagSort=0;
// 		Noeud[i].clesHaplo_1 = 0; //"0.1";
// 		Noeud[i].clesHaplo_2 = 0; //"0.1";
// 	}

// 	std::unordered_map<int, haplotype*> hapRef; // empty unordered_map
// 	hapRef.reserve(lNIndividu);
// 	haplotype *hapVide = new haplotype();
// 	hapVide->hap  = 0;
// 	hapVide->pos  = BP_len;
// 	hapVide->fixe = 1;
// 	hapRef[0]=hapVide;

// 	//label the nodes that are probands
// 	for(i=0;i<lNProposant;i++){
// 		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
// 	}

// 	int cleFixe = 1; // haplotype keys 
// 	//identifier et etiqueter les points de departs et les haplos ancetres (starting points and ancestor haplotypes)
// 	for(i=0;i<lNAncetre;i++)
// 	{		
// 		NoeudAnc[i]->allele = 0;
// 		NoeudAnc[i]->etat=GENDEPART;
// 	    NoeudAnc[i]->clesHaplo_1 = cleFixe++;
// 	    NoeudAnc[i]->clesHaplo_2 = cleFixe++;

// 		haplotype *tmp1 = new haplotype();//[1];
// 		tmp1->hap  = NoeudAnc[i]->nom;
// 		tmp1->pos  = BP_len;
// 		tmp1->fixe = 1;
// 		hapRef[NoeudAnc[i]->clesHaplo_1] = tmp1;
		
// 		haplotype *tmp2 = new haplotype();//[1];
// 		tmp2->hap  = -1 * NoeudAnc[i]->nom;
// 		tmp2->pos  = BP_len;
// 		tmp2->fixe = 1;
// 		hapRef[NoeudAnc[i]->clesHaplo_2] = tmp2;
// 	}

// 	//identifier et marque les noeuds utile et ceux inutile a la recherche
// 	for(i=0;i<lNAncetre;i++) ExploreArbre(NoeudAnc[i]);
		
// 	//create the order of traversal and calculate the jumps (Jumps are unecessary, only for speeding up allele calculations, will test if it works without)
// 	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
// 	NOrdre=0;

// 	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
// 	for(i=0;i<lNAncetre;i++)StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut); // les infos de NoeudAnc sont pointes par Ordre dans "le bon ordre".

// 	//create mersenne twister generator to use for random  distributions 
// 	std::mt19937 my_rng = std::mt19937(seed);
// 	//initialize crossover class
// 	Crossovers crossovers;
  
// 	void (Crossovers::*SampleCO)(const int&, double*, double*, int&, std::mt19937&, double*);
	
// 	if 		(model==1) SampleCO=&Crossovers::Poisson_CO;
// 	else if (model==2) SampleCO=&Crossovers::Poisson_ZT;
// 	else if (model==3){
// 		SampleCO=&Crossovers::Gamma_CO;
// 		crossovers.init_gamma(probRecomb[0], probRecomb[1], Morgan_Len[0], Morgan_Len[1]);
// 	}

// 	void (*convert_dist)(int&, double*, const double&, const int&, int*, double*, int*); //(nbrecomb, CO_array, Morgan_len, bp_len, bp_map, cm_map, precision)

// 	if 		(convert == 0) convert_dist = &no_convert;  // no genetic/physical map is specified, then no need to scale wrt physical distance (assumed 1:1)
// 	else if (convert == 1) convert_dist = &convert1; 	// genetic/physical map is specified
	
// 	double pHap1, pHap2;
// 	double CO_arrayF[20]; //hold crossover positions for Father's chromosome
// 	double CO_arrayM[20]; //hold crossover positions for Mother's chromsome

// 	int BP_CO_arrayF[20]; //crossover positions in BP
// 	int BP_CO_arrayM[20];

// 	int nbRecomb1 =0;
// 	int nbRecomb2 =0;
	
// 	std::uniform_real_distribution<> u_dist(0, 1);

// 	//run simulation
// 	int clesSim = cleFixe;
// 	for(int i=0;i<NOrdre;i++) {
// 		//Simulate meiosis in the parents, store the location of crossovers (in genetic distance scaled to [0,1]) in CO_array.
// 		(crossovers.*SampleCO)(1, probRecomb, Morgan_Len, nbRecomb1, my_rng, CO_arrayF);
// 		(crossovers.*SampleCO)(2, probRecomb, Morgan_Len, nbRecomb2, my_rng, CO_arrayM);

// 		// locations of crossovers in CO_array will be converted to physical distance. 
// 		convert_dist(nbRecomb1, CO_arrayF, Morgan_Len[0], BP_len, bp_map_FA ,cm_map_FA, BP_CO_arrayF);
// 		convert_dist(nbRecomb2, CO_arrayM, Morgan_Len[1], BP_len, bp_map_MO ,cm_map_MO, BP_CO_arrayM);

// 		pHap1 = u_dist(my_rng); 
// 		makeRecombF(Ordre[i], &hapRef, pHap1, nbRecomb1, BP_CO_arrayF, clesSim, BP_len);

// 		pHap2 = u_dist(my_rng);
// 		makeRecombM(Ordre[i], &hapRef, pHap2, nbRecomb2, BP_CO_arrayM, clesSim, BP_len);
// 	}
// 	//calculate prop. IBD matrix from simulated proband haplotypes
// 	//write to the matrix in R
// 	haplotype *h1_1, *h1_2, *h2_1, *h2_2;
// 	int total_IBD;
// 	for(int i=0;i<lNProposant;i++){
// 		h1_1 = hapRef.find(NoeudPro[i]->clesHaplo_1)->second;
// 		h1_2 = hapRef.find(NoeudPro[i]->clesHaplo_2)->second;

// 		int j = 0;
// 		for( ; j < i; j++){
// 			total_IBD = 0;
// 			h2_1 = hapRef.find(NoeudPro[j]->clesHaplo_1)->second;
// 			h2_2 = hapRef.find(NoeudPro[j]->clesHaplo_2)->second;
// 			total_IBD = total_IBD + check_p_IBD(h1_1, h2_1, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_1, h2_2, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_2, h2_1, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_2, h2_2, BP_len);
// 			R_matrix[i * lNProposant + j] = total_IBD/4;
// 		}
// 		for( ; j < lNProposant; j++){
// 			R_matrix[i * lNProposant + j] = 0;
// 		}
// 	}

// 	//read in the map file, get positions and # markers
// 	std::vector<int> map_pos;
// 	read_map_file(map_filepath, map_pos);
// 	int n_markers = map_pos.size();

// 	//read in founder haplotypes into hash map
// 	// +indID is key for first (paternal) copy of ind, -indID is key for maternal copy
// 	genotype_map founder_genotypes;
// 	founder_genotypes.reserve(2 * lNProposant + 1);
// 	//read ped file
// 	read_ped_file(ped_filepath, founder_genotypes, n_markers);

// 	//create output file
// 	std::ofstream outfile(out + ".ped");
// 	std::vector<char> hap1(n_markers), hap2(n_markers);
// 	for(int i=0; i<lNProposant; i++){
// 		//put the first 6 col of pedfile
// 		outfile << "0\t" << NoeudPro[i]->nom << "\t0\t0\t-9\t-9\t";
// 		for (int j = 0; j<n_markers; j++){
// 			outfile << hap1[j] << '\t' << hap2[j] << '\t';
// 		}
// 		outfile << std::endl;
// 	}	
// 	outfile.close();


// 	//delete haplotypes 
// 	//can delete the internal haplotypes earlier to free up some memory first
// 	for(int i=0; i<clesSim; i++) {
// 		haplotype* tmp = hapRef.find(i)->second;//hapKey.second;
// 		while(tmp->next_segment != NULL) {
// 			haplotype* tmp_back = tmp;
// 			tmp = tmp->next_segment;
// 			delete tmp_back;
// 		}
// 		delete tmp;
// 	}

// 	} catch(std::exception &ex) {
// 	forward_exception_to_r(ex);
// 	} catch(...){
// 	::Rf_error("c++ exception (unknown reason)"); 
// 	} 
// };

// void pIBD_matrix2(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
// 				double* probRecomb, double* Morgan_Len, int BP_len, int model,  
// 				int convert, double* cm_map_FA, double* cm_map_MO, int* bp_map_FA, int* bp_map_MO, double* R_matrix, std::string out, int seed) 
// {
// 	try{
// 	//CREATION DE TABLEAU D'INDIVIDU
// 	int lNIndividu;
// 	CIndSimul *Noeud=NULL;
// 	LoadGenealogie(Genealogie, GTRUE, &lNIndividu, &Noeud);

// 	//CREATION D'UN VECTEUR DE PROPOSANT
// 	CIndSimul **NoeudPro=NULL;
// 	LoadProposant(plProposant,lNProposant,&NoeudPro);

// 	//CREATION OF AN ANCESTOR VECTOR *** with the reference haplotypes for the chosen ancestors. ***
// 	CIndSimul **NoeudAnc=NULL;
// 	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

// 	//Creation des tableau
// 	INITGESTIONMEMOIRE;
// 	CIndSimul** Ordre = (CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));

// 	//Pour le sort spcial		
// 	int*	OrdreSaut	= (int*) memalloc(lNIndividu,sizeof(int*));				
// 	int NOrdre;
	
// 	int i;

// 	//Initialize all the nodes
// 	for(i=0;i<lNIndividu;i++)
// 	{
// 		Noeud[i].allele = 0;
// 		Noeud[i].etat=GENNONEXPLORER;
// 		Noeud[i].bFlagSort=0;
// 		Noeud[i].clesHaplo_1 = 0; //"0.1";
// 		Noeud[i].clesHaplo_2 = 0; //"0.1";
// 	}

// 	std::unordered_map<int, haplotype*> hapRef; // empty unordered_map
// 	hapRef.reserve(lNIndividu);
// 	haplotype *hapVide = new haplotype();
// 	hapVide->hap  = 0;
// 	hapVide->pos  = BP_len;
// 	hapVide->fixe = 1;
// 	hapRef[0]=hapVide;

// 	//label the nodes that are probands
// 	for(i=0;i<lNProposant;i++){
// 		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
// 	}

// 	int cleFixe = 1; // haplotype keys 
// 	//identifier et etiqueter les points de departs et les haplos ancetres (starting points and ancestor haplotypes)
// 	for(i=0;i<lNAncetre;i++)
// 	{		
// 		NoeudAnc[i]->allele = 0;
// 		NoeudAnc[i]->etat=GENDEPART;
// 	    NoeudAnc[i]->clesHaplo_1 = cleFixe++;
// 	    NoeudAnc[i]->clesHaplo_2 = cleFixe++;

// 		haplotype *tmp1 = new haplotype();//[1];
// 		tmp1->hap  = NoeudAnc[i]->nom;
// 		tmp1->pos  = BP_len;
// 		tmp1->fixe = 1;
// 		hapRef[NoeudAnc[i]->clesHaplo_1] = tmp1;
		
// 		haplotype *tmp2 = new haplotype();//[1];
// 		tmp2->hap  = -1 * NoeudAnc[i]->nom;
// 		tmp2->pos  = BP_len;
// 		tmp2->fixe = 1;
// 		hapRef[NoeudAnc[i]->clesHaplo_2] = tmp2;
// 	}

// 	//identifier et marque les noeuds utile et ceux inutile a la recherche
// 	for(i=0;i<lNAncetre;i++) ExploreArbre(NoeudAnc[i]);
		
// 	//create the order of traversal and calculate the jumps (Jumps are unecessary, only for speeding up allele calculations, will test if it works without)
// 	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
// 	NOrdre=0;

// 	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
// 	for(i=0;i<lNAncetre;i++)StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut); // les infos de NoeudAnc sont pointes par Ordre dans "le bon ordre".

// 	//create mersenne twister generator to use for random  distributions 
// 	std::mt19937 my_rng = std::mt19937(seed);
// 	//initialize crossover class
// 	Crossovers crossovers;
  
// 	void (Crossovers::*SampleCO)(const int&, double*, double*, int&, std::mt19937&, double*);
	
// 	if 		(model==1) SampleCO=&Crossovers::Poisson_CO;
// 	else if (model==2) SampleCO=&Crossovers::Poisson_ZT;
// 	else if (model==3){
// 		SampleCO=&Crossovers::Gamma_CO;
// 		crossovers.init_gamma(probRecomb[0], probRecomb[1], Morgan_Len[0], Morgan_Len[1]);
// 	}

// 	void (*convert_dist)(int&, double*, const double&, const int&, int*, double*, int*); //(nbrecomb, CO_array, Morgan_len, bp_len, bp_map, cm_map, precision)

// 	if 		(convert == 0) convert_dist = &no_convert;  // no genetic/physical map is specified, then no need to scale wrt physical distance (assumed 1:1)
// 	else if (convert == 1) convert_dist = &convert1; 	// genetic/physical map is specified
	
// 	double pHap1, pHap2;
// 	double CO_arrayF[20]; //hold crossover positions for Father's chromosome
// 	double CO_arrayM[20]; //hold crossover positions for Mother's chromsome

// 	int BP_CO_arrayF[20]; //crossover positions in BP
// 	int BP_CO_arrayM[20];

// 	int nbRecomb1 =0;
// 	int nbRecomb2 =0;
	
// 	std::uniform_real_distribution<> u_dist(0, 1);

// 	//run simulation
// 	int clesSim = cleFixe;
// 	for(int i=0;i<NOrdre;i++) {
// 		//Simulate meiosis in the parents, store the location of crossovers (in genetic distance scaled to [0,1]) in CO_array.
// 		(crossovers.*SampleCO)(1, probRecomb, Morgan_Len, nbRecomb1, my_rng, CO_arrayF);
// 		(crossovers.*SampleCO)(2, probRecomb, Morgan_Len, nbRecomb2, my_rng, CO_arrayM);

// 		// locations of crossovers in CO_array will be converted to physical distance. 
// 		convert_dist(nbRecomb1, CO_arrayF, Morgan_Len[0], BP_len, bp_map_FA ,cm_map_FA, BP_CO_arrayF);
// 		convert_dist(nbRecomb2, CO_arrayM, Morgan_Len[1], BP_len, bp_map_MO ,cm_map_MO, BP_CO_arrayM);

// 		pHap1 = u_dist(my_rng); 
// 		makeRecombF(Ordre[i], &hapRef, pHap1, nbRecomb1, BP_CO_arrayF, clesSim, BP_len);

// 		pHap2 = u_dist(my_rng);
// 		makeRecombM(Ordre[i], &hapRef, pHap2, nbRecomb2, BP_CO_arrayM, clesSim, BP_len);
// 	}
// 	//calculate prop. IBD matrix from simulated proband haplotypes
// 	//write to the matrix in R
// 	haplotype *h1_1, *h1_2, *h2_1, *h2_2;
// 	int total_IBD;
// 	for(int i=0;i<lNProposant;i++){
// 		h1_1 = hapRef.find(NoeudPro[i]->clesHaplo_1)->second;
// 		h1_2 = hapRef.find(NoeudPro[i]->clesHaplo_2)->second;

// 		int j = 0;
// 		for( ; j < i; j++){
// 			total_IBD = 0;
// 			h2_1 = hapRef.find(NoeudPro[j]->clesHaplo_1)->second;
// 			h2_2 = hapRef.find(NoeudPro[j]->clesHaplo_2)->second;
// 			total_IBD = total_IBD + check_p_IBD(h1_1, h2_1, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_1, h2_2, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_2, h2_1, BP_len);
// 			total_IBD = total_IBD + check_p_IBD(h1_2, h2_2, BP_len);
// 			R_matrix[i * lNProposant + j] = total_IBD/4;
// 		}
// 		for( ; j < lNProposant; j++){
// 			R_matrix[i * lNProposant + j] = 0;
// 		}
// 	}

// 	std::ofstream outfile("qqq.GENLIB", std::ofstream::binary);
//     outfile.write((char*) &numSim, sizeof(int));
// 	outfile.write((char*) &numPro, sizeof(int));

// 	for(int i=0; i<lNProposant; i++){
// 		write_binary(hapRef.find(NoeudPro[i]->clesHaplo_1)->second, outfile, NoeudPro[i]->nom);
// 		write_binary(hapRef.find(NoeudPro[i]->clesHaplo_2)->second, outfile, -1 * NoeudPro[i]->nom);
// 	}

// 	//delete haplotypes 
// 	//can delete the internal haplotypes earlier to free up some memory first
// 	for(int i=0; i<clesSim; i++) {
// 		haplotype* tmp = hapRef.find(i)->second;//hapKey.second;
// 		while(tmp->next_segment != NULL) {
// 			haplotype* tmp_back = tmp;
// 			tmp = tmp->next_segment;
// 			delete tmp_back;
// 		}
// 		delete tmp;
// 	}

// 	} catch(std::exception &ex) {
// 	forward_exception_to_r(ex);
// 	} catch(...){
// 	::Rf_error("c++ exception (unknown reason)"); 
// 	} 
// };

// void write_binary(haplotype *tmp, std::ofstream& outfile, int name){
// 	bool done = false;
// 	int filepos1, filepos2, blocklen, counter;
// 	outfile.write ((char*)& name, sizeof(int));
// 	filepos1 = outfile.tellp();
// 	outfile.write ((char*)& blocklen, sizeof(int));
// 	while(not done){
// 		counter = 1;
// 		outfile.write((char*) & (tmp->pos), sizeof(int));
// 		outfile.write((char*) & (tmp->hap), sizeof(int));
// 		if (tmp->next_segment == NULL){
// 			done = true;
// 			filepos2 = outfile.tellp();
// 			outfile.seekp(filepos1);
// 			outfile.write( (char*) &counter, sizeof(int));
// 			outfile.seekp(filepos2);
// 		}
// 		else {
// 			tmp = tmp->next_segment;
// 			counter++;
// 		}
// 	}
// };


// //have the binary proband info in terms of segments and positions
// //have founder info file... 
// //is best way to still convert founder info to unordered map of vectors
// void convert_binary_pedmap(std::ifstream& bfile, std::string& map_filepath, std::string& ped_filepath, std::vector<int>& pro_vec, std::vector<int>& bool shuffle){

// 	//read npro from bfile
// 	uint32_t nSim, nPro, nAnc;
// 	bfile.read(reinterpret_cast<char *>(&nSim), sizeof(nSim));
// 	bfile.read(reinterpret_cast<char *>(&nPro), sizeof(nPro));
// 	bfile.read(reinterpret_cast<char *>(&nAnc), sizeof(nAnc));
	
// 	// assert nPro == pro_vec.size()

// 	std::vector<int> map_pos;
// 	read_map_file(map_filepath, map_pos);
// 	int n_markers = map_pos.size();

// 	//hoqw do we get the number of founders?
// 	genotype_map founder_genotypes;
// 	founder_genotypes.reserve(2 * nAnc + 1);
// 	//only keep the founder in founder vector
// 	//return message of additional founders in the file
// 	read_ped_file(ped_filepath, founder_genotypes, n_markers);

// 	for (int i=0; i < nSim; i++){

// 	}


// 	std::vector<char> hap1(n_markers), hap2(n_markers);
// 	//throw error if an IID ius not in pro_vec
// 	for (int i=0; i < nPro; i++){
// 		read_binary_line(bfile, hap1, map_pos, founder_genotypes);
// 		read_binary_line(bfile, hap2, map_pos, founder_genotypes);
// 		write_ped_line(pedfile, hap1, hap2, ID);
// 	}	
// };

// void write_ped_line(std::ofstream& pedfile, std::vector<char>& hap1, std::vector<char>& hap2, int& IID){
// 	int n_markers = hap1.size();
// 	//assert hap1 size = hap2 size
// 	//put the first 6 col of pedfile
// 	pedfile << "0\t" << IID << "\t0\t0\t-9\t-9\t";
// 	for (int j = 0; j<n_markers; j++){
// 			pedfile << hap1[j] << '\t' << hap2[j] << '\t';
// 	}
// 	pedfile << "\n";
// }

// //for the binary gene drop output
// //if pass the ifstream by ref it will keep its buffer location
// void read_binary_line(std::ifstream& bfile, std::vector<char>& geno_vec, std::vector<int>& map, const genotype_map& founder_genotypes){
// 	int IID, blocklen, segID, segPos;
// 	bfile.read(reinterpret_cast<char *>(&IID), 4);
// 	bfile.read(reinterpret_cast<char *>(&blocklen), 4);	

// 	int marker_count = 0;
// 	for (int i=0; i < blocklen; i++){
// 		bfile.read(reinterpret_cast<char *>(&segID), 4);
// 		bfile.read(reinterpret_cast<char *>(&segPos), 4);			
// 		while (map[marker_count] <= segPos){
// 			geno_vec[marker_count] = founder_genotypes.find(segID)->second.data()[marker_count] ;
// 			marker_count++;
// 		}
// 	}
// }