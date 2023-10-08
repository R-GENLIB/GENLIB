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