#include <sdsl/suffix_arrays.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>

//#define K 5

using namespace sdsl;

int main(int argc, char** argv) {
  // Default read file is 0, can be changed with first command line argument
  std::string read_num = "0";
  // Default edit distance is 10, can change with second command line argument
  int K = 10;
  if(argc>1) read_num = argv[1];
  if(argc>2) K = atoi(argv[2]);
  std::string in_file = "reads_"+read_num;
  std::string out_file = "assembly_"+read_num;
  
  csa_wt<> fm_index;
  std::ifstream genome_in("ref_genome");
  std::string genome; 
  genome_in >> genome;
  construct_im(fm_index, genome.c_str(), 1);
  
  std::ifstream reads_in(in_file);
  std::string read;

  std::ofstream assembly_out(out_file);
  std::string result;
  result = genome;
  
  std::vector<int> A(1000000,0);
  std::vector<int> T(1000000,0);
  std::vector<int> C(1000000,0);
  std::vector<int> G(1000000,0);

  int no_match_count = 0;
  int multiple_matches = 0;
  
  while(reads_in >> read){
    std::set<uint64_t> final_loc; 
    sdsl::int_vector<64> locations = locate(fm_index, read.begin(), read.end());
    if(locations.size() == 0) {
      int seg_len = round(read.size()/((double)(K+1)));
      for(int i=0;i<(K+1);i++){
	      int start = i*seg_len;
	      int end = std::min((i+1)*seg_len, (int)read.length());
	      std::string sub_str = read.substr(start, end-start);
	      locations = locate(fm_index, sub_str.begin(), sub_str.end());
	      for(uint64_t loc : locations){
	        if((int)loc-start<0 || ((int)loc-start+read.length()) > genome.length()){
	          continue;
	        }
	        int mismatches = 0;
	        for(int j=0; j<start && mismatches<=K; j++){
	          if (read.at(j) != genome.at(loc - start + j)) mismatches++;
	        }
	        for(int j=end; j<read.length() && mismatches<=K; j++){
	          if (read.at(j) != genome.at(loc - start + j)) mismatches++;
	        }
	        if(mismatches<=K) final_loc.insert(loc-start); 
	      }
      }
    }
    else {
      for(uint64_t loc : locations) {
	      final_loc.insert(loc);
      }
    }
    if(final_loc.size()==0) no_match_count++;
    if(final_loc.size()>1) multiple_matches += final_loc.size()-1;
    for(uint64_t loc: final_loc){
      for(int i=0; i<read.size(); i++){
	      if(read.at(i)=='A'){
	        A[loc+i]++;
	      }
	      else if(read.at(i)=='T'){
	        T[loc+i]++;
	      }
	      else if(read.at(i)=='C'){
	        C[loc+i]++;
        }
	      else if(read.at(i)=='G'){
	        G[loc+i]++;
	      }
      }
    }
  }
  for(int i=0; i<result.size(); i++){
    if(A[i]>=T[i]&&A[i]>=C[i]&&A[i]>=G[i]){
      result[i] = 'A';
    }
    else if(T[i]>=A[i]&&T[i]>=C[i]&&T[i]>=G[i]){
      result[i] = 'T';
    }
    else if(C[i]>=A[i]&&C[i]>=T[i]&&C[i]>=G[i]){
      result[i] = 'C';
    }
    else if(G[i]>=A[i]&&G[i]>=T[i]&&G[i]>=C[i]){
      result[i] = 'G';
    }
  }
  assembly_out << result;
  std::cout << "Assembly written to " << out_file << "\n";
  std::cout << "Reads not matched: " << no_match_count << "\n";
  std::cout << "Number of excess matches: " << multiple_matches << "\n";
  genome_in.close();
  reads_in.close();
  assembly_out.close();
  return 0;
}
