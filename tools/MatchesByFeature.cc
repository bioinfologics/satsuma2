#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iomanip>
#include <string>

class Chromosome{
  public:
    std::string name;
    uint64_t start,end;
    std::vector<uint8_t> featuresv;
    Chromosome(std::string _name, uint64_t _start, uint64_t _end){
      name=_name;
      start=_start;
      end=_end;
      featuresv=std::vector<uint8_t>(end-start+1,0);
    }
    ~Chromosome(){
      featuresv.clear();
    }
    Chromosome operator=(const Chromosome & rhs){
      name=rhs.name;
      start=rhs.start;
      end=rhs.end;
      featuresv=rhs.featuresv;
    }
};

class MatchesByFeatureTracker{
  public:
    MatchesByFeatureTracker(char * _gff3_filename, std::vector<std::string> _features, std::vector<std::string> _matches_filenames);
    void process_all_files(int np);// uses async to process multiple files in one go   
    void print_results();
  private:
    std::vector<uint64_t> process_file(int file_id);//processes a single file, returns results in the results vector, deep copies the chromosomes
    std::vector<std::string> features;
    std::vector<std::string> matches_filenames;
    std::vector<std::vector<uint64_t>> results;
    std::vector<uint64_t> feature_totals;
    std::map<std::string,Chromosome> chromosomes;
};

MatchesByFeatureTracker::MatchesByFeatureTracker(char * _gff3_filename, std::vector<std::string> _features, std::vector<std::string> _matches_filenames){
  std::ifstream infile;
  std::string line,sline[5];
  features=_features;
  matches_filenames=_matches_filenames;
  //TODO: copies features
  //TODO: initializes Chromosomes
  //takes a gff file, loads a uint8 with bit-field values for presence for up to 8 features being present or not per position.
  std::cout<<"Creating chromosome structure... "<<std::flush;
  infile.open(_gff3_filename);
  //for each chromosome line
  while (getline (infile,line) ) {
    if ('#'==line[0]) continue;
    std::stringstream ss(line);
    for (auto i=0; i<5; i++){
    //create an entry on the chromosome list with name, start, size
      ss>>sline[i];
    //TODO: filter small chromosomes or whatever? can be done beforehand!
    }
    if (sline[2]=="chromosome" && sline[0].size()<4) {
      //std::cout<<sline[0]<<" "<<sline[3]<<"-"<<sline[4]<<std::endl;
      Chromosome c= Chromosome(sline[0],atol(sline[3].c_str()),atol(sline[4].c_str()));
      chromosomes.insert(std::make_pair(sline[0], c));
    }
  }
  infile.close();
  std::cout<<"DONE!"<<std::endl;
  
  infile.open(_gff3_filename);
  std::cout<<"Loading features... "<<std::flush;
  while (getline (infile,line)) {
    if ('#'==line[0]) continue;
    std::stringstream ss(line);
    for (auto i=0; i<5; i++){
      ss>>sline[i];
    }
    std::map<std::string,Chromosome>::iterator c=chromosomes.find(sline[0]);
    if ( chromosomes.end() != c ){//XXX: this "eliminates" hanging-in-the-air features
      for (auto i=0;i<features.size();i++) {
        if (features[i]==sline[2]) {
          for (auto j=atol(sline[3].c_str());j <= atol(sline[4].c_str()); j++)
            c->second.featuresv[j-c->second.start]|=(1<<i);
          break;
        }
      }
    }
  }
  std::cout<<"DONE!"<<std::endl;
  infile.close();

  std::cout<<"computing feature totals... "<<std::flush;
  uint64_t t[128];
  for (auto i=0;i<128;i++) t[i]=0;
  for (auto & c : chromosomes){
    for (auto & r: c.second.featuresv) t[r]++;
  }
  feature_totals.resize(features.size()+1);
  for (auto i=0;i<features.size();i++){
    feature_totals[i]=0;
    for (auto j=0;j<128;j++){
      if ((1<<i)&j) {
        if (128&j) feature_totals[i]+=t[j];
      }
    }
  }
  feature_totals[features.size()]=t[0];
  std::cout<<"DONE!"<<std::endl;
  //TODO:initialize results
  results.resize(matches_filenames.size());
}

std::vector<uint64_t> MatchesByFeatureTracker::process_file(int _id){
  //TODO: DEEP COPY of the chromosomes into chrom
  std::map<std::string,Chromosome> chrom=chromosomes;
  std::cout<<"Loading alignement coords... "<<std::flush;
  std::ifstream infile;
  std::string line,sline[5];
  infile.open( matches_filenames[_id].c_str());
  while (getline (infile,line)){
    std::stringstream ss(line);
    for (auto i=0; i<3; i++) ss>>sline[i];
    std::map<std::string, Chromosome>::iterator c=chrom.find(sline[0]);
    if ( chrom.end() != c ){//XXX: this "eliminates" hanging-in-the-air features
      for (auto j=atol(sline[1].c_str());j < atol(sline[2].c_str()); j++)///TODO: check it's 0-based in-between
        c->second.featuresv[j]|=128;
    }
  }
  std::cout<<"DONE!"<<std::endl;
  infile.close();
  //prints a stat for covered% for each of the features across the genome
  uint64_t h[256];
  //std::cout<<"Counting bp coverage... "<<std::flush;
  for (auto i=0;i<256;i++) h[i]=0;//better safe than sorry
  for (auto & c : chrom){
    uint8_t * bc_data=c.second.featuresv.data();
    for (unsigned long i=0; i<c.second.featuresv.size(); i++){
      h[bc_data[i]]++;
    }
    //std::cout<<"+"<<std::flush;
  }
  //std::cout<<"DONE!"<<std::endl;
  //for (auto i=0;i<256;i++) std::cout <<i<<": "<<h[i]<<std::endl;
  //std::cout<<"Breakdown by feature:"<<std::endl;
  std::vector<uint64_t> covered;
  covered.resize(features.size()+1);
  for (auto i=0;i<features.size();i++){
    covered[i]=0;
    for (auto j=0;j<256;j++){
      if ((1<<i)&j) {
        if (128&j) covered[i]+=h[j];
      }
    }
  }
  covered[features.size()]=h[128];
  return covered;
}
void MatchesByFeatureTracker::process_all_files(int np){
  for (auto i=0; i<matches_filenames.size(); i++)
    results[i]=process_file(i);

}
void MatchesByFeatureTracker::print_results(){
  std::cout<<"-- TOTAL BP --"<<std::endl;
  std::cout<<"order,filename";
  for (auto f: features) std::cout<<","<<f;
  std::cout<<",None"<<std::endl;
  for (auto i=0;i<results.size();i++){
    std::cout<<i<<","<<matches_filenames[i];
    for (auto r: results[i]) std::cout<<","<<r;
    std::cout<<std::endl;
  }
  std::cout<<std::endl<<std::endl<<"-- % BP --"<<std::endl;
  std::cout<<"order,filename";
  for (auto f: features) std::cout<<","<<f;
  std::cout<<",None"<<std::endl;
  for (auto i=0;i<results.size();i++){
    std::cout<<i<<","<<matches_filenames[i];
    for (auto j=0; j< results[i].size();++j) std::cout<<","<<((double) results[i][j])/feature_totals[j];
    std::cout<<std::endl;
  }
}


int main(int argc, char *argv[]){
  char* gff3_filename=argv[1];
  std::vector<std::string> features,filenames;
  for (auto i=2; i<9; i++) {
    if (std::string(argv[i])!="-") features.push_back(std::string(argv[i]));
  }
  for (auto i=9; i<argc; i++)
    filenames.push_back(std::string(argv[i]));

  MatchesByFeatureTracker mbft=MatchesByFeatureTracker(gff3_filename,features,filenames);
  mbft.process_all_files(8);
  mbft.print_results();
}
