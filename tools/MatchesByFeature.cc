#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iomanip>

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
};

int main(int argc, char *argv[]){
  std::ifstream infile;
  std::string line,sline[5];
  std::map<std::string,std::shared_ptr<Chromosome>> chromosomes;
  std::vector<std::string> features;
  //TODO: create a list of the features
  for (auto i=2; i<argc-1; i++)
    features.push_back(std::string(argv[i]));
  //takes a gff file, loads a uint8 with bit-field values for presence for up to 8 features being present or not per position.
  //for each chromosome line
  std::cout<<"Creating chromosome structure... "<<std::flush;
  infile.open(argv[1]);
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
      std::shared_ptr<Chromosome> pc(new Chromosome(sline[0],atol(sline[3].c_str()),atol(sline[4].c_str())));
      chromosomes.insert(std::make_pair(sline[0], pc));
    }
  }
  
  infile.close();
  std::cout<<"DONE!"<<std::endl;
  
  infile.open(argv[1]);

  std::cout<<"Loading features... "<<std::flush;
  while (getline (infile,line) ) {
    if ('#'==line[0]) continue;
    std::stringstream ss(line);
    for (auto i=0; i<5; i++){
      ss>>sline[i];
    }
    std::map<std::string, std::shared_ptr<Chromosome>>::iterator c=chromosomes.find(sline[0]);
    if ( chromosomes.end() != c ){//XXX: this "eliminates" hanging-in-the-air features
      for (auto i=0;i<features.size();i++) {
        if (features[i]==sline[2]) {
          for (auto j=atol(sline[3].c_str());j <= atol(sline[4].c_str()); j++)
          c->second->featuresv[j-c->second->start]|=(1<<i);
          break;
        }
      }
    }
  }
  std::cout<<"DONE!"<<std::endl;
  infile.close();
  std::cout<<"Loading alignement coords... "<<std::flush;
  infile.open(argv[argc-1]);
  while (getline (infile,line) ) {
    std::stringstream ss(line);
    for (auto i=0; i<3; i++) ss>>sline[i];
    std::map<std::string, std::shared_ptr<Chromosome>>::iterator c=chromosomes.find(sline[0]);
    if ( chromosomes.end() != c ){//XXX: this "eliminates" hanging-in-the-air features
          for (auto j=atol(sline[1].c_str());j < atol(sline[2].c_str()); j++)///TODO: check it's 0-based in-between
          c->second->featuresv[j]|=128;
    }
    
  }
  std::cout<<"DONE!"<<std::endl;
  //prints a stat for covered% for each of the features across the genome
  uint64_t h[256];
  std::cout<<"Counting bp coverage... "<<std::flush;
  for (auto i=0;i<256;i++) h[i]=0;//better safe than sorry
  for (auto c=chromosomes.begin();c!=chromosomes.end();++c){
    /*for (auto basecov=c->second->featuresv.begin();basecov!=c->second->featuresv.end();++basecov){
      h[*basecov]++;
    }*/
    uint8_t * bc_data=c->second->featuresv.data();
    for (unsigned long i=0; i<c->second->featuresv.size(); i++){
      h[bc_data[i]]++;
    }
    std::cout<<"+"<<std::flush;
  }
  std::cout<<"DONE!"<<std::endl;
  for (auto i=0;i<256;i++) std::cout <<i<<": "<<h[i]<<std::endl;
  std::cout<<"Breakdown by feature:"<<std::endl;
  for (auto i=0;i<features.size();i++) {
    unsigned long int covered=0,total=0;
    for (auto j=0;j<256;j++){
      if ((1<<i)&j) {
        if (128&j) covered+=h[j];
        total+=h[j];
      }
    }
    std::cout<<features[i]<<": "<<std::fixed<<(100.0*(double)covered/(double)total)<<"% ("<<covered<<"/"<<total<<")"<<std::endl;
    
  }
  std::cout<<"None: "<<std::fixed<<(100.0*(double)h[128]/(double)(h[128]+h[0]))<<"% ("<<h[128]<<"/"<<(h[128]+h[0])<<")"<<std::endl;
  
}
