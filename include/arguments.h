#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

struct options{

	std::string GENOTYPE_FILE_PATH;
	std::string OUTPUT_PATH;
	int max_iterations ; 
	int num_of_evec ;
	bool getaccuracy ;
	bool debugmode;
	bool var_normalize;
	int accelerated_em;
	int l;
	double convergence_limit;
	bool memory_efficient;
	bool fast_mode;
	bool missing;
	bool text_version;
};

template<typename T, typename U>
struct is_same {
    static const bool value = false; 
};

template<typename T>
struct is_same<T,T> { 
   static const bool value = true; 
};


extern options command_line_opts;


void exitWithError(const std::string &error) {
	std::cout << error;
	std::cin.ignore();
	std::cin.get();

	exit(EXIT_FAILURE);
}

class Convert{
public:
	template <typename T>
	static std::string T_to_string(T const &val){
		std::ostringstream ostr;
		ostr << val;
		return ostr.str();
	}
		
	template <typename T>
	static T string_to_T(std::string const &val){
		std::istringstream istr(val);
		T returnVal;
		if(is_same<T,bool>::value){
			if (!(istr >> std::boolalpha >> returnVal))
				exitWithError("CFG: Not a valid bool received!\n");
			return returnVal;
		}
		else{
			if (!(istr >> returnVal))
				exitWithError("CFG: Not a valid " + (std::string)typeid(T).name() + " received!\n");
			return returnVal;
		}
	}

	static std::string string_to_T(std::string const &val){
		return val;
	}
};


class ConfigFile{
private:
	std::map<std::string, std::string> contents;
	std::string fName;

	void removeComment(std::string &line) const{
		if (line.find('#') != line.npos)
			line.erase(line.find('#'));
	}

	bool onlyWhitespace(const std::string &line) const{
		return (line.find_first_not_of(' ') == line.npos);
	}
	bool validLine(const std::string &line) const{
		std::string temp = line;
		temp.erase(0, temp.find_first_not_of("\t "));
		if (temp[0] == '=')
			return false;

		for (size_t i = temp.find('=') + 1; i < temp.length(); i++)
			if (temp[i] != ' ')
				return true;

		return false;
	}

	void extractKey(std::string &key, size_t const &sepPos, const std::string &line) const{
		key = line.substr(0, sepPos);
		if (key.find('\t') != line.npos || key.find(' ') != line.npos)
			key.erase(key.find_first_of("\t "));
	}
	void extractValue(std::string &value, size_t const &sepPos, const std::string &line) const{
		value = line.substr(sepPos + 1);
		value.erase(0, value.find_first_not_of("\t "));
		value.erase(value.find_last_not_of("\t ") + 1);
	}

	void extractContents(const std::string &line){
		std::string temp = line;
		temp.erase(0, temp.find_first_not_of("\t "));
		size_t sepPos = temp.find('=');

		std::string key, value;
		extractKey(key, sepPos, temp);
		extractValue(value, sepPos, temp);

		if (!keyExists(key))
			contents.insert(std::pair<std::string, std::string>(key, value));
		else
			exitWithError("CFG: Can only have unique key names!\n");
	}

	void parseLine(const std::string &line, size_t const lineNo){
		if (line.find('=') == line.npos)
			exitWithError("CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n");

		if (!validLine(line))
			exitWithError("CFG: Bad format for line: " + Convert::T_to_string(lineNo) + "\n");

		extractContents(line);
	}

	void ExtractKeys(){
		std::ifstream file;
		file.open(fName.c_str());
		if (!file)
			exitWithError("CFG: File " + fName + " couldn't be found!\n");

		std::string line;
		size_t lineNo = 0;
		while (std::getline(file, line)){
			lineNo++;
			std::string temp = line;

			if (temp.empty())
				continue;

			removeComment(temp);
			if (onlyWhitespace(temp))
				continue;

			parseLine(temp, lineNo);
		}

		file.close();
	}
public:
	ConfigFile(const std::string &fName){
		this->fName = fName;
		ExtractKeys();
	}

	bool keyExists(const std::string &key) const{
		return contents.find(key) != contents.end();
	}

	template <typename ValueType>
	ValueType getValueOfKey(const std::string &key, ValueType const &defaultValue = ValueType()) const{
		if (!keyExists(key))
			return defaultValue;

		return Convert::string_to_T<ValueType>(contents.find(key)->second);
	}
};

void parse_args(int argc, char const *argv[]){
	
	// Setting Default Values
	command_line_opts.num_of_evec=2;
	command_line_opts.max_iterations=2*command_line_opts.num_of_evec;
	command_line_opts.getaccuracy=false;
	command_line_opts.debugmode=false;
	command_line_opts.OUTPUT_PATH = "";
	bool got_genotype_file=false;
	command_line_opts.var_normalize=false;
	command_line_opts.l=2;
	command_line_opts.accelerated_em=0;
	command_line_opts.convergence_limit= -1.0;
	command_line_opts.memory_efficient=false;
	command_line_opts.fast_mode=true;
	command_line_opts.missing=false;
	command_line_opts.text_version = false;
	

	if(argc<3){
		cout<<"Correct Usage is "<<argv[0]<<" -p <parameter file>"<<endl;
		exit(-1);
	}

	if(strcmp(argv[1],"-p")==0){

		std::string cfg_filename = std::string(argv[2]);
		ConfigFile cfg(cfg_filename);
		got_genotype_file=cfg.keyExists("genotype");
		command_line_opts.max_iterations = cfg.getValueOfKey<int>("max_iterations",5);
		command_line_opts.num_of_evec=cfg.getValueOfKey<int>("num_evec",2);
		command_line_opts.getaccuracy=cfg.getValueOfKey<bool>("accuracy",false);
		command_line_opts.debugmode=cfg.getValueOfKey<bool>("debug",false);
		command_line_opts.l=cfg.getValueOfKey<int>("l",0);
		command_line_opts.OUTPUT_PATH = cfg.getValueOfKey<string>("output_path",string(""));
		command_line_opts.GENOTYPE_FILE_PATH = cfg.getValueOfKey<string>("genotype",string(""));
		command_line_opts.convergence_limit =cfg.getValueOfKey<double>("convergence_limit",-1.0);
		command_line_opts.var_normalize = cfg.getValueOfKey<bool>("var_normalize",false);
		command_line_opts.accelerated_em = cfg.getValueOfKey<int>("accelerated_em",0);
		command_line_opts.memory_efficient = cfg.getValueOfKey<bool>("memory_efficient",false);	
		command_line_opts.fast_mode = cfg.getValueOfKey<bool>("fast_mode",true);
		command_line_opts.missing = cfg.getValueOfKey<bool>("missing",false);	
		command_line_opts.text_version = cfg.getValueOfKey<bool>("text_version",false);							
	}
	else{
		bool got_max_iter = false;
		for (int i = 1; i < argc; i++) { 
			if (i + 1 != argc){
				if(strcmp(argv[i],"-g")==0){
					command_line_opts.GENOTYPE_FILE_PATH = string(argv[i+1]);
					got_genotype_file=true;
					i++;
				}
				else if(strcmp(argv[i],"-o")==0){
					command_line_opts.OUTPUT_PATH = string(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-k")==0){
					command_line_opts.num_of_evec = atoi(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-m")==0){
					command_line_opts.max_iterations = atoi(argv[i+1]);
					got_max_iter = true;
					i++;
				}
				else if(strcmp(argv[i],"-l")==0){
					command_line_opts.l = atoi(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-cl")==0){
					command_line_opts.convergence_limit = atof(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-aem")==0){
					command_line_opts.accelerated_em = atof(argv[i+1]);
					i++;
				}
				else if(strcmp(argv[i],"-v")==0)
					command_line_opts.debugmode=true;
				else if(strcmp(argv[i],"-vn")==0)
					command_line_opts.var_normalize=true;
				else if(strcmp(argv[i],"-a")==0)
					command_line_opts.getaccuracy=true;
				else if(strcmp(argv[i],"-mem")==0)
					command_line_opts.memory_efficient=true;
				else if(strcmp(argv[i],"-miss")==0)
					command_line_opts.missing=true;
				else if(strcmp(argv[i],"-nfm")==0)
					command_line_opts.fast_mode=false;
				else if(strcmp(argv[i],"-txt")==0)
					command_line_opts.text_version=true;
				
				else{
					cout<<"Not Enough or Invalid arguments"<<endl;
					cout<<"Correct Usage is "<<argv[0]<<" -g <genotype file> -k <num_of_evec> -m <max_iterations> -v (for debugmode) -a (for getting accuracy)"<<endl;
					exit(-1);
				}
			}
			else if(strcmp(argv[i],"-v")==0)
			command_line_opts.debugmode=true;
			else if(strcmp(argv[i],"-a")==0)
				command_line_opts.getaccuracy=true;
			else if(strcmp(argv[i],"-vn")==0)
					command_line_opts.var_normalize=true;
			else if(strcmp(argv[i],"-mem")==0)
					command_line_opts.memory_efficient=true;
			else if(strcmp(argv[i],"-nfm")==0)
					command_line_opts.fast_mode=false;
			else if(strcmp(argv[i],"-miss")==0)
					command_line_opts.missing=true;
			else if(strcmp(argv[i],"-txt")==0)
					command_line_opts.text_version=true;
		}
		if(!got_max_iter)
			command_line_opts.max_iterations = 2*command_line_opts.num_of_evec;

	}
	
	if(got_genotype_file==false){
		cout<<"Genotype file missing"<<endl;
		cout<<"Correct Usage is "<<argv[0]<<" -g <genotype file> -k <num_of_evec> -m <max_iterations> -v (for debugmode) -a (for getting accuracy)"<<endl;
		exit(-1);
	}

}

#endif