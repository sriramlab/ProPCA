#include <bits/stdc++.h>
#include "genotype.h"
#include "storage.h"

using namespace std;

void genotype::init_means(bool is_missing){

	columnmeans.resize(Nsnp);
	for(int i=0;i<Nsnp;i++){
		double sum = columnsum[i];
		if(is_missing)
			columnmeans[i] = sum*1.0/(Nindv-not_O_i[i].size());
		else
			columnmeans[i] = sum*1.0/Nindv;
	}
}

void genotype::read_txt_naive (std::string filename,bool allow_missing){
	
	ifstream ifs (filename.c_str(), ios::in);                                       
	
	std::string line;
	std::getline(ifs, line);
    std::istringstream iss(line);
    if (!(iss >> Nsnp >> Nindv)) { 
		cout<<"ERROR: Header with number of SNPs and individuals not present"<<endl; 
		exit(-1);
	}
	
	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}

	int i=0;

	vector <bool> m;
	vector <bool> l;
	while(std::getline(ifs,line)){
		int sum=0;
		for(int j=0;j<line.size();j++){
			int val = int(line[j]-'0');	
			if(val==0){
				l.push_back(false);
				m.push_back(false);
			}
			else if(val==1){
				sum+=1;
				l.push_back(true);
				m.push_back(false);
			}
			else if(val==2){
				sum+=2;
				l.push_back(false);
				m.push_back(true);
			}
			else if(val==9 && allow_missing){
				not_O_i[i].push_back(j);
				not_O_j[j].push_back(i);
				l.push_back(false);
				m.push_back(false);
			}
			else{
				cout<<"Invalid entry in Genotype Matrix"<<endl;
				cout<<"If there is Missing data, run with -miss flag"<<endl;
				exit(-1);				
			}
		}
		i++;
		columnsum.push_back(sum);
		msb.push_back(m);
		lsb.push_back(l);
		m.clear();
		l.clear();
	}
	init_means(allow_missing);
}


void genotype::read_txt_mailman (std::string filename,bool allow_missing){
   	ifstream ifs (filename.c_str(), ios::in);                                       
	
	// Calculating the sizes and other stuff for genotype matrix
	std::string line;
	std::getline(ifs, line);
    std::istringstream iss(line);
    if (!(iss >> Nsnp >> Nindv)) { 
		cout<<"ERROR: Header with number of SNPs and individuals not present"<<endl; 
		exit(-1);
	}
	segment_size_hori = ceil(log(Nindv)/log(3));
	segment_size_ver = ceil(log(Nsnp)/log(3));
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	Nsegments_ver = ceil(Nindv*1.0/(segment_size_ver*1.0));
	p.resize(Nsegments_hori,std::vector<int>(Nindv));

	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}

	int i=0;

	while(std::getline(ifs,line)){
		int horiz_seg_no = i/segment_size_hori ;
		int sum=0;
		for(int j=0;j<line.size();j++){
			int val = int(line[j]-'0');
			if(val==0 || val==1 || val==2){
				sum+=val;
				p[horiz_seg_no][j] = (3 * p[horiz_seg_no][j]) + val ;
			}	
			else if(val==9 && allow_missing){
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j] ;
				not_O_i[i].push_back(j);
				not_O_j[j].push_back(i);
			}
			else{
				cout<<"ERROR: Invalid character in genotype file"<<endl;
				cout<<"If there is Missing data, run with -miss flag"<<endl;				
				exit(-1);
			}			
		}
		i++;
		columnsum.push_back(sum);
	}
	init_means(allow_missing);	
}


template<typename T>
static std::istream & binary_read(std::istream& stream, T& value){
	return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template <class T>
inline void printvector(vector<T> &t, string delim = " ", bool newline = false){
		for (int i = 0; i < t.size(); i++)
				cout << t[i] << delim;
        if (newline)
            cout << endl;
}

template <class T>
inline void printvectornl(vector<T> &t, string delim = " "){
    printvector (t, delim, true);
}


void genotype::read_bim (string filename){
	ifstream inp(filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int j = 0 ;
	int linenum = 0 ;
	while(std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		if (line.empty())
			continue;
		j++;
	}
	Nsnp = j;
	inp.close();
}

void genotype::read_fam (string filename){
	ifstream inp(filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int j = 0 ;
	int linenum = 0 ;
	while(std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		if (line.empty())
			continue;
		j++;
	}
	Nindv = j;
	inp.close();
}

void genotype::set_metadata() { 
	wordsize = sizeof(char) * 8;
	unitsize = 2;
	unitsperword = wordsize/unitsize;
	mask = 0;
	for (int i = 0 ; i < unitsize; i++)
		mask = mask |(0x1<<i);
    nrow = Nsnp;
    ncol = ceil(1.0*Nindv/unitsperword);
}

void genotype::read_bed_mailman_nomissing (string filename )  {
   	ifstream ifs (filename.c_str(), ios::in|ios::binary);                                       
   	char magic[3];
	set_metadata ();

    gtype =  new unsigned char[ncol];

   	binary_read(ifs,magic);

	segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	p.resize(Nsegments_hori,std::vector<int>(Nindv));
	
	int sum=0;
	columnsum.resize (Nsnp);    

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	vector<int> v (Nindv);
	int y[4];
	for (int i = 0 ; i < Nsnp; i++){
		int horiz_seg_no = i/segment_size_hori ;
	   	ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));   
    	for (int k = 0 ;k < ncol ; k++) {
        	unsigned char c = gtype [k];
			// Extract PLINK genotypes
        	y[0] = (c)&mask;
        	y[1] = (c>>2)&mask;
        	y[2] = (c>>4)&mask;
        	y[3] = (c>>6)&mask;
			int j0 = k * unitsperword;
			// Handle number of individuals not being a multiple of 4
			int lmax = 4;
			if (k == ncol - 1)  {
				lmax = Nindv%4;
				lmax = (lmax==0)?4:lmax;
			}	
			// Note  : Plink uses different values for coding genotypes
			// Note  : Does not work for missing values
			// To handle missing data it is recommended to write a separate function. This is easy to do.
			// This will avoid the performance hit of checking for and handling missing values
			for ( int l = 0 ; l < lmax; l++){
				int j = j0 + l ;
				int ver_seg_no = j/segment_size_ver ;
				// Extract  PLINK coded genotype and convert into 0/1/2
				// PLINK coding: 
				// 00->0
				// 01->missing
				// 10->1
				// 11->2
				int val = y[l];
				val-- ; 
				val =  (val < 0 ) ? 0 :val ;
				sum += val;
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j]  + val;
			}

    	}
		columnsum[i] = sum;
		sum = 0 ;
	}
	init_means(false);	

	delete[] gtype;
}

void genotype::read_bed_mailman_missing (string filename )  {

	ifstream ifs (filename.c_str(), ios::in|ios::binary);                                       
   	char magic[3];
	set_metadata ();
    gtype =  new unsigned char[ncol];

   	binary_read(ifs,magic);

	segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	p.resize(Nsegments_hori,std::vector<int>(Nindv));
	not_O_i.resize(Nsnp);
	not_O_j.resize(Nindv);	
	
	int sum=0;
	columnsum.resize (Nsnp);    

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	vector<int> v (Nindv);
	int y[4];
	for (int i = 0 ; i < Nsnp; i++){
		int horiz_seg_no = i/segment_size_hori ;
	   	ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));   
    	for (int k = 0 ;k < ncol ; k++) {
        	unsigned char c = gtype [k];
			// Extract PLINK genotypes
        	y[0] = (c)&mask;
        	y[1] = (c>>2)&mask;
        	y[2] = (c>>4)&mask;
        	y[3] = (c>>6)&mask;
			int j0 = k * unitsperword;
			// Handle number of individuals not being a multiple of 4
			int lmax = 4;
			if (k == ncol - 1)  {
				lmax = Nindv%4;
				lmax = (lmax==0)?4:lmax;
			}	
			// Note  : Plink uses different values for coding genotypes
			// Note  : Does not work for missing values
			// To handle missing data it is recommended to write a separate function. This is easy to do.
			// This will avoid the performance hit of checking for and handling missing values
			for ( int l = 0 ; l < lmax; l++){
				int j = j0 + l ;
				int ver_seg_no = j/segment_size_ver ;
				// Extract  PLINK coded genotype and convert into 0/1/2
				// PLINK coding: 
				// 00->0
				// 01->missing
				// 10->1
				// 11->2
				int val = y[l];
				if(val==1){
					not_O_i[i].push_back(j);
					not_O_j[j].push_back(i);
				}
				val-- ; 
				val =  (val < 0 ) ? 0 :val ;
				sum += val;
				p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j]  + val;
			}
    	}
		columnsum[i] = sum;
		sum = 0 ;
	}
	init_means(true);	

	delete[] gtype;
}

void genotype::read_bed_naive (string filename, bool allow_missing)  {

	ifstream ifs (filename.c_str(), ios::in|ios::binary);                                       
   	char magic[3];
	set_metadata ();
    gtype =  new unsigned char[ncol];

   	binary_read(ifs,magic);

	msb.resize(Nsnp,std::vector<bool>(Nindv));
	lsb.resize(Nsnp,std::vector<bool>(Nindv));

	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}
	
	int sum=0;
	columnsum.resize (Nsnp);    

	// Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
	// allele to code a SNP as 0 or 1.
	// This flipping does not matter for results.
	vector<int> v (Nindv);
	int y[4];
	for (int i = 0 ; i < Nsnp; i++){
	   	ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));   
    	for (int k = 0 ;k < ncol ; k++) {
        	unsigned char c = gtype [k];
			// Extract PLINK genotypes
        	y[0] = (c)&mask;
        	y[1] = (c>>2)&mask;
        	y[2] = (c>>4)&mask;
        	y[3] = (c>>6)&mask;
			int j0 = k * unitsperword;
			// Handle number of individuals not being a multiple of 4
			int lmax = 4;
			if (k == ncol - 1)  {
				lmax = Nindv%4;
				lmax = (lmax==0)?4:lmax;
			}	
			// Note  : Plink uses different values for coding genotypes
			// Note  : Does not work for missing values
			// To handle missing data it is recommended to write a separate function. This is easy to do.
			// This will avoid the performance hit of checking for and handling missing values
			for ( int l = 0 ; l < lmax; l++){
				int j = j0 + l ;
				// Extract  PLINK coded genotype and convert into 0/1/2
				// PLINK coding: 
				// 00->0
				// 01->missing
				// 10->1
				// 11->2
				int val = y[l];
				if(val==1 && allow_missing){
					not_O_i[i].push_back(j);
					not_O_j[j].push_back(i);
				}
				val-- ; 
				val =  (val < 0 ) ? 0 :val ;
				sum += val;
				if(val==0){
					lsb[i][j] = false;
					msb[i][j]= false;
				}
				else if(val==1){
					lsb[i][j]= true;
					msb[i][j]= false;
				}
				else if(val==2){
					lsb[i][j]= false;
					msb[i][j]= true;
				}
				else{
					cout<<"Invalid entry in Genotype Matrix"<<endl;
					cout<<"If there is Missing data, run with -miss flag"<<endl;
					exit(-1);				
				}
			}
    	}
		columnsum[i] = sum;
		sum = 0 ;
	}
	init_means(allow_missing);	

	delete[] gtype;
}


void genotype::read_bed (string filename, bool allow_missing, bool mailman_mode )  {
	if(mailman_mode){
		if (allow_missing)
			read_bed_mailman_missing (filename);
		else
			read_bed_mailman_nomissing (filename);
	}
	else{
		read_bed_naive(filename,allow_missing);
	}
}

void genotype::read_plink(std::string filenameprefix, bool allow_missing,bool mailman_mode)  { 
	
	std::stringstream f1;  
	f1 << filenameprefix << ".bim";
	read_bim (f1.str());
	std::stringstream f2;  
	f2 << filenameprefix << ".fam";
	read_fam (f2.str());
	std::stringstream f3;  
	f3 << filenameprefix << ".bed";
	read_bed (f3.str(), allow_missing,mailman_mode);
}


// Accessor Functions

double genotype::get_geno(int snpindex,int indvindex,bool var_normalize=false){
	double m = msb[snpindex][indvindex];
	double l = lsb[snpindex][indvindex];
	double geno = (m*2.0+l) - get_col_mean(snpindex);
	if(var_normalize)
		return geno/get_col_std(snpindex);
	else
		return geno;
}

double genotype::get_col_mean(int snpindex){
	double temp = columnmeans[snpindex];
	return temp;
}

double genotype::get_col_sum(int snpindex){
	double temp = columnsum[snpindex];
	return temp;
}

double genotype::get_col_std(int snpindex){
	double p_i = get_col_mean(snpindex);
	if(p_i == 0 || p_i == 2)
		return 1.0;
	double temp = sqrt(p_i*(1-(0.5*p_i))) ; 
	return temp;
}

void genotype::generate_eigen_geno(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &geno_matrix,bool var_normalize){
	for(int i=0;i<Nsnp;i++){
		for(int j=0;j<Nindv;j++){
			double m = msb[i][j];
			double l = lsb[i][j];
			double geno = (m*2.0+l) - get_col_mean(i);
			if(var_normalize)
				geno_matrix(i,j) = geno/get_col_std(i);
			else
				geno_matrix(i,j) = geno;
		}
	}
}

// Modifier Functions

void genotype::update_col_mean(int snpindex,double value){
	columnmeans[snpindex] = value;
}



/* Redundant Function
void genotype::read_genotype_eff (std::string filename,bool allow_missing){
	FILE* fp;
	fp= fopen(filename.c_str(),"r");
	int j=0;
	int i=0;
	char ch;
	// Calculating the sizes and other stuff for genotype matrix
	int rd = fscanf(fp,"%d %d\n",&Nsnp,&Nindv);
	segment_size_hori = ceil(log(Nindv)/log(3));
	segment_size_ver = ceil(log(Nsnp)/log(3));
	Nsegments_hori = ceil(Nsnp*1.0/(segment_size_hori*1.0));
	Nsegments_ver = ceil(Nindv*1.0/(segment_size_ver*1.0));
	Nbits_hori = ceil(log2(pow(3,segment_size_hori)));
	Nbits_ver = ceil(log2(pow(3,segment_size_ver)));
	Nelements_hori = floor( (Nindv * Nbits_hori *1.0) / 32) + 1;
	Nelements_ver = floor( (Nsnp * Nbits_ver*1.0) / 32) + 1;
	cout<<Nbits_hori<<"  "<<Nbits_ver<<"  "<<Nelements_hori<<"  "<<Nelements_ver<<endl;
	p_eff.resize(Nsegments_hori,std::vector<unsigned>(Nelements_hori));
	q_eff.resize(Nsegments_ver,std::vector<unsigned>(Nelements_ver));
	int sum=0;
	if(allow_missing){
		not_O_i.resize(Nsnp);
		not_O_j.resize(Nindv);	
	}

    do{
		int rd = fscanf(fp,"%c",&ch);
		if(ch=='\n'){
			i++;
			columnsum.push_back(sum);
			sum=0;
			j=0;
		}
		else{
			int val = int(ch-'0');
			int horiz_seg_no = i/segment_size_hori ;
			int ver_seg_no = j/segment_size_ver ;
			if(val==0){
				int temp = 3* extract_from_arr(j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(temp,j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(3*extract_from_arr(i,Nbits_ver,q_eff[ver_seg_no]),i,Nbits_ver,q_eff[ver_seg_no]);
			}
			else if(val==1){
				sum+=1;
				int temp = 3* extract_from_arr(j,Nbits_hori,p_eff[horiz_seg_no]) + 1;
				add_to_arr(temp,j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(3*extract_from_arr(i,Nbits_ver,q_eff[ver_seg_no]) + 1,i,Nbits_ver,q_eff[ver_seg_no]);
			}
			else if(val==2){
				sum+=2;
				int temp = 3* extract_from_arr(j,Nbits_hori,p_eff[horiz_seg_no]) + 2;
				add_to_arr(temp,j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(3*extract_from_arr(i,Nbits_ver,q_eff[ver_seg_no]) + 2,i,Nbits_ver,q_eff[ver_seg_no]);
			}
			else if(val==9 && allow_missing){
				int temp = 3* extract_from_arr(j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(temp,j,Nbits_hori,p_eff[horiz_seg_no]);
				add_to_arr(3*extract_from_arr(i,Nbits_ver,q_eff[ver_seg_no]),i,Nbits_ver,q_eff[ver_seg_no]);
				not_O_i[i].push_back(j);
				not_O_j[j].push_back(i);				
			}
			else{
				cout<<"Invalid entry in Genotype Matrix"<<endl;
				cout<<"If there is Missing data, run with -miss flag"<<endl;
				exit(-1);
			}
			j++;
		}
	}while(!feof(fp));
	i--;
	init_means(allow_missing);	
}
*/