#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "storage.h"
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/SVD"


class genotype {
	std::vector< std::vector <bool> > msb;
	std::vector< std::vector <bool> > lsb;
	std::vector<int> columnsum;
	std::vector<double> columnmeans;
	
	public:	

		unsigned char mask;
    	int wordsize;
	    unsigned int unitsperword;
    	int unitsize;
	 	int nrow, ncol;
		unsigned char *gtype;

		int Nsnp, Nindv, Nsegments_hori, segment_size_hori, segment_size_ver, Nsegments_ver;
		int Nbits_hori, Nbits_ver;
		int Nelements_hori, Nelements_ver;
		std::vector< std::vector<int> > p;

		// std::vector< std::vector<unsigned> > p_eff;
		// std::vector< std::vector<unsigned> > q_eff;

		std::vector< std::vector<int> > not_O_j;
		std::vector< std::vector<int> > not_O_i;
		
		void init_means(bool is_missing);

		float get_observed_pj(const std::string &line);
		float get_observed_pj(const unsigned char* line);

		void read_txt_naive(std::string filename,bool allow_missing);
		
		void read_txt_mailman (std::string filename,bool allow_missing);

		void read_plink (std::string filenameprefix, bool allow_missing,bool mailman_mode);
		void read_bed (std::string filename, bool allow_missing, bool mailman_mode )  ;
		void read_bed_mailman (std::string filename,bool allow_missing )  ;
		void read_bed_naive (std::string filename ,bool allow_missing)  ;

		int countlines(std::string filename);
		void set_metadata () ;


		double get_geno(int snpindex,int indvindex,bool var_normalize);
		double get_col_mean(int snpindex);
		double get_col_sum(int snpindex);		
		double get_col_std(int snpindex);		
		void update_col_mean(int snpindex,double value);

		void generate_eigen_geno(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &geno_matrix,bool var_normalize);


};

#endif
