#ifndef MAILMAN_H
#define MAILMAN_H

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "storage.h"
#include <assert.h>
#include <emmintrin.h>
#include <vector>

namespace mailman {

	/* Compute Y=A*X where A is a m X n matrix represented in mailman form. X is a n X k matrix
 * m : number of rows
 * n : number of columns
 * k : batching factor: how many vectors to operate on 
 * p : matrix A represented in mailman form
 * x : matrix X
 * yint : intermediate computation
 * c : intermediate computation
 * y : result
 	*/
inline 	void fastmultiply_normal(int m, int n , int k, std::vector<int> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y){
		for (int i = 0 ; i < n; i++)  {
			int l = p[i]  ;
			for (int j = 0 ; j < k ; j ++)
				yint[l*k + j] += x(i,j);

		}

		int d = pow(3,m);
		for (int j  = 0 ;  j < m ; j++)  {
			d =d /3;
			for (int l = 0; l < k ; l++)
				c [l] = 0 ; 
			for (int i = 0 ; i < d; i++) { 
				for (int l = 0; l < k ; l++){
					double z1 = yint[l + (i + d)*k];
					double z2 = yint[l + (i +2*d)*k];
					yint[l+(i+d)*k] = 0;
					yint[l+(i+2*d)*k] = 0 ;
					yint[l+i*k] = yint[l+i*k] + z1 + z2;
					c[l] += (z1 + 2*z2);
				}
			}
			for (int l = 0; l < k ; l++)
				y[j][l] = c[l];
		}
		for (int l = 0; l < k ; l++)
			yint[l] = 0;
	}

	/* Compute Y=X*A + Y_0 where A is a m X n matrix represented in mailman form. X is a k X m matrix. Y_0 is a k X n matrix.
 * X is specified as a matrix k X m_0 matrix X_0 (m_0 > m) 
 * and a index: start \in {1..m_0 }
 * so that X = X_0 [start:(start+m-1),]
 *
 * m : number of rows
 * n : number of columns
 * k : batching factor: how many vectors to operate on 
 * start: index into X_0
 * p : matrix A represented in mailman format
 * x : matrix X_0 
 * yint : intermediate computation
 * c : intermediate computation
 * y : result. also contains Y_0 that is updated.
 	*/
inline	void fastmultiply_pre_normal(int m, int n , int k, int start, std::vector<int> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y){
		int size1 = pow(3.,m);
		memset (yint, 0, size1* sizeof(double));

		int prefix = 1 ;
		for (int i  = m - 1 ; i >= 0 ; i--) { 
			int i1 = start + i;
			for (int j = 0 ; j < prefix; j++) {
				int offset0 = j*k;
				int offset1 = (prefix + j )*k ;
				int offset2 = (2 * prefix + j )*k ;
				for (int l = 0 ; l < k ; l++){ 
					yint[offset1  + l] = yint[offset0 + l] + x(i1,l);
					yint[offset2 + l] = yint[offset0 + l] + 2 * x(i1 ,l);
				}
			}
			prefix *= 3;
		}

		for (int i = 0 ; i < n; i++){
			for (int l = 0 ; l < k  ; l++) {
				y[i][l] += yint[l + p[i]*k];
				// yint[l+ p[i]*k] = 0 ;
			}
		}
	}


	/* Compute Y=A*X where A is a m X n matrix represented in mailman form. X is a n X k matrix
 * m : number of rows
 * n : number of columns
 * k : batching factor: how many vectors to operate on 
 * p : matrix A represented in mailman form
 * x : matrix X
 * yint : intermediate computation
 * c : intermediate computation
 * y : result
 	*/
	#ifdef __SSE__
		// k must be a multiple of 10
	void fastmultiply_sse (int m, int n , int k, std::vector<int> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y){
		__m128d x0, x2, x4, x6, x8;
		__m128d y0, y2, y4, y6, y8;
		__m128d z0, z2, z4, z6, z8;

		assert( k%10 == 0 && "k should be a multiple of 10");

		int blocksize = 10;

		for (int i = 0 ; i < n; i++)  {
			int l = p[i]  * k ;

			for (int j = 0 ; j < k ; j+=blocksize){
				int j1 = j + l ;
				x0 = _mm_loadu_pd (&x(i,j));
				x2 = _mm_loadu_pd (&x(i,j + 2));
				x4 = _mm_loadu_pd (&x(i,j + 4));
				x6 = _mm_loadu_pd (&x(i,j + 6));
				x8 = _mm_loadu_pd (&x(i,j + 8));

				y0 = _mm_loadu_pd  (yint+j1);
				y2 = _mm_loadu_pd  (yint+j1+2);
				y4 = _mm_loadu_pd  (yint+j1+4);
				y6 = _mm_loadu_pd  (yint+j1+6);
				y8 = _mm_loadu_pd  (yint+j1+8);

				y0 = _mm_add_pd ( x0, y0);
				y2 = _mm_add_pd ( x2, y2);
				y4 = _mm_add_pd ( x4, y4);
				y6 = _mm_add_pd ( x6, y6);
				y8 = _mm_add_pd ( x8, y8);

				_mm_storeu_pd (&yint[j1], y0);
				_mm_storeu_pd (&yint[j1+2], y2);
				_mm_storeu_pd (&yint[j1+4], y4);
				_mm_storeu_pd (&yint[j1+6], y6);
				_mm_storeu_pd (&yint[j1+8], y8);
			}
		}

		int d = pow(3,m);
		for (int j  = 0 ;  j < m ; j++)  {
			d =d /3;
			for (int l = 0; l < k ; l++)
				c [l] = 0 ; 

			for (int i = 0 ; i < d; i++) { 

				int o1 = i*k; int o2 = (i+d)*k; int o3 = (i+2*d)*k;


				for (int t = 0 ; t < k; t+= blocksize) { 

					int p1 = o1 + t; 
					int p2 = o2 + t;
					int p3 = o3 + t;

					y0 = _mm_load_pd  (&yint[p1]);
					y2 = _mm_load_pd  (&yint[2+p1]);
					y4 = _mm_load_pd  (&yint[4+p1]);
					y6 = _mm_load_pd  (&yint[6+p1]);
					y8 = _mm_load_pd  (&yint[8+p1]);

					x0 = _mm_load_pd  (&yint[p2]);
					x2 = _mm_load_pd  (&yint[2+p2]);
					x4 = _mm_load_pd  (&yint[4+p2]);
					x6 = _mm_load_pd  (&yint[6+p2]);
					x8 = _mm_load_pd  (&yint[8+p2]);

					z0 = _mm_load_pd  (&yint[p3]);
					z2 = _mm_load_pd  (&yint[2+p3]);
					z4 = _mm_load_pd  (&yint[4+p3]);
					z6 = _mm_load_pd  (&yint[6+p3]);
					z8 = _mm_load_pd  (&yint[8+p3]);

					y0 = _mm_add_pd ( y0, x0);
					y2 = _mm_add_pd ( y2, x2);
					y4 = _mm_add_pd ( y4, x4 );
					y6 = _mm_add_pd ( y6, x6 );
					y8 = _mm_add_pd ( y8, x8 );
					y0 = _mm_add_pd ( y0, z0);
					y2 = _mm_add_pd ( y2, z2);
					y4 = _mm_add_pd ( y4, z4 );
					y6 = _mm_add_pd ( y6, z6 );
					y8 = _mm_add_pd ( y8, z8 );

					z0 = _mm_add_pd ( z0, z0);
					z2 = _mm_add_pd ( z2, z2);
					z4 = _mm_add_pd ( z4, z4 );
					z6 = _mm_add_pd ( z6, z6 );
					z8 = _mm_add_pd ( z8, z8 );
					z0 = _mm_add_pd ( z0, x0);
					z2 = _mm_add_pd ( z2, x2);
					z4 = _mm_add_pd ( z4, x4 );
					z6 = _mm_add_pd ( z6, x6 );
					z8 = _mm_add_pd ( z8, x8 );

					x0 = _mm_load_pd  (&c[t+0]);
					x2 = _mm_load_pd  (&c[t+2]);
					x4 = _mm_load_pd  (&c[t+4]);
					x6 = _mm_load_pd  (&c[t+6]);
					x8 = _mm_load_pd  (&c[t+8]);
					z0 = _mm_add_pd ( z0, x0);
					z2 = _mm_add_pd ( z2, x2);
					z4 = _mm_add_pd ( z4, x4 );
					z6 = _mm_add_pd ( z6, x6 );
					z8 = _mm_add_pd ( z8, x8 );

					_mm_store_pd (&yint[p1+0], y0);
					_mm_store_pd (&yint[p1+2], y2);
					_mm_store_pd (&yint[p1+4], y4);
					_mm_store_pd (&yint[p1+6], y6);
					_mm_store_pd (&yint[p1+8], y8);

					_mm_store_pd (&c[t+0], z0);
					_mm_store_pd (&c[t+2], z2);
					_mm_store_pd (&c[t+4], z4);
					_mm_store_pd (&c[t+6], z6);
					_mm_store_pd (&c[t+8], z8);
				
				}
				for (int l = 0; l < k ; l++){
                	yint[l+(i+d)*k] = 0;
					yint[l+(i+2*d)*k] = 0 ;
                }

			}
			for (int l = 0; l < k ; l++)
				y[j][l] = c[l];
		}
		for (int l = 0; l < k ; l++)
			yint[l] = 0;
	}

	


	void fastmultiply_pre_sse (int m, int n , int k, int start, std::vector<int> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y){
		
//		std::cout << "(" << start << "," << start + m << ")" << std::endl;

		int size1 = pow(3.,m);
		memset (yint, 0, size1* sizeof(double));
		__m128d x0, x2, x4, x6, x8;
		__m128d y0, y2, y4, y6, y8;
		__m128d z0, z2, z4, z6, z8;

		assert( k%10 == 0 && "k should be a multiple of 10");
		
		int blocksize = 10;

		int prefix = 1 ;
		for (int i  = m - 1 ; i >= 0 ; i--) { 
			int i1 = start + i;
			for (int j = 0 ; j < prefix; j++) {
				int offset0 = j*k;
				int offset1 = (prefix + j )*k ;
				int offset2 = (2 * prefix + j )*k ;

				for (int l = 0 ; l < k ; l+= blocksize){
					x0 = _mm_loadu_pd (&x(i1,l));
					x2 = _mm_loadu_pd (&x(i1,l + 2));
					x4 = _mm_loadu_pd (&x(i1,l + 4));
					x6 = _mm_loadu_pd (&x(i1,l + 6));
					x8 = _mm_loadu_pd (&x(i1,l + 8));

					y0 = _mm_loadu_pd  (yint+offset0+l);
					y2 = _mm_loadu_pd  (yint+offset0+l+2);
					y4 = _mm_loadu_pd  (yint+offset0+l+4);
					y6 = _mm_loadu_pd  (yint+offset0+l+6);
					y8 = _mm_loadu_pd  (yint+offset0+l+8);
	
 
					y0 = _mm_add_pd ( x0, y0);
					y2 = _mm_add_pd ( x2, y2);
					y4 = _mm_add_pd ( x4, y4);
					y6 = _mm_add_pd ( x6, y6);
					y8 = _mm_add_pd ( x8, y8);

					_mm_storeu_pd (&yint[offset1+l], y0);
					_mm_storeu_pd (&yint[offset1+l+2], y2);
					_mm_storeu_pd (&yint[offset1+l+4], y4);
					_mm_storeu_pd (&yint[offset1+l+6], y6);
					_mm_storeu_pd (&yint[offset1+l+8], y8);


					y0 = _mm_add_pd ( x0, y0);
					y2 = _mm_add_pd ( x2, y2);
					y4 = _mm_add_pd ( x4, y4);
					y6 = _mm_add_pd ( x6, y6);
					y8 = _mm_add_pd ( x8, y8);

					_mm_storeu_pd (&yint[offset2+l], y0);
					_mm_storeu_pd (&yint[offset2+l+2], y2);
					_mm_storeu_pd (&yint[offset2+l+4], y4);
					_mm_storeu_pd (&yint[offset2+l+6], y6);
					_mm_storeu_pd (&yint[offset2+l+8], y8);


				}
			}
			prefix *= 3;
		}


		for (int i = 0 ; i < n; i++){
			int offset  = p[i]*k;
			for (int l = 0 ; l < k  ; l+=blocksize) {
				y0 = _mm_load_pd  (&yint[l+offset]);
				y2 = _mm_load_pd  (&yint[2+l+offset]);
				y4 = _mm_load_pd  (&yint[4+l+offset]);
				y6 = _mm_load_pd  (&yint[6+l+offset]);
				y8 = _mm_load_pd  (&yint[8+l+offset]);

				x0 = _mm_load_pd  (&y[i][l]);
				x2 = _mm_load_pd  (&y[i][2+l]);
				x4 = _mm_load_pd  (&y[i][4+l]);
				x6 = _mm_load_pd  (&y[i][6+l]);
				x8 = _mm_load_pd  (&y[i][8+l]);

				x0 = _mm_add_pd ( y0, x0);
				x2 = _mm_add_pd ( y2, x2);
				x4 = _mm_add_pd ( y4, x4 );
				x6 = _mm_add_pd ( y6, x6 );
				x8 = _mm_add_pd ( y8, x8 );

				_mm_store_pd (&y[i][l], x0);
				_mm_store_pd (&y[i][l+2], x2);
				_mm_store_pd (&y[i][l+4], x4);
				_mm_store_pd (&y[i][l+6], x6);
				_mm_store_pd (&y[i][l+8], x8);
				// y[i][l] += yint[l + p[i]*k];
			}
		}
	}

	#endif

	/* Redundant Function 
	void fastmultiply_memory_eff(int m, int n , int k, std::vector<unsigned> &p, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &x, double *yint, double *c, double **y,int Nbits){
		
		for (int i = 0 ; i < n; i++)  {
			int l = extract_from_arr(i,Nbits,p);
			for (int j = 0 ; j < k ; j ++)
				yint[l*k + j] += x(i,j);

		}

		int d = pow(3,m);
		for (int j  = 0 ;  j < m ; j++)  {
			d =d /3;
			for (int l = 0; l < k ; l++)
				c [l] = 0 ; 
			for (int i = 0 ; i < d; i++) { 
				for (int l = 0; l < k ; l++){
					double z1 = yint[l + (i + d)*k];
					double z2 = yint[l + (i +2*d)*k];
					yint[l+(i+d)*k] = 0;
					yint[l+(i+2*d)*k] = 0 ;
					yint[l+i*k] = yint[l+i*k] + z1 + z2;
					c[l] += (z1 + 2*z2);
				}
			}
			for (int l = 0; l < k ; l++)
				y[j][l] = c[l];
		}
		for (int l = 0; l < k ; l++)
			yint[l] = 0;
	}
	*/

}


#endif
