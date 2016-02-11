#include "Matrix.h"
#include <vector>
#include <utility>
#include <stdlib.h>
#include <algorithm> 
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cblas.h>
#include <thread>

class Graph{
	public:
		Graph(int N, int M, std::vector<std::pair<int,int>>& edges){
			V_SZ = N;
			E_SZ = M;
			m_adj.construct(V_SZ, V_SZ, false, false);
			for (auto & e : edges){
				int u = e.first; int v = e.second;
				m_adj.set(u, v, true);
				m_adj.set(v, u, true);
			}
		}

		Graph(int N, int M, Matrix<bool>& adj){
			V_SZ = N;
			E_SZ = M;
			m_adj = adj;
		}

		Graph(){
			V_SZ = 0;
			E_SZ = 0;
		}

		Graph(const Graph& g2){
			copy(g2);
		}

		Graph &operator=(const Graph &g2) {
			copy(g2);
			return *this;
		}

		void copy(const Graph &g2){
			V_SZ = g2.V_SZ;
			E_SZ = g2.E_SZ;
			m_adj = g2.getAdj();
		}

		void genGraph(int N, int M){
			std::vector<std::pair<int,int>> all_edges;
			for (int i = 0; i < N; i++){
				for (int j = i + 1; j < N; j++){
					all_edges.push_back({i, j});
				}
			}
			std::random_shuffle(all_edges.begin(), all_edges.end());
			
			std::vector<std::pair<int,int>> edges;
			for (int i = 0; i < M; i++){
				edges.push_back(all_edges[i]);
			}
			copy(Graph(N, M, edges));
		}

		void genRandGraph(int N, int type = 0){
			int total_edges = (N * (N-1))/2;
			int ten_percent = total_edges/10;
			srand (time(NULL));
			int e_rand = rand() % ten_percent;
			switch (type){
				case 1: // Sparse
					break;
				case 2: // Dense
					e_rand += 0.9 * total_edges;
					break;
				default:
					e_rand = rand() % total_edges;
			}
			genGraph(N, e_rand);
		}

/**************** FLOYD WARSHALL *********************/

		// Complexity O(N^3)
		Matrix<int>& FloydWarshall(){
			setup_dist();

			for (int k = 0; k < V_SZ; k++){
				for (int i = 0; i < V_SZ; i++){
					for (int j = 0; j < V_SZ; j++){
						int new_dist = m_dist(i, k) + m_dist(k, j);
						if (m_dist(i, j) > new_dist){
							m_dist(i, j) = new_dist;
						}
					}
				}
			}

			return m_dist;
		}

/*********** FLOYD WARSHALL PARALLEL***************/

		// Complexity O(N^3)
		Matrix<int>& FloydWarshallParallel(){
			setup_dist();

			for (int k = 0; k < V_SZ; k++){
				cilk_for (int i = 0; i < V_SZ; i++){
					for (int j = 0; j < V_SZ; j++){
						int new_dist = m_dist(i, k) + m_dist(k, j);
						if (m_dist(i, j) > new_dist){
							m_dist(i, j) = new_dist;
						}
					}
				}
			}

			return m_dist;
		}

/******************* MIN PLUS **********************/

		// Complexity O(N^3)
		void MatrixMultiplication(Matrix<int>& dest, const Matrix<int>& orig){
			for (int i = 0; i < V_SZ; i++){
				for (int j = 0; j < V_SZ; j++){
					for (int k = 0; k < V_SZ; k++){
						dest(i, j) = std::min(dest(i,j), orig(i, k) + orig(k, j));
					}
				}
			}
		}

		// Complexity O(N^3 * log N)
		Matrix<int>& MinPlus(){
			setup_dist();

			int log_base2 = 0;
			int N = V_SZ;
			while (N >>= 1) log_base2++;

			for (int k = 0; k < log_base2; k++){
				Matrix<int> aux = m_dist;
				MatrixMultiplication(m_dist, aux);
			}

			return m_dist;
		}

/***************** MIN PLUS PARALLEL *****************/


		// Complexity O(N^3)
		void MatrixMultiplicationParallel(Matrix<int>& dest, const Matrix<int>& orig){
			cilk_for (int i = 0; i < V_SZ; i++){
				for (int j = 0; j < V_SZ; j++){
					for (int k = 0; k < V_SZ; k++){
						dest(i, j) = std::min(dest(i,j), orig(i, k) + orig(k, j));
					}
				}
			}
		}

		// Complexity O(N^3 * log N)
		Matrix<int>& MinPlusParallel(){
			setup_dist();

			int log_base2 = 0;
			int N = V_SZ;
			while (N >>= 1) log_base2++;

			for (int k = 0; k < log_base2; k++){
				Matrix<int> aux = m_dist;
				MatrixMultiplicationParallel(m_dist, aux);
			}

			return m_dist;
		}

/***************** MIN PLUS OPTIMIZED *****************/

	int g_cacheBlockSize = 32;
		// Matrix Multiplication with loop unrolling x5 and blocking.
	void MatrixMultiplicationOptimized(Matrix<int>& dest, const Matrix<int>& orig){
		int i, j, k, l;
		int N = dest.dim1();
		int limit0 = N; 			// Index i limit 
		int limit1 = N; 			// Index j limit
		int limit2 = N; 			// Index k limit
		int aux_i, aux_j, aux_k;
		int aux_limit_i; 	 			// Block index limit i
		int aux_limit_j; 	 			// Block index limit j
		int aux_limit_k; 	 			// Block index limit k
		int unroll_factor = 5;
		int unroll_limit; 	 			// Loop unroll index limit

		
		for (i = 0; i < limit0; i += g_cacheBlockSize) {
		    // Blocking index i limit
		    aux_limit_i = std::min((i+g_cacheBlockSize), limit0);
		    
			for (j = 0; j < limit1; j += g_cacheBlockSize) {
			    // Blocking index j limit
			    aux_limit_j = std::min((j+g_cacheBlockSize), limit1);
			    
				for (k = 0; k < limit2; k += g_cacheBlockSize) {
				    // Blocking index k limit
				    aux_limit_k = std::min((k+g_cacheBlockSize), limit2);
				    
		            unroll_limit = aux_limit_k - (unroll_factor-1); // Unrolling by factor of 5
		            
		          	for(aux_i = i; aux_i < aux_limit_i; ++aux_i) {
		            	for(aux_j = j; aux_j < aux_limit_j; ++aux_j) {

		                	int acc0 = INF; int acc1 = INF; int acc2 = INF; int acc3 = INF; int acc4 = INF;
						
											// Unrolling for k loop
		                	for(aux_k = k; aux_k < unroll_limit; aux_k+=unroll_factor) {
		                    	acc0 = std::min(acc0, (orig(aux_i, aux_k) + orig(aux_k, aux_j)));
		                    	acc1 = std::min(acc1, (orig(aux_i, aux_k+1) + orig(aux_k+1, aux_j)));
		                    	acc2 = std::min(acc2, (orig(aux_i, aux_k+2) + orig(aux_k+2, aux_j)));
		                    	acc3 = std::min(acc3, (orig(aux_i, aux_k+3) + orig(aux_k+3, aux_j)));
		                    	acc4 = std::min(acc4, (orig(aux_i, aux_k+4) + orig(aux_k+4, aux_j)));
		                    } 
		                    
		                    // Gather possible uncounted elements
		                    for (; aux_k < aux_limit_k; ++aux_k)
		                    	dest(aux_i, aux_j) = std::min(dest(aux_i, aux_j), (orig(aux_i, aux_k) + orig(aux_k, aux_j)));   
		                   
		                    // Sum up everything
		                    dest(aux_i, aux_j) = std::min(dest(aux_i, aux_j), std::min(acc0, std::min(acc1, std::min(acc2, std::min(acc3, acc4))))); 
						
		            	}  
		        	}   
				}
			}
		}
		return;
	}

			// Complexity O(N^3 * log N)
		Matrix<int>& MinPlusOptimized(){
			setup_dist();

			int log_base2 = 0;
			int N = V_SZ;
			while (N >>= 1) log_base2++;

			for (int k = 0; k < log_base2; k++){
				Matrix<int> aux = m_dist;
				MatrixMultiplicationOptimized(m_dist, aux);
			}

			return m_dist;
		}


		void setup_dist(){
			Matrix<int> new_dist (V_SZ, V_SZ);
			for (int i = 0; i < V_SZ; i++){
				for (int j = 0; j < V_SZ; j++){
					new_dist(i, j) = (m_adj(i, j) == true) ? 1 : INF;
				}
			}
			m_dist = new_dist;
		}

		const Matrix<bool>& getAdj() const{
			return m_adj;
		}

		const Matrix<int>& getDist() const{
			return m_dist;
		}

		const int getVert() const{
			return V_SZ;
		}

	private:
		int V_SZ, E_SZ;
		Matrix<bool> m_adj;
		Matrix<int> m_dist;
};
