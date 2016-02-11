#include "Graph.h"
#include <string>
#include <fstream>

class GraphIO{
	public:
		void read_file(std::string filename){
			std::ifstream input;
			input.open(filename, std::ifstream::in);

			int N = 0, M = 0;
			std::vector<std::pair<int,int>> edges;
			if (input.is_open()){
				input >> N >> M;
				for (int i = 0; i < M; i++){
					int u, v;
					input >> u >> v;
					edges.push_back({u, v});
				}
			}
			input.close();
			G = Graph(N, M, edges);
		}

		void read_file_adj(std::string filename){
			std::ifstream input;
			input.open(filename, std::ifstream::in);

			int N = 0, M = 0;
			Matrix<bool> adj;
			if (input.is_open()){
				input >> N;
				adj.construct(N, N);
				for (int i = 0; i < N; i++){
					for (int j = 0; j < N; j++){
						bool v; 
						input >> v;
						adj(i, j) = v;
						if (v) M++;
					}
				}
			}
			input.close();
			G = Graph(N, M, adj);
		}

		void write_file_adj(std::string filename){
			std::ofstream output;
			output.open(filename);
			if (output.is_open()) {
				int N = G.getVert();
				Matrix<bool> adj = G.getAdj();
				output << N << '\n';
				for (int i = 0; i < N; i++){
					for (int j = 0; j < N; j++){
						output << adj(i, j) << " ";
					}
				}
			}
			output.close();
		}

		void setGraph(const Graph& g2){
			G = g2;
		}

		Graph& getGraph(){
			return G;
		}

	private:
		Graph G;

};
