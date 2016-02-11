#include "GraphIO.h"
#include <iostream>

int main(){	
	GraphIO R;
	R.read_file_adj("test_adj.txt");
	
	Graph G = R.getGraph();
	//std::cout << G.getAdj() << std::endl;
	__cilkrts_end_cilk();  

  __cilkrts_set_param("nworkers", "2");
  __cilkrts_init();

	for (int i = 0; i < 1000; i++){
		G.FloydWarshallParallel();
	}

/*
	Matrix<int> d1, d2, d3, d4, d5, d6, d7;
	G.FloydWarshall();
	d1 = G.getDist();
	G.FloydWarshallParallel();
	d2 = G.getDist();
	G.MinPlus();
	d3 = G.getDist();
	G.MinPlusParallel();
	d4 = G.getDist();
	G.MinPlusOptimized();
	d5 = G.getDist();
	printf("%d\n", (d1 == d2 && d1 == d3 && d1 == d4 && d1 == d5));*/

	//std::cout << d6 << std::endl;

		/*
	G6 = G5;
	G7 = G5;
	G8 = G5;
	std::cout << G5.getAdj() << std::endl;

	std::cout << G5.FloydWarshall() << std::endl;
	std::cout << G6.MinPlus() << std::endl;

	G7.FloydWarshallParallel();
	G8.MinPlusParallel();

	Matrix<int> d1, d2, d3, d4;
	d1 = G5.getDist();
	d2 = G6.getDist();
	d3 = G7.getDist();
	d4 = G8.getDist();



//	std::cout << (d1 == d2 && d1 == d3 && d1 == d4) << std::endl;*/

	return 0;
}
