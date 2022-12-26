#include "DumpReader.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

vector<vector<double>> distanceMatrix(int natoms, vector<vector<double>> atomCoords)
	{
	 std::vector<std::vector<double>> distances(natoms, std::vector<double>(natoms));
	/*#pragma omp parallel num_threads(4)*/
	{
		for (int atoma=0; atoma<natoms; atoma++)
			{
           		for (int atomb=0; atomb<natoms; atomb++)
                		{
                  		distances[atoma][atomb] = std::sqrt(pow(atomCoords[atoma][0]-atomCoords[atomb][0],2) + pow(atomCoords[atoma][0]-atomCoords[atomb][0],2) +pow(atomCoords[atoma][0]-atomCoords[atomb][0],2));
                		}
        		}
	}
	return distances;
  	}


int main(int argc, char *argv[]){
	if (argc ==0){
	cout << "Please give a .dump file as argument!";
	}
	auto start = high_resolution_clock::now();
	Dump NewDump(argv[1]);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	std::cout << "Read a total of "<< NewDump.nsnapshots << " snapshots in " << duration.count() << " seconds." << std::endl;

	start = high_resolution_clock::now();
	for (int i=0; i< 10; i++)
	{
		vector<vector<double>> Distances = distanceMatrix(NewDump.snapshots[0].natoms, NewDump.snapshots[0].atomCoords);
		stop = high_resolution_clock::now();
	}
	auto duration2 = duration_cast<seconds>(stop - start);
	std::cout << "Calculated 10 distance matrices in  " << duration2.count() << " seconds." << std::endl;
/*
	for (int i = 0; i < NewDump.snapshots[0].natoms; i++)
		{
		for (int j = 0; j < NewDump.snapshots[0].natoms; j++)
		{
			cout << "Element at array[" << i << "][" << j << "]: ";
			cout << NewDump.snapshots[0].distances[i][j]<<endl;
		}
	}
*/
}
