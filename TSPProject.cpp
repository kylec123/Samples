#Solving the traveling salesman problem with simulated annealing and genetic algorithm.
#include <string>
using std::string;
#include <cstring>
#include <iostream>
using std::cout;
using std::cin;
using std::scientific;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Globals
int** dist; //pointer to distance matrix
string DISTANCEMATRIXFILE = "distMatrix.txt"; //matrix file name
int numOfCities = 9;

//GA globals
const int populationSize = 5000;
double epsilon = 0.01;
double mutationChance = 0.05;
int pathLen2[populationSize];//for pathLength function, global to be able to access this to calculate average of pathlengths

//SA Globals
double startTemp = 1e50;
double endTemp = 0.001;
double coolingFactor = 0.9999;

void getDistMatrix(string fileName)//reads in matrix, works for arbritrarily sized nxn input matrices. 
{
	ifstream file;
	file.open(fileName.c_str());
	int size=0;
	string line;//to temp store the line read using getline
	getline(file,line);
	char* line1 = new char[line.length()+1];//convert to c string 
	strcpy(line1,line.c_str());
	for(int i=0;i<strlen(line1);i++)
	{
		if(line1[i] == ' ')//each time a blankspace is reached increment the column size of the matrix
			size++;
	}
	file.clear();//clear EOF flag
	file.seekg(0);//return to beginning of file to start reading in the matrix.
	dist = new int*[size];//declare memory to hold matrix
	for(int i=0;i<size;i++)
		dist[i] = new int[size];
	for(int i=0;i<size;i++)//reads in matrix element by element
		for(int j=0;j<size; j++)
		{
			file >> dist[i][j];
		}
	file.close();
}
int pathLength(int** pop, int i)//path length function for GA
{
	int curDist=0;
	for(int k=0;k<numOfCities-1;k++)
	{
		curDist += dist[pop[i][k]][pop[i][k+1]];
	}
	return curDist;
}

double* normFunction(int** pop)//Normalization function for selection
{
	for(int i=0;i<populationSize;i++)
	{
		pathLen2[i] = pathLength(pop, i);
	}
	double a[populationSize];//1/pathLen[] for each member of population
	for(int i=0;i<populationSize;i++)
		a[i] = 1.0/pathLen2[i];
	double b = 0;  //holds the sum of the inverse distances
	for(int i=0;i<populationSize;i++)
		b+= a[i];
	static double normScore[populationSize];
	for(int i=0;i<populationSize;i++)
		normScore[i] = a[i]/b;
	return normScore;
}

void GeneticAlgorithm()
{
	///////GENETIC ALGORITHM//////////
	int** population = new int* [populationSize];//declare memory for population
	for(int i=0;i<populationSize;i++)
		population[i] = new int[numOfCities];
	for(int i=0;i<populationSize;i++)//initialize all population to 1, 2, 3, 4, etc
	{
		for(int j=0;j<numOfCities;j++)
		{
			population[i][j] = j;
		}
	}
	for(int i=0;i<populationSize;i++)//randomize population
		std::random_shuffle(&population[i][1], &population[i][numOfCities]);
	int shortestDistance = 2147483647;//largest int
	for(int i=0; i<populationSize;i++)
	{
		if(pathLength(population, i)<shortestDistance)
		{
			shortestDistance = pathLength(population, i);
		}
	}
	double* normScore = normFunction(population);//initialize pathLen2 by calling normFunction
	double prevAverage = 0;//get initial average
	for(int i=0;i<populationSize;i++)
		prevAverage+=pathLen2[i];
	prevAverage = prevAverage/(double)populationSize;
	double curAverage = 0;
	cout <<"Shortest path distance computed before GA was: " << shortestDistance << "\n";
	ofstream graphData;
	graphData.open("GraphDataGA.csv");
	clock_t start, end;
	start = clock();
	int iters = 0;
	while(fabs(curAverage-prevAverage) > epsilon)//converge when this condition breaks.
	{
		cout << "Convergence: " << fabs(curAverage-prevAverage) << " vs. " << epsilon << "\n";
		//STEP 1: GENERATE FITNESS SCORES
		normScore = normFunction(population);
		prevAverage = curAverage;
		for(int i=0;i<populationSize;i++)
			curAverage+=pathLen2[i];
		curAverage = curAverage/(double)populationSize;
		graphData << iters << ", " << curAverage << "\n";//output average at each iteration to file, for graph generation.
		//END: STEP 1
		//STEP 2: SELECTION
		
		double sum=0;
		double startingPoint[populationSize];
		startingPoint[0] = 0;
		//std::sort(&normScore[0], &normScore[populationSize-1], std::greater<double>());
		for(int i=1;i<populationSize;i++)//fit normScores to a line between 0 and 1.
		{
			sum+=normScore[i-1];
			startingPoint[i] = sum;
		}
		int selectedPopulation[populationSize][numOfCities];
		for(int i=0;i<populationSize;i++)//selection loop
		{
			double selected = rand()/(double)RAND_MAX;//random number between 0-1, to select a random population member
			int j=0;
			while(j<populationSize-1 && startingPoint[j+1] < selected)
			{
				j++;
			}
			//j now holds the index of the selected population member
			for(int k=0;k<numOfCities;k++)
			{
				selectedPopulation[i][k]=population[j][k];
			}
		}
		//Now use selectedPopulation[][] to preform crossover and mutation
		//END: STEP 2
		//STEP 3: CROSSOVER (using Partially Mapped Crossover)(PMX)
		int crossoverPosition = rand() % (numOfCities-2) + 2;//random crossover point from 2-(n-1)
		for(int i=0;i<populationSize;i+=2)//crossover loop
		{
			int tempChrome[numOfCities];//save original chromosome 1 to use on second crossover
			for(int k=0;k<numOfCities;k++)
				tempChrome[k]=selectedPopulation[i][k];
			for(int k=1;k<crossoverPosition;k++)//first offspring changes
			{
				for(int l=k+1;l<numOfCities;l++)
				{
					if(selectedPopulation[i+1][k]==selectedPopulation[i][l])
					{	
						int temp = selectedPopulation[i][k];
						selectedPopulation[i][k] = selectedPopulation[i][l];
						selectedPopulation[i][l] = temp;
					}
				}
			}
			for(int k=1;k<crossoverPosition;k++)//second offspring changes
			{
				for(int l=k+1;l<numOfCities;l++)
				{
					if(tempChrome[k]==selectedPopulation[i+1][l])
					{
						int temp = selectedPopulation[i+1][k];
						selectedPopulation[i+1][k] = selectedPopulation[i+1][l];
						selectedPopulation[i+1][l] = temp;
					}
				}
			}
			crossoverPosition = rand() % (numOfCities-2) + 2;//choose another random crossover point for next two selected parents
		}
		//END: STEP 3
		//STEP 4: MUTATION
		for(int i=0;i<populationSize;i++)
		{
			double mutate = rand()/(double)RAND_MAX;//randomly generated number between 0-1.
			if(mutate<=mutationChance)//mutate if this is true (similar to the fit to line method used for selection, except only one value to fit)
			{
				int mutatedGene = rand() % (numOfCities-2) + 2;//Choose random position (in the range of 2 to numOfCities-1) in path sequence and then flip it with the previous position value
				int temp = selectedPopulation[i][mutatedGene];
				selectedPopulation[i][mutatedGene]=selectedPopulation[i][mutatedGene-1];
				selectedPopulation[i][mutatedGene-1] = temp;
			}
		}
		//END STEP 4
		for(int i=0;i<populationSize;i++)//update population
		{
			for(int j=0;j<numOfCities;j++)
				population[i][j] = selectedPopulation[i][j];
		}
		iters++;
	}
	graphData.close();
	end = clock();
	cout << "Genetic Algorithm for TSP using "<< iters <<" iterations took " << (end-start)/(double)((CLOCKS_PER_SEC))<< " seconds.\n";
	shortestDistance = 2147483647;//largest int
	int k=0;
	for(int i=0; i<populationSize;i++)
	{	
		if(pathLength(population, i)<shortestDistance)
		{
			shortestDistance = pathLength(population, i);
			k=i;
		}
	}
	cout <<"Shortest path distance computed was: " << shortestDistance << "\n";
	cout << "The path is: \n";
	for(int i=0;i<numOfCities;i++)
		cout << population[k][i] << " ";
	cout << "\n";
	for(int i=0;i<populationSize;i++)
		delete [] population[i];
}

int pathLengthV2(int* tour)
{
	int curDist=0;

	for(int k=0;k<numOfCities-1;k++)
	{
		curDist += dist[tour[k]][tour[k+1]];
	}
	return curDist;
}
void SimulatedAnnealing()
{
	//INITIALIZATION SECTION
	int numOfIterations=0;
	int currentSol[numOfCities];//holds the inital solution at first, then holds the current solution during algorithm execution
	for(int i=0;i<numOfCities;i++)
	{
		currentSol[i] = i;
	}
	std::random_shuffle(&currentSol[1], &currentSol[numOfCities]);//Randomly generated initial tour generated
	int curDist = pathLengthV2(currentSol);
	cout << "Shortest path distance computed before SA: " << curDist << "\n";
	int neighborSol[numOfCities];//array to hold neighbor
	double currentTemp = startTemp;//initialize to starting temp
	clock_t start, end;
	start = clock();
	ofstream graphData;
	graphData.open("GraphDataSA.csv");
	while(currentTemp > endTemp)//Annealing loop
	{
		if(numOfIterations%1000 == 0)
			graphData << numOfIterations << ", " << pathLengthV2(currentSol) << "\n";//output pathLength at each 1000th iteration to file, for graph generation.
		int randSwap1 = rand() % (numOfCities-1) + 1;//random value to swap with, to generate a random neighbor.
		int randSwap2 = rand() % (numOfCities-1) + 1;//random value to swap with, to generate a random neighbor.
		while(randSwap1 == randSwap2)//enforce two different random values
		{
			randSwap1 = rand() % (numOfCities-1) + 1;
		}
		for(int i=0;i<numOfCities;i++)//reset neighbor to current
		{
			neighborSol[i] = currentSol[i];
		}
		neighborSol[randSwap1] = currentSol[randSwap2];//swap values in the neighbor  "generate neighbor"
		neighborSol[randSwap2] = currentSol[randSwap1];
		int neighborDist = pathLengthV2(neighborSol);
		int delE = neighborDist - curDist;
		double probability = exp(-delE/currentTemp);
		double randVal = rand()/(double)RAND_MAX;//randomly generated number between 0-1.
		if(neighborDist < curDist || randVal < probability)//if neighbor is better, accept it.  If not better, accept with probability e^(change in dist/tempature)
		{
			currentSol[randSwap1] = neighborSol[randSwap1];
			currentSol[randSwap2] = neighborSol[randSwap2];
			curDist = pathLengthV2(currentSol);
		}
		
		currentTemp *= coolingFactor;//decrease temp
		numOfIterations++;
	}
	end = clock();
	graphData.close();
	cout << "Simulated Annealing completed after " << numOfIterations << " iterations in "<< (end-start)/(double)((CLOCKS_PER_SEC))
	<< " seconds\n";
	cout <<"Shortest path distance computed was: " << pathLengthV2(currentSol) <<"\n";
	cout << "The path is: \n";
	for(int i=0;i<numOfCities;i++)
		cout << currentSol[i] << " ";
	
}

int main()
{
	srand(time(0));//initialize random number generator with unique seed value.
	getDistMatrix(DISTANCEMATRIXFILE);
	cout << "\n\nGenetic Algorithm\n\n";
	GeneticAlgorithm();
	cout << "\n\nSimulated Annealing\n\n";
	SimulatedAnnealing();
	return 0;
}