//Brian Valenzi, bv457, 4776793, assignment3

#include <iostream>
#include <fstream>
#include <math.h>

const int FILE_NAME_SIZE = 20;
const char fileName[FILE_NAME_SIZE] = "ass3.txt";
const int TOTAL_NUM_VERTICES = 25;
const int ROOT_INDEX = 1;

using namespace std;
struct Coordinates;

//struct for coordinates used in array to store x, y, hueristic values
//Coordinate[0] == a 
//Coordinate[1] == b etc...
struct Coordinates
{
	int x;
	int y;
	int hueristic;
	Coordinates()
	{
		this-> x = 0;
		this-> y = 0;
		this-> hueristic = 0;
	}
	void setCordinates(int x, int y)
	{
		this-> x = x;
		this-> y = y;
	}
	void setHueristic(int hueristic)
	{
		this->hueristic = hueristic;
	}
};
//banshee doesnt assign default 0 values with arr[x]={};
void assignBoolDefaultValue(bool *arr, int size)
{
	for (int i = 0; i < size; ++i)
	{
  		arr[i] = 0;
	}
}
//banshee doesnt assign default 0 values with arr[x]={};
void assignDefaultValue(int *arr, int size)
{
	for (int i = 0; i < size; ++i)
	{
  		arr[i] = 0;
	}
}
//swap int values used in heap
void swap(int &x, int &y)
{
	int z;

	z = x;
	x = y;
	y = z;
}
//insert root node into heap
void insertHeapRoot(int *heap, int index)
{
	heap[ROOT_INDEX] = index;
}
//insert nodes into heap
void insertIntoHeap(int *heap, int index, int &heapEnding, int totalHeapLength)
{
	if(heapEnding < totalHeapLength)
	{
		heap[heapEnding] = index;
		heapEnding++;
	}
}
//heap sift up for either dijkstra or A*
void siftUp(int *heap, int *distance, int position, Coordinates *gCoords)
{
	if(position == 1 || position == 0)
	{
		return;
	}
	else
	{
		int temp = position / 2;
		if(gCoords == NULL)
		{
			if(distance[heap[temp]] < distance[heap[position]])
			{
				return;
			}
			else
			{
				swap(heap[temp], heap[position]);
				siftUp(heap, distance, temp, NULL);
			}
		}
		else
		{
			if((distance[heap[temp]]+gCoords[heap[temp]].hueristic) < (distance[heap[position]]+gCoords[heap[position]].hueristic))
			{
				return;
			}
			else
			{
				swap(heap[temp], heap[position]);
				siftUp(heap, distance, temp, gCoords);
			}
		}
	}
}
//heap sift down for either dijkstra or A*
void siftDown(int *heap, int *distance, int position, int heapEnd, Coordinates *gCoords)
{
	int temp = position * 2;
	if(temp > heapEnd || temp+1 > heapEnd)
	{
		return;
	}
	if(gCoords == NULL)
	{
		if(distance[heap[temp]] > distance[heap[temp+1]])
		{
			temp +=1;
		}
		if(distance[heap[position]] > distance[heap[temp]])
		{
			swap(heap[position], heap[temp]);
			siftDown(heap, distance, temp, heapEnd, NULL);
		}
	}
	else
	{
		if((distance[heap[temp]]+gCoords[heap[temp]].hueristic) >(distance[heap[temp+1]]+gCoords[heap[temp]].hueristic))
		{
			temp +=1;
		}
		if((distance[heap[position]]+gCoords[heap[position]].hueristic) > (distance[heap[temp]]+gCoords[heap[temp]].hueristic))
		{
			swap(heap[position], heap[temp]);
			siftDown(heap, distance, temp, heapEnd, gCoords);
		}
	}
}

//calculate hueristic with euclidean straight line distance
int calculateHueristic(Coordinates gCoords[], int index, int destination)
{
	int xCoords =(gCoords[index].x - gCoords[destination].x)*2;
	int yCoords = (gCoords[index].y - gCoords[destination].y)*2;
	return sqrt(abs(xCoords + yCoords));
}

//calculates and stores all nodes hueristic from current node to end goal node read from text file
void setNodeHueristic(Coordinates graphCoordinates[], int nVertices, int destination)
{
	int h = 0;
	for(int i = 0; i < nVertices; i++)
	{
		h = calculateHueristic(graphCoordinates, i, destination);
		graphCoordinates[i].setHueristic(h);
	}
}

void readFile(int adjMatrix[][TOTAL_NUM_VERTICES], Coordinates graphCoordinates[], int &nVertices, char &source, char &destination)
{
	int nEdges = 0;
	int totalFileReads = 0;

	int x =0;
	int y =0;

	int weight = 0;
	int count = 0;

	fstream fileIn;

	fileIn.open(fileName);
	if(!fileIn)
	{
		cerr << "Unable to open "<<fileName<<endl;
		exit(1);
	}
	fileIn>>nVertices>>nEdges;
	totalFileReads = nVertices + nEdges;

	while(!fileIn.eof())
	{
		if(count < nVertices)
		{
			//store node and its x y coordiantes
			fileIn>> source >> x >> y;
			graphCoordinates[source-97].setCordinates(x, y);
			count++;
		}
		else if(count < totalFileReads)
		{
			//store weights in adjacency matrix
			fileIn>> source >> destination >>weight;
			adjMatrix[source-97][destination-97] = weight;
			count++;
		}
		else
		{
			//store start and goal nodes
			fileIn>> source >> destination;
		}
	}
	fileIn.close();
}

void printPath(int parent[], int source, int destination) 
{ 
	//recursive approach since parent array from dijkstra algorithm store previous node, recursilvely traverse backwards
    if (destination == -1) 
    {
        return; 
    }
  
    printPath(parent, source, parent[destination]); 
    if(destination < 25 && destination > 0)
    {
    	cout<<(char)(parent[destination]+97)<<" ";
    } 
} 
void printStats(int parent[], int source, int destination, int distance, int nodeCount)
{
	cout<<"\tPath: ";
	printPath(parent, source, destination);
	//printPath does not print end vertice must implicitly print it 
    cout<<(char)(destination+97)<<endl;
    cout<<"\tDistance of the path: "<<distance<<endl;
    cout<<"\tTotal vertices visited: "<<nodeCount<<endl;
}

void dijkstra(int adjMatrix[][TOTAL_NUM_VERTICES], int nVertices, int source, int destination, int *parent, int &dist, int &nodesTraversed)
{
	int visited[nVertices];
	//banshee doesnt assign default 0 values with arr[x]={};
	assignDefaultValue(visited, nVertices);

	int nodeCount = nVertices;

	int previousNode = 0;

	int heap[nVertices];
	int heapEnd = ROOT_INDEX;

	int distance[nVertices];
	//if false, havent visited
	bool candidateSet[nVertices];
	//banshee doesnt assign default 0 values with arr[x]={};
	assignBoolDefaultValue(candidateSet, nVertices);

	//int parent[nVertices];
	int currentVertice = 0;

	for(int i = 0; i < nVertices; i++)
	{
		//assigned max int for easier value comparisons
		if(adjMatrix[source][i] == 0)
		{
			distance[i] = INT_MAX;	
		}
		else
		{
			//store weights to other nodes from current index
			distance[i] = adjMatrix[source][i];
			//insert into heap to have smallest weight at root
			insertIntoHeap(heap, i, heapEnd, nVertices);
			siftUp(heap, distance, heapEnd-1, NULL);
		}
		parent[i] = 0;
	}
	visited[0] = source;
	//start node
	nodeCount--;
	while(nodeCount>0)
	{
		//smallest edge weight at heap root
		currentVertice = heap[ROOT_INDEX];
		nodeCount--;
		visited[nVertices-nodeCount] = currentVertice;
		candidateSet[currentVertice] = true;
		//maintain min heap property
		swap(heap[ROOT_INDEX], heap[heapEnd-1]);
		heapEnd--;
		siftDown(heap, distance, ROOT_INDEX, heapEnd-1, NULL);

		if(currentVertice == destination)
		{
			parent[source] = -1;
			//cout<<"Dijkstra shortest path: "<<endl;
			dist = distance[currentVertice];
			nodesTraversed = nVertices - nodeCount;
			//printStats(parent, source, destination, distance, currentVertice, nVertices - nodeCount);
			return;
		}
		for(int i = 0; i < nVertices; i++)
		{
			if(!candidateSet[i])
			{
				if(adjMatrix[currentVertice][i] != 0)
				{
					if(distance[i] > distance[currentVertice] + adjMatrix[currentVertice][i])
					{
						distance[i] = distance[currentVertice] + adjMatrix[currentVertice][i];
						insertIntoHeap(heap, i, heapEnd, nVertices);
						siftUp(heap, distance, heapEnd-1, NULL);
						parent[i] = currentVertice;
						previousNode = currentVertice;
					}
				}
			}
		}
	}
}

void A_Star(int adjMatrix[][TOTAL_NUM_VERTICES], int nVertices, int source, int destination, Coordinates graphCoordinates[], int *parent, int &dist, int &nodesTraversed)
{
	int visited[nVertices];
	assignDefaultValue(visited, nVertices); 
	int nodeCount = nVertices;

	int heap[nVertices];
	int heapEnd = ROOT_INDEX;

	//0 == infinite
	int distance[nVertices];
	assignDefaultValue(distance, nVertices);
	//if false, havent visited
	bool candidateSet[nVertices];
	assignBoolDefaultValue(candidateSet, nVertices); 
	//int parent[nVertices];
	int currentVertice = 0;

	for(int i = 0; i < nVertices; i++)
	{
		//assigned max int for easier value comparisons
		if(adjMatrix[source][i] == 0)
		{
			distance[i] = INT_MAX;	
		}
		else
		{
			//store weights to other nodes from current index
			distance[i] = adjMatrix[source][i];
			//insert into heap to have smallest weight at root
			insertIntoHeap(heap, i, heapEnd, nVertices);
			siftUp(heap, distance, heapEnd-1, graphCoordinates);
		}
		parent[i] = 0;
	}
	visited[0] = source;
	nodeCount--;
	while(nodeCount>0)
	{
		//smallest edge weight at heap root
		currentVertice = heap[ROOT_INDEX];
		nodeCount--;
		visited[nVertices-nodeCount] = currentVertice;
		candidateSet[currentVertice] = true;
		//maintain min heap property
		swap(heap[ROOT_INDEX], heap[heapEnd-1]);
		heapEnd--;
		siftDown(heap, distance, ROOT_INDEX, heapEnd-1, graphCoordinates);

		if(currentVertice == destination)
		{
			parent[source] = -1;
			dist = distance[currentVertice];
			nodesTraversed = nVertices - nodeCount;
			return;
		}
		for(int i = 0; i < nVertices; i++)
		{
			if(!candidateSet[i])
			{
				if(adjMatrix[currentVertice][i] != 0)
				{
					if(distance[i] > distance[currentVertice] + adjMatrix[currentVertice][i])
					{
						distance[i] = distance[currentVertice] + adjMatrix[currentVertice][i];
						insertIntoHeap(heap, i, heapEnd, nVertices);
						siftUp(heap, distance, heapEnd-1, graphCoordinates);
						parent[i] = currentVertice;
					}
				}
			}
		}
	}
}
//calculates second shortest path by removing an edge from the matrix, find shortest path store it then compares it to next iterationm if next iteration is shorter save its values.
//also writes removed edges back into matrix after each path finding iteration
void calcSecondShortestPath(int second_adjMatrix[TOTAL_NUM_VERTICES][TOTAL_NUM_VERTICES], int source, int destination, int nVertices, int shortestPathParent[], Coordinates graphCoordinates[], char flag)
{
	int secondShortestPathParent[nVertices];
	int secondShortestPathDist = 0;

	int count = 0;
	int index = destination;
	int placeHolder = INT_MAX;
	int arrBuffer[nVertices];
	int nodeCountBuffer=0;
	int nodeCount = 0;

	int buffer[nVertices];
	buffer[0] = index;

	//get the path from shortest path nodes traversed and store sequentially in an array to easily iterate through in algorithm below
	for(int i = 1; i < nVertices; i++)
	{
		if(shortestPathParent[index] == -1){count = i; break;}
		buffer[i] = shortestPathParent[index];
		index = shortestPathParent[index];
	}

	int weightBuffer = -1;
	for(int i = 0; i < count-1; i++)
	{
		if(weightBuffer != -1)
		{
			//replace the edge from previous iteration
			second_adjMatrix[buffer[i]][buffer[i+1]] = weightBuffer;
			second_adjMatrix[buffer[i+1]][buffer[i]] = weightBuffer;
		}
		weightBuffer = second_adjMatrix[buffer[i]][buffer[i+1]];
		//remove the edge
		second_adjMatrix[buffer[i]][buffer[i+1]] = 0;
		second_adjMatrix[buffer[i+1]][buffer[i]] = 0;
		if(flag == 'D')
		{
			dijkstra(second_adjMatrix, nVertices, source, destination, secondShortestPathParent, secondShortestPathDist, nodeCount);
		}
		else
		{
			A_Star(second_adjMatrix, nVertices, source, destination, graphCoordinates, secondShortestPathParent, secondShortestPathDist, nodeCount);
		}
		//check if values are current second shortest path if so store if not keep prior values
		if(secondShortestPathDist < placeHolder && secondShortestPathDist != 0)
		{
			placeHolder = secondShortestPathDist;
			memcpy(arrBuffer, secondShortestPathParent, sizeof(secondShortestPathParent));
			nodeCountBuffer = nodeCount;
		}
	}
	if(flag == 'D')
	{
		cout<<"\nDijkstra's second shortest path: "<<endl;
	}
	else
	{
		cout<<"\nA* shortest path: "<<endl;
	}
	printStats(arrBuffer, source, destination, placeHolder, nodeCountBuffer);
}

int main()
{
	int nVertices = 1;
	char source;
	char destination;

	int shortestPathDist = 0;
	int nodeCount = 0;

	int adjMatrix[TOTAL_NUM_VERTICES][TOTAL_NUM_VERTICES]={};
	int second_adjMatrix[TOTAL_NUM_VERTICES][TOTAL_NUM_VERTICES]={};
	Coordinates graphCoordinates[TOTAL_NUM_VERTICES];

	readFile(adjMatrix, graphCoordinates, nVertices, source, destination);
	cout<<endl<<"Start and end vertices: "<<source<< " "<<destination<<endl;
	//copy the matrix for second shortest paths (removing edges and replace after the iteration)
	memcpy(second_adjMatrix, adjMatrix, sizeof(adjMatrix));

	int shortestPathParent[nVertices];

	source-=97;
	destination-=97;
	//set the heuristics for struct Coordinates array containing all nodes x,y,heuristic values
	setNodeHueristic(graphCoordinates, nVertices, destination);

	//shortest path Dijkstra
	dijkstra(adjMatrix, nVertices, source, destination, shortestPathParent, shortestPathDist, nodeCount);
	cout<<"\nDijkstra's shortest path: "<<endl;
	printStats(shortestPathParent, source, destination, shortestPathDist, nodeCount);
	//second shortest path Dijkstra
	calcSecondShortestPath(adjMatrix, source, destination, nVertices, shortestPathParent, graphCoordinates, 'D');

	//A* shortest path
	A_Star(second_adjMatrix, nVertices, source, destination, graphCoordinates, shortestPathParent, shortestPathDist, nodeCount);
	cout<<"\nA* shortest path: "<<endl;
	printStats(shortestPathParent, source, destination, shortestPathDist, nodeCount);
	//second shortest path A*
	calcSecondShortestPath(second_adjMatrix, source, destination, nVertices, shortestPathParent, graphCoordinates, 'A');
}
/*
Adjacency matrix (2D array) was used to store the graph structure. An array of struct Coordinates was used to store the position of each vertices x, y coordinates and the straight line heuristic from that node to the goal node.
Where a = 0, b = 1 etc.. as the array index of the struct array.
A min heap was used to track next node to traverse to keeping the smallest edge weight for dijkstra and smallest edge weight + heuristic for A* at the root node.
To calculate the second shortest path, the shortest path nodes are stored and each loop iteration removes one of the edges from the cloned adjacency matrix and the shortest path is then evaluated with the edge removed then writes the edge back in after
for next edge removal iteration. 
The path and distance is saved if it is the current shortest until all the edges from the original shortest path have all been exhausted, after this linear edge cancelling the second shortest path should have been found. 
This method is used for both Dijkstra and A*.
#QUESTIONS#
It can work for some cases although! issue is if the graph has many paths from source to destination with same total weight

1.No, the second shortest path could still traverse the same amount of accumulated weight just different path in some test cases
2.If the graph contains cycles provided they are not negative cycles it can still find same accumulated weight path
3.even undirected it can still find same accumulated weight as shortest path

a solution could be an added condition checking if the path is greater than first shortest path and the least greatest of all second shortest paths calculated
*/
