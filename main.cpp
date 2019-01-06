
/*
Student Name: Barış Ege Sevgili
Student Number: 2015400084
Compile Status: Compiling
Program Status: Working
Notes: I implemented the project as in the first approach (the bonus part not included)
*/



#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#define NUM_PIXELS 200				//as we are given 200*200 images, the NUM_PIXELS is set to 200


using namespace std;


int main(int args, char* argv[]) {

	int _size;				//#of processes (defined from command line)
	int N;					
	int rows_per_slave;
	int world_rank;
    	
	// Starting up MPI
	MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
	
	N = _size-1; 				//number of slaves 'N'
	rows_per_slave = NUM_PIXELS / N;	//number of rows assigned to each slave
	
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	double beta = stod(argv[3]);		//we take the value from command line
        double pi = stod(argv[4]);		//we take the value from command line
        double delta_E;
        double gamma = 0.5 * log((1 - pi) / pi); //calculations using the formula

	// if current process is the master process
	if(world_rank == 0) {
		
		// X denotes the initial image (noisy), we set the array according to the given input file, by arguments
		int X[NUM_PIXELS][NUM_PIXELS];
			
		//in order to take the input file (noisy image), storing it in X
		ifstream in_file(argv[1], ios::in);
	    	for (int i = 0; i < NUM_PIXELS; i++) {
	    		for (int j = 0; j < NUM_PIXELS; j++) {
	    		        in_file >> X[i][j];
	    		}
	    	}

	    	// Sending the parts of the noisy image to the related slaves
	    	for(int i = 0 ; i < N ; i++) { // n turns 0; i< N
	    		//MPI_Send(address of data start, size of data, MPI_INT, target, tag, MPI_COMM_WORLD)
	    		MPI_Send(X[(i)*rows_per_slave], rows_per_slave * NUM_PIXELS, MPI_INT, i+1, 0, MPI_COMM_WORLD);
	    	}

		// Z denotes the image that we make our changes on
		int Z[NUM_PIXELS][NUM_PIXELS];

		// master process receives the Z_segments that are coming from slaves, after the method is applied
		//collecting these segments in Z array
	    	for(int i = 1 ; i <= N ; i++) {
	    		//MPI_Recv(address of data start, size of data, MPI_INT, source, tag, MPI_COMM_WORLD)
	    		MPI_Recv(Z[rows_per_slave*(i - 1)], rows_per_slave * NUM_PIXELS, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    	}
		
		//printing out the converted text (denoised text(which denotes image)) to the output file
	    	ofstream out_file(argv[2], ios::out);
	    	for(int i = 0; i < NUM_PIXELS; i++){
	        	for(int j = 0; j < NUM_PIXELS; j++){
	            	out_file << Z[i][j] << " ";
	        	}
	        out_file << endl;
	    	
		}
    	}	
    	else{	//current process is a slave
    	
    		int X_segment[rows_per_slave][NUM_PIXELS];		//array to keep the segment of the initial image(text), sent from master process
    		int Z_segment[rows_per_slave][NUM_PIXELS];		//the array which we apply our method on and we make changes on it accordingly
		//creating arrays to handle upper and lower rows of the slaves (neighbour rows) 
	        int neighbour_up[NUM_PIXELS] ;
		int neighbour_down[NUM_PIXELS];

		//slave receives the segment of the X array, which is sent from master process, to each slave specifically
		//keeps that data in the X_segment array
		MPI_Recv(X_segment, rows_per_slave * NUM_PIXELS, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		//initializing Z_segment array, the same as X_segment (kind of copying X to Z)
		//since we need X_segment as it is first initialized later on, we created Z_segment array to be able to make changes on
        	for(int i = 0; i <rows_per_slave; i++){ 
        		for(int k=0; k<NUM_PIXELS; k++){
        			Z_segment[i][k] = X_segment[i][k];        	    
        		}
        	}
        	
        
        	int T = 500000; 	//defines on how many pixels we apply our method
        	int rand_I, rand_J;	//initializing variables for random pixels
        	srand(time(NULL));	//feeding the random number generator

		//iterates T times(picks a random pixel, may flip it or not)

        	for(int i = 0; i < T; i++) {
	
			//arranges the communication between slaves
        		if(world_rank != N && world_rank!= 1){	// if the slave is not an outer slave (neither uppermost nor undermost) it sends its outer
								// rows to neighbour slaves, receives neighbour rows (neighbour_up&neighbour_down) from neighbour slaves
        		        MPI_Recv(neighbour_up, NUM_PIXELS, MPI_INT, world_rank -1 , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	        	MPI_Send(Z_segment[rows_per_slave-1], NUM_PIXELS, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
        	        	MPI_Recv(neighbour_down, NUM_PIXELS, MPI_INT, world_rank + 1 , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	        	MPI_Send(Z_segment[0], NUM_PIXELS, MPI_INT,world_rank - 1, 0, MPI_COMM_WORLD);
               
        		}else if (world_rank == 1){		//if its the uppermost slave, sends info to down, recieves info from down
        	        	MPI_Send(Z_segment[rows_per_slave-1], NUM_PIXELS, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
        		        MPI_Recv(neighbour_down, NUM_PIXELS, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    	}else{					////if its the undermost slave, sends info to up, recieves info from up
				MPI_Recv(neighbour_up, NUM_PIXELS, MPI_INT, world_rank - 1 , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        	    		
				MPI_Send(Z_segment[0],NUM_PIXELS, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD);
        		       
        	        	
        	    	}
        	    	
        	    	
        	    	// picking a random pixel
        		  
			
			rand_I = rand() % rows_per_slave;
			rand_J = rand() % NUM_PIXELS;

			int row_st= rand_I-1;		//the starting point of vertical iterations
			int row_end= rand_I+1;		//the end point of vertical iterations
			int col_st= rand_J-1;		//the starting point of horizontal iterations
			int col_end=rand_J+1;		//the end point of horizontal iterations

			if(rand_I == 0)				//those if statements are used to prevent outofbounds error, at the edges
		   		row_st=0;
		  	if(rand_I == rows_per_slave-1)
		   		row_end=rows_per_slave-1;
			if(rand_J == 0)
		   		col_st=0;
		   	if(rand_I == NUM_PIXELS-1)
		   		col_end=NUM_PIXELS-1;

			int row_st_init= row_st;		//created to initialize the values in the for loops, to the initial values of corresponding variables
			int row_end_init= row_end;
			int col_st_init= col_st;
			int col_end_init=col_end;
			   		
			//sum of all the (-1) and (+1) values around the selected random pixel (neighbour pixels)   	
		   	int sum=0;
			   	
			//does the iterative summation through the Z_segment array's corresponding entries
		   	for(row_st=row_st_init; row_st <= row_end; row_st++){
		   		for(col_st=col_st_init; col_st<=col_end; col_st++){
				    	sum+= Z_segment[row_st][col_st];				
		   		}
		   	}
			 	
			//if selected pixel is in the first row of a slave, neighbour_up row's corresponding pixels are also taken into account, and added to sum value
			//(in the case of the slave is not the first slave (uppermost portion of the image))
		   	if(world_rank!=1 && rand_I ==0){
		   		for (col_st=col_st_init; col_st <= col_end; col_st++) {
					sum += neighbour_up[col_st];
				}
			}
			//if selected pixel is in the last row of a slave, neighbour_down row's corresponding pixels are also taken into account, and added to sum value			
			//(in the case of the slave is not the last slave (undermost portion of the image))		   	
			if(world_rank!=N && rand_I == (rows_per_slave- 1)){
		   		 for (col_st=col_st_init; col_st <= col_end; col_st++) {
				 	sum += neighbour_down[col_st];
				 }
			}

			//calculating the necessary values, using the formula, if the condition inside the if statement is satisfied, the pixel is flipped
			delta_E = (-2.0) * gamma * X_segment[rand_I][rand_J] * Z_segment[rand_I][rand_J] - 2 * beta * Z_segment[rand_I][rand_J] * (sum - Z_segment[rand_I][rand_J]);
	       	    	if (log((double) rand() / RAND_MAX) < delta_E)
	       	        Z_segment[rand_I][rand_J] = -Z_segment[rand_I][rand_J];
			        
		}
		
		//sending Z to the master 
		MPI_Send(Z_segment, rows_per_slave * NUM_PIXELS, MPI_INT, 0, 0, MPI_COMM_WORLD);
    	}	
    	MPI_Finalize();
    	return 0;
}
