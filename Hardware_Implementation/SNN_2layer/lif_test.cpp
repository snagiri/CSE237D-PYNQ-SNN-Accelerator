#include <bits/stdc++.h>
using namespace std;
#include<cstring>
#include "lif.hpp"

//#define T 100


float Prest_i=0;
float Pmin_i=-1;
float D_i=0.5; //decay


int tim[101];//
int t_back=-20;
int t_fore=20;
int epoch=1;
int num_of_images=6;
float w_max=0.5;
float w_min=-0.5;

axis_in *train;
axis_in *spike_count;


//void LIF_Network(float train[784][101], int spike_count[8]);

int main()
{
	for(int i=1;i<101;i++)
			tim[i]=i;

	axis_in *train;
	train = (axis_in*)malloc(79184*sizeof(axis_in));
//	axis_in input[157584];
//	train = input;
//	axis_in inp = *train;

	axis_out *spike_count;
//	axis_out output[8];
//	spike_count = output;
	spike_count = (axis_out*)malloc(8*sizeof(axis_out));
//	axis_out out = *spike_count;
//	axis_in *temp = train;

//	(*spike_count).data = 1;

//	printf("Hey!");
	for(int num_e=0;num_e<epoch;num_e++)
	{
		for(int i=1;i<7;i++) //i stands for images
		{
			axis_in *temp = train;
			for(int j=0;j<784;j++)
			{
				for(int k=0;k<101;k++)
				{
					(*temp++).data=0;
				}
			}
			temp=train;

//			printf("Train zero");
			string fname = "100strain"+to_string(i) + ".txt";
			ifstream f(fname);
    		if(f.is_open())
    		{
    			for(int ind_a = 0; ind_a < 784; ind_a++)
        		{
        			for(int ind_b=0;ind_b<101;ind_b++)
        			{
        				f >> (*temp++).data;
//        				printf("Hi!");
        			}
        		}
    		}
    		temp=train;
//    		printf("Train load");
//    		for(int i=0; i<784; i++)
//    		{
//    			for(int j=0; j<201; j++)
//    			{
//    				cout << (*temp++).data << " ";
//    			}
//    			cout << endl;
//    		}
//
			LIF_Network(train, spike_count);
////
			for(int i=0;i<8;i++)
				cout<<(*spike_count++).data<<" ";
			cout<<endl;

		}
	}
}
