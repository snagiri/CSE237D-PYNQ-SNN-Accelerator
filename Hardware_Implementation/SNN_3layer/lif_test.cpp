#include <bits/stdc++.h>
using namespace std;
#include "lif.hpp"

int main(){

	axis_in *inp;
	inp = (axis_in*) malloc(23520*sizeof(axis_in));

	axis_out *spike_sum;
	spike_sum = (axis_out*) malloc(10*sizeof(axis_out));

	int epochs = 1, n_images = 2;
	for(int i=0; i < epochs; i++)
	{
		for(int j=1; j < n_images; j++)
		{
			axis_in *temp = inp;
			for(int j=0;j<30;j++)
			{
				for(int k=0;k<784;k++)
				{
					(*temp++).data=0;
				}
			}
			temp=inp;

			string fname = "st"+to_string(j) + ".txt";
			ifstream f(fname);
			if(f.is_open())
			{
				for(int ind_a = 0; ind_a < 30; ind_a++)
				{
					for(int ind_b=0;ind_b<784;ind_b++)
					{
						f >> (*temp++).data;
					}
				}
			}
			temp = inp;

			LIF_Network(inp, spike_sum);

			for(int i=0;i<10;i++)
				cout<<(*spike_sum++).data<<" ";
			cout<<endl;
		}
	}
	return 0;
}
