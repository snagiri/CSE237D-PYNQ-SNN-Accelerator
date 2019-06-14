#include "lif.hpp"
#include <bits/stdc++.h>
using namespace std;
#include<cstring>

int m=784;
int n=8;

#define aptype float

class neuron
{
	public:
		int n_pth = 150;
		int t_ref = 4;
		int t_rest = -1;
		float P = 0;
		float Prest = 0;
		float D = 0.5;
		float Pmin = -1;
	int check()
	{
		if(P >= n_pth)
		{
			P = 0;
			//cout<<"lala"<<endl;
			return 1;
		}
		else if (P < Pmin)
		{
			P = Prest;
			return 0;
		}
		else
			return 0;

	}
	void inhibit()
	{
		P = Pmin;
	}
	void initial()
	{
		t_rest = -1;
		P = Prest;
	}
};

float dot_done(float (&arr1)[784], float (&arr2)[784], int len)
{
	float mult = 0;
	for (int i = 0; i < len; i++)
	{
		mult+=arr1[i]*arr2[i];
	}
	return mult;
}

void LIF_Network(axis_in *train, axis_out *spike_count)
{
#pragma HLS INTERFACE axis depth=8 port=spike_count
#pragma HLS INTERFACE axis depth=79184 port=train
	int f_spike=0;
	neuron layer2[8];
	aptype active_pot[8];
	int temp_spike[8];
	float temp_train_complete[101][784];

	for(int i=0; i<8; i++)
		temp_spike[i] = 0;

	for(int k=0;k<8;k++)
	{
		active_pot[k]=0;
		layer2[k].initial();
	}

//	axis_in *temp_ptr = train;
//	axis_in *temp_next;

	for(int i=0; i< 784; i++)
	{
		for(int j=0; j<101; j++)
		{
			temp_train_complete[j][i] = (*train++).data;
		}
	}

	for(int t=1;t<101;t++)
	{
//#pragma HLS UNROLL factor=10
//		temp_ptr = train;
		for(int j=0;j<n;j++)
		{
			if(layer2[j].t_rest<t)
			{
				//cout<<"yes i am less than 0"<<" and my j is "<<j<<"my t is "<<t<<endl;
				float train_temp[784];
				for(int ne=0;ne<m;ne++)
				{
//					train_temp[ne]=(*temp_ptr).data;
					train_temp[ne]=temp_train_complete[t][ne];
//					 if(t==3&&j==1)
//					 {
//					 	cout<< train_temp[ne]<< " ";
//					 }
//					 temp_ptr = temp_ptr + 101;
				}
//				temp_ptr = train;
//				train++;
				//cout<<endl<<"end of train"<<endl;
				float syntemp[784];
				//cout<<"t is "<<t<<"and j is "<<j<<endl;
				for(int ind_c=0;ind_c<m;ind_c++)
				{
					syntemp[ind_c] = weight_matrix[j][ind_c];
					// if(t==1 && j==1)
					// 	cout<<syntemp[ind_c]<<" ";
				}
				float dot_buf = dot_done(syntemp,train_temp,784);
				layer2[j].P += dot_buf; //np.dot(synapse[j], train[:,t])
//				if(t ==1 && j==1)
//						cout<<"dot prod"<<dot_buf<<endl;
				if(layer2[j].P>layer2[j].Prest)
					layer2[j].P -= layer2[j].D;
				active_pot[j] = layer2[j].P;
			}
		}
		if(f_spike==0)
		{
			float high_pot;
			aptype* highpot_ind;
			aptype* poi = active_pot;
			highpot_ind = max_element(active_pot,active_pot+n);
			long na;
			na = highpot_ind-poi;
			//cout<<"heyyy"<<endl;
			//cout<<na<<"index"<<endl;
			high_pot = float(*highpot_ind);

			if(high_pot > 150)
			{
//					cout<<"high_pot"<<high_pot<<endl;
				f_spike=1;
				const int N = sizeof(active_pot) / sizeof(int);
				float winner = distance(active_pot, max_element(active_pot, active_pot + N));
//					cout<<"winner isss"<<" "<<winner<<endl;
				for(int ind_d=0;ind_d<n;ind_d++)
				{
					if(ind_d!=winner)
					{
						layer2[ind_d].P = layer2[ind_d].Pmin;
					}
				}

			}
		}

		for(int ind_e=0;ind_e<n;ind_e++)
		{
			int s= layer2[ind_e].check();
			if(s==1)
			{
				temp_spike[ind_e] += 1;
				//cout<<endl;
				//cout<<"spikecount"<<spike_count[j]<<endl;
				layer2[ind_e].t_rest = t + layer2[ind_e].t_ref;
			}
		}
//		train++;

//		for(int i=0; i<8; i++)
//		{
//			cout << active_pot[i] << " ";
//		}
//		cout << endl;
	}
	for(int i=0; i<8; i++)
	{
//		(*spike_count++).data = temp_spike[i];
		axis_out out;
		out.last = 0;
		out.data = temp_spike[i];
		if(i == 7)
			out.last = 1;
//		cout<<out.last<<" ";
		*spike_count++ = out;
	}
//	spike_count = temp;
}

