#include "lif.hpp"
#include <bits/stdc++.h>
using namespace std;

int layers[3] = {784, 64, 10};
int T=30;
float pot_rest = 0;
float pot_thr = 1;
float tau_ref = 0.002;

float nodes0[784] = {0};
float nodes1[64] = {0};
float nodes2[10] = {0};

float refractory_time0[784] = {0.002};
float refractory_time1[64] = {0.002};
float refractory_time2[10] = {0.002};

float expm1f(float x){
	return expf(x)-1;
}

float log1pf(float x){
       return logf(1+x);
}

float clip(float n, float lower, float upper){
	return std::max(lower, std::min(n, upper));
}

void m_multiply1(float train_new[64], float train[784])
{
	float res[64];
	for(int i=0; i < 64; i++){
		res[i] = 0;
		for(int j=0; j<784; j++)
		{
#pragma HLS PIPELINE
			res[i] += train[j] * weights1[i][j];
		}
	}
	for(int i=0; i<64; i++)
	{
		train_new[i] = res[i];
	}
}

void m_multiply2(float train_new2[10], float train[64])
{
	float res[10];
	for(int i=0; i < 10; i++){
		res[i] = 0;
		for(int j=0; j<64; j++)
		{
#pragma HLS PIPELINE
			res[i] += train[j] * weights2[i][j];
		}
	}
	for(int i=0; i<64; i++)
	{
		train_new2[i] = res[i];
	}
}

void update_voltage(float* train, float* nodes, int n, float t, float* refractory_time)
{
	float tau_rc = 0.02;
    float tau_ref = 0.005;
    float delta,nu;
    float spike_mask;
    float t_spike;
//	float temp[n];

//	for(int i=0;i<n;i++)
//		temp[i]=nodes[i];

	for(int i=0; i < n; i++){
		refractory_time[i] -= t;
        nu = t - refractory_time[i];
        delta = clip(nu,0,t);
		nodes[i] -= (train[i]-nodes[i])*expm1f(-delta/tau_rc);
                //cout << t<< "h" << expm1(-t/tau_rc) << " ";
		if(nodes[i] > 1) train[i] = 1/t;
        //spike_mask = temp[i] > 1;
        float eps = pow(10, -9);
        //t_spike = t + tau_rc * log1p(-1*(temp[i]-1)/(train[i]-1)+eps);
        t_spike = t + tau_rc * log1pf(-1*(nodes[i]-1)/(train[i]-1));
        if(nodes[i] < 0) nodes[i] =0;

        refractory_time[i] = tau_ref + t_spike;
	}

//	for(int i=0; i < n; i++){
//		nodes[i] = temp[i];
//	}
}

void LIF_Network(axis_in *inp, axis_out *spike_sum)
{
#pragma HLS INTERFACE axis depth=10 port=spike_sum
#pragma HLS INTERFACE axis depth=23520 port=inp

	int spikes[30][10];
	float train0[784];
#pragma HLS ARRAY_PARTITION variable=train0 cyclic factor=4 dim=1
	float train1[64];
#pragma HLS ARRAY_PARTITION variable=train1 cyclic factor=4 dim=1
	float train2[10];
//	float inp_internal[30][784];

//	for(int i=0; i<30; i++)
//	{
//		for(int j=0; j<784; j++)
//		{
////			inp_internal[i][j] = (*inp++).data;
//			cout << (*inp++).data << " ";
//		}
//		cout << endl;
//	}

	for(int t=0; t < 30; t++)
	{
		for(int ind5=0;ind5<784;ind5++)
		{
//			train0[ind5]=inp_internal[t][ind5];
			train0[ind5] = (*inp++).data;
//			cout << train0[ind5] << " ";
		}
//		cout <<endl;

		update_voltage(train0, nodes0, 784, dt[t], refractory_time0);

		for(int sp=0; sp < 784; sp++)
			if(nodes0[sp] > pot_thr)
				nodes0[sp] = pot_rest;

		m_multiply1(train1, train0);
//		for(int i=0; i<784; i++)
//		{
//			cout << train1[i] << " ";
//		}
//		cout << endl;
		update_voltage(train1, nodes1, 64, dt[t], refractory_time1);
		for(int j=0; j < 64; j++)
			if(nodes1[j] > pot_thr)
				nodes1[j] = pot_rest;

//		for(int i=0; i<784; i++)
//		{
//			cout << nodes2[i] << " ";
//		}
//		cout << endl;

		m_multiply2(train2, train1);
		update_voltage(train2, nodes2, 10, dt[t], refractory_time2);

//		for(int i=0; i<784; i++)
//		{
//			cout << nodes2[i] << " ";
//		}
//		cout << endl;

		for(int k=0; k < 10; k++)
		{
			if(nodes2[k] > pot_thr)
			{
				nodes2[k] = pot_rest;
				spikes[t][k] = 1;
			}
			else
				spikes[t][k] = 0;
		}
	}

//	for(int i=0; i<30; i++)
//	{
//		for(int j=0; j<10; j++)
//		{
//			cout << spikes[i][j] << " ";
//		}
//		cout << endl;
//	}

	for(int sn=0; sn < 10; sn++)
	{
		int tot = 0;
		for(int sm=0; sm < 30; sm++)
			tot += spikes[sm][sn];

		axis_out out;
		out.last=0;
		out.data = tot;
		if(sn==9)
			out.last=1;
		*spike_sum++ = out;
	}
}
