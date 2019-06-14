#include </Users/srinithya/Documents/CSE237D/c_code/stdc++.h>
using namespace std;
#include<cstring>


double Pref_i;
double Prest_i;
double Pmin_i;
double Pth_i;
double D_i; //decay



int T;
int* tim;//
int t_back;
int t_fore;
int m;
int n;
int epoch;
int num_of_images;
double w_max;
double w_min;


void neuron_initialise();
void param_initialise();
double dot_done(double (&arr1)[784], double (&arr2)[784], int len);
void learned_weights(string filename, double (&weight_m)[8][784]);

void neuron_initialise()
{
	Pref_i = 0;
	Prest_i = 0;
	Pmin_i = -1;
	Pth_i = 150;
	D_i = 0.5;
	//cout<<"lol";
}

void param_initialise()
{
	T = 200;
	tim = new int[T+1];
	for(int i=1;i<T+1;i++)
		tim[i]=i;
	t_back = -20;
	t_fore = 20;
	m = 784;
	n = 8;
	epoch = 1;
	num_of_images = 6;
	w_max = 0.5;
	w_min = -0.5;
}


double dot_done(double (&arr1)[784], double (&arr2)[784], int len)
{
	// cout<<"ready"<<endl;
	// for(int ind_f=0;ind_f<len;ind_f++)
	// 	cout<<arr2[ind_f]<<" ";
	// cout<<endl;
	double mult = 0;
	for (int i = 0; i < len; i++)
	{
		mult+=arr1[i]*arr2[i];
	}
	return mult;
}

void learned_weights(string filename, double (&weight_m)[8][784])
{
	//cout<<filename;
	ifstream file(filename);
    if(file.is_open())
    {

        for(int i = 0; i < n; i++)
        {
        	//cout<<i<<"tadaaaaa";
        	for(int j=0;j<m;j++)
        	{
        		file >> weight_m[i][j];
        		//cout<<"hey hiii"<<endl;
        		//cout<<(weights[i][j])<<" ";
        	}
        }
    }
}


class neuron
{
	public:
		int n_pth = Pth_i;
		int t_ref = 4;
		int t_rest = -1;
		double P = 0;
		double Prest = 0;
		double D = 0.5;
		double Pmin = -1;
	int check()
	{
		if(P >= n_pth)
		{
			P = Pref_i;
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



int main()
{
	neuron_initialise();
	param_initialise();
	neuron layer2[n]; //create n number of neuron objs
	double synapse[n][m];
	for (int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
			synapse[i][j] =0;
	}
	string fname = "weights.txt";
	double weight_matrix[8][784];
	learned_weights(fname,weight_matrix);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			synapse[i][j]=weight_matrix[i][j];
		}
			
	}

	for(int num_e=0;num_e<epoch;num_e++)
	{
		for(int i=1;i<7;i++) //i stands for images
		{
			int spike_count[n];
			double active_pot[n];
			//intialise output spikes and potentials of neurons
			for(int k=0;k<n;k++)
			{
				spike_count[k]=0;
				active_pot[k]=0;
				layer2[k].initial();
				//cout<<" initial potential and time rest are "<<layer2[k].P<<"&"<<" "<<layer2[k].t_rest<<endl;
				//cout<< "Initial pot of "<<k<<"is "<<layer2[k].P<<endl;
			}
			cout<<endl;
			//train load
			double train[784][201];
			string fname = "strain"+to_string(i) + ".txt";
			ifstream f(fname);
    		if(f.is_open())
    		{
    			for(int ind_a = 0; ind_a < 784; ind_a++)
        		{
        			for(int ind_b=0;ind_b<201;ind_b++)
        			{
        				f >> train[ind_a][ind_b];
        			}
        		}
    		}
    		//train load end 
			int f_spike=0;
			for(int t=1;t<T+1;t++)
			{
				for(int j=0;j<n;j++)
				{
					if(layer2[j].t_rest<t)
					{
						//cout<<"yes i am less than 0"<<" and my j is "<<j<<"my t is "<<t<<endl;
						double train_temp[784];
						for(int ne=0;ne<m;ne++)
						{
							train_temp[ne]=train[ne][t];
							// if(t==1&&j==1)
							// {
							// 	cout<< train_temp[ne]<< " ";
							// }
						}
						//cout<<endl<<"end of train"<<endl;
						double syntemp[784];
						//cout<<"t is "<<t<<"and j is "<<j<<endl;
						for(int ind_c=0;ind_c<m;ind_c++)
						{
							syntemp[ind_c] = synapse[j][ind_c];
							// if(t==1 && j==1)
							// 	cout<<syntemp[ind_c]<<" ";
						}
						double dot_buf = dot_done(syntemp,train_temp,784);
						layer2[j].P += dot_buf; //np.dot(synapse[j], train[:,t])
						// if(t ==1 && j==1)
						// 	cout<<"dot prod"<<dot_buf<<endl;
						if(layer2[j].P>layer2[j].Prest)
							layer2[j].P -= layer2[j].D;
						active_pot[j] = layer2[j].P;
					}
				}
				if(f_spike==0)
				{
					double high_pot;
					double* highpot_ind;
					double* poi = active_pot;
					highpot_ind = max_element(active_pot,active_pot+n);
					long na;
					na = highpot_ind-poi;
					//cout<<"heyyy"<<endl;
					//cout<<na<<"index"<<endl;
					high_pot = double(*highpot_ind);
					

					if(high_pot > Pth_i)
					{
						// cout<<"high_pot"<<high_pot<<endl;
						f_spike=1;
						const int N = sizeof(active_pot) / sizeof(int);
						double winner = distance(active_pot, max_element(active_pot, active_pot + N));
						// cout<<"winner isss"<<" "<<winner<<endl;
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
						spike_count[ind_e] += 1;
						//cout<<endl;
						//cout<<"spikecount"<<spike_count[j]<<endl;
						layer2[ind_e].t_rest = t + layer2[ind_e].t_ref;
					}
				}
			}

			for(int i=0;i<n;i++)
				cout<<spike_count[i]<<" ";
			cout<<endl;

		}
	}
	cout<<endl;
}





