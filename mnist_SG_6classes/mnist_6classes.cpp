#include </Users/srinithya/Documents/CSE237D/Spiking-Neural-Network-master/c_code/stdc++.h>
using namespace std;
#include<cstring>


double Pref;
double Prest;
int Pmin;
double Pth;
double D; //decay



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
double** learned_weights();
double dot_done(double* arr1, double* arr2, int len);
double random_uniform(double w_min, double w_max);
//double** encode_stochastic(int** img);
//int** matrix_mul(int** mat1,int** mat2, int mat1_nrows,int mat1_ncols, int mat2_nrows, int mat2_ncols);

void neuron_initialise()
{
	Pref = 0;
	Prest = 0;
	Pmin = -1;
	Pth = 140;
	D = 0.5;
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
	Pth = 150;
	m = 784;
	n = 8;
	epoch = 1;
	num_of_images = 6;
	w_max = 0.5;
	w_min = -0.5;
}

// int** matrix_mul(int** mat1,int** mat2, int mat1_nrows,int mat1_ncols, int mat2_nrows, int mat2_ncols)
// {

// 	int** mult = 0;
// 	mult = new int*[mat1_nrows];
// 	for (int i = 0; i < mat1_nrows; i++)
// 	{
// 		mult[i] = new int[mat2_ncols];

//         for (int j = 0; j < mat2_ncols; j++)
//         {
//                   // fill in some initial values
//                   // (filling in zeros would be more logic, but this is just for the example)
//         	for(int k=0;k<mat1_ncols;k++)
//         	{
//         		mult[i][j] += mat1[i][k]*mat2[k][j];
//         	}
//         }
//     }
//     return mult;
// }
double dot_done(double* arr1, double* arr2, int len)
{
	double mult = 0;
	for (int i = 0; i < len; i++)
	{
		mult+=arr1[i]*arr2[i];
	}
	return mult;
}

double** learned_weights(string filename)
{
	double** weights;
	weights = new double*[6];
	cout<<filename;
	ifstream file(filename);
    if(file.is_open())
    {

        for(int i = 0; i < 6; i++)
        {
        	//cout<<i<<"tadaaaaa";
        	weights[i] = new double[m];
        	for(int j=0;j<m;j++)
        	{
        		file >> weights[i][j];
        		//cout<<"hey hiii"<<endl;
        		//cout<<(weights[i][j])<<" ";
        	}
        }
    }
    return weights;
}

double random_uniform(double w_min, double w_max)
{
	std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> distribution(w_min, w_max);
    double d=distribution(mt);
    //cout<<endl;
    //cout<<d;
    return d;
}


// double** encode_stochastic(int** pot)
// {

// 	double** train;
// 	train = new double*[784];
// 	int n=0;
// 	for(int i=0;i<28;i++)
// 	{
// 		for(int j=0;j<28;j++)
// 		{
// 			int temp[T+1];
// 			int spike_temp[T+1];
// 			train[n] = new double[T+1];
// 			for(int k=0;k<T+1;k++)
// 			{
// 				temp[k]=random_uniform(0, 1.0);
// 				if(temp[k]<pot[i][j])
// 					train[n][k]=1;
// 				else
// 					train[n][k]=0;
// 			}
// 			n=n+1;
// 		}

// 	}
// 	return train;
// }

class neuron
{
	public:
		int t_ref = 4;
		int t_rest = -1;
		double P = 0;
		double Prest = 0;
		double D = 0.5;
		int Pmin = -1;
	int check()
	{
		if(P > Pth)
		{
			P = Pref;
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
		int t_rest = -1;
		double P = Prest;
	}
};



int main()
{
	neuron_initialise();
	param_initialise();
	//cout<<"mmmm"<<endl;
	neuron layer2[n]; //create n number of neuron objs
	double synapse[n][m];
	//cout<<"heyy"<<endl;
	for (int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
			synapse[i][j] =0;
	}
	string fname = "wt.txt";
	//cout<<"heyyy2"<<endl;
	double **weight_matrix = learned_weights(fname);
	//cout<<"heyyy3"<<endl;
	for(int i=0;i<num_of_images;i++)
	{
		for(int j=0;j<m;j++)
		{
			synapse[i][j]=weight_matrix[i][j];
		}
			
	}
	//cout<<"heyyy4"<<endl;
	for(int k=num_of_images;k<n;k++)
	{
		//cout<<endl<<endl<<endl;
		//cout<<"whats ur problem?";
		//cout<<"k is"<<k;
		for(int j=0;j<m;j++)
		{
			synapse[k][j]=random_uniform(w_min,w_max);
		}
	}
	//cout<<endl<<"heyyy5";
	for(int i=0;i<epoch;i++)
	{
		for(int j=1;j<7;j++)
		{
			int spike_count[n];
			double active_pot[n];
			for(int k=0;k<n;k++)
			{
				spike_count[k]=0;
				active_pot[k]=0;
				layer2[k].initial();
			}
			//cout<<"heyyy6";

			//read the image here
			//int** img = 


			// string filename = "st"+to_string(j) + ".txt";
			// //cout<<filename<<endl;
			// int** img;
			// img = new int*[28];
			// ifstream file(filename);
   //  		if(file.is_open())
   //  		{
   //  			for(int i = 0; i < 28; i++)
   //      		{
   //      			img[i] = new int[28];
   //      			for(int j=0;j<28;j++)
   //      			{
   //      				file >> img[i][j];
   //      			}
   //      		}
   //      		//cout<<endl<<"lollol";
   //  		}
			//int** pot=rf(img);
			//double** train = encode_stochastic(img);


			
			double** train;
			train = new double*[784];
			string fname = "strain"+to_string(j) + ".txt";
			ifstream f(fname);
    		if(f.is_open())
    		{
    			for(int i = 0; i < 784; i++)
        		{
        			train[i] = new double[201];
        			for(int j=0;j<201;j++)
        			{
        				f >> train[i][j];
        			}
        		}
        		//cout<<endl<<"lollol again"<<endl;
    		}
    		//cout<<endl<<"just checking"<<train[547][90]<<endl;
			int f_spike=0;
			for(int t=0;t<T+1;t++)
			{
				for(int j=0;j<n;j++)
				{
					if(layer2[j].t_rest<t)
					{
						double* train_temp = 0;
						train_temp = new double[m];
						for(int ne=0;ne<m;ne++)
						{
							train_temp[ne]=train[ne][t];
						}
						double *syntemp;
						syntemp = new double[m];
						for(int i=0;i<m;i++)
							syntemp[i] = synapse[j][i];
						layer2[j].P += dot_done(syntemp,train_temp,T+1); //np.dot(synapse[j], train[:,t])
						//cout<<"dot prod"<<layer2[j].P<<endl;
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
					

					if(high_pot > Pth)
					{
						cout<<"high_pot"<<high_pot<<endl;
						f_spike=1;
						const int N = sizeof(active_pot) / sizeof(int);
						double winner = distance(active_pot, max_element(active_pot, active_pot + N));
						cout<<"winner isss"<<" "<<winner<<endl;
						for(int i=0;i<n;i++)
						{
							if(i!=winner)
							{
								layer2[i].P = layer2[i].Pmin;
							}
						}

					}
				}
				for(int j=0;j<n;j++)
				{
					int s= layer2[j].check();
					if(s==1)
					{
						spike_count[j] += 1;
						//cout<<endl;
						//cout<<"spikecount"<<spike_count[j]<<endl;
						layer2[j].t_rest = t + layer2[j].t_ref;
					}
				}
			}

			for(int i=0;i<n;i++)
				cout<<spike_count[i]<<" ";
			cout<<endl;

		}
	}
}





