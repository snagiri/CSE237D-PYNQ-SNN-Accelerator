#include </Users/srinithya/Documents/CSE237D/stdc++.h>
using namespace std;

void learned_weights(string filename1, string filename2, double (&weight_1)[64][784], double (&weight_2)[10][64])
{
	cout<<filename1;
	ifstream file(filename1);
    if(file.is_open())
    {

        for(int i = 0; i < 64; i++)
        {
        	//cout<<i<<"tadaaaaa";
        	for(int j=0;j<784;j++)
        	{
        		file >> weight_1[i][j];
        		//cout<<"hey hiii"<<endl;
        		//cout<<(weights[i][j])<<" ";
        	}
        }
    }
    file.close();
    cout<<filename2;
    ifstream f(filename2);
    if(f.is_open())
    {

        for(int i = 0; i < 10; i++)
        {
        	//cout<<i<<"tadaaaaa";
        	for(int j=0;j<64;j++)
        	{
        		f >> weight_2[i][j];
        		//cout<<"hey hiii"<<endl;
        		//cout<<(weights[i][j])<<" ";
        	}
        }
    }
    f.close();
}

void read_spike_train(string fname, double (&inp)[30][784])
{
	ifstream fi(fname);
	if(fi.is_open())
	{
    	for(int ind_a = 0; ind_a <30; ind_a++)
        {
        	for(int ind_b=0;ind_b<784;ind_b++)
        		fi >> inp[ind_a][ind_b];
    	}
    }
    fi.close();
}

void m_multiply1(double (&train_new)[64], double train[784], double w[64][784])
{
	double res[64]={0};

	//for(int i=0;i<784;i++)
	//	cout<<res[i]<<" ";
	for(int i=0; i < 64; i++){
		res[i] = 0;
		for(int j=0; j<784; j++)
		{
			res[i] += train[j] * w[i][j];
		}
	}

	for(int k=0; k<64; k++)
	{
		train_new[k] = res[k];
		//cout<<train[j]<<" ";
	}
}
void m_multiply2(double (&train_new2)[10],double train[64], double w2[10][64])
{
	double res2[10]={0};

	//for(int i=0;i<784;i++)
	//	cout<<res[i]<<" ";
	for(int i=0; i < 10; i++){
		res2[i] = 0;
		for(int j=0; j<64; j++)
		{
			res2[i] += train[j] * w2[i][j];
		}
	}

	for(int k=0; k<10; k++)
	{
		train_new2[k] = res2[k];
		//cout<<train[j]<<" ";
	}
}

double expm1(double x){
	return exp(x)-1;
}
double log1p(double x){
       return log(1+x);
}
       
double clip(double n, double lower, double upper){
	return std::max(lower, std::min(n, upper));
}

void update_voltage(double* train, double* nodes, int n,double t, double* refractory_time)
{
	double tau_rc = 0.02;
    double tau_ref = 0.005;
    double delta,nu;
    double spike_mask;
    double t_spike;
	double temp[n];
	for(int i=0;i<n;i++)
		temp[i]=nodes[i];

	for(int i=0; i < n; i++){
		refractory_time[i] -= t;
        nu = t - refractory_time[i];
        delta = clip(nu,0,t);
		temp[i] -= (train[i]-temp[i])*expm1(-delta/tau_rc);
                //cout << t<< "h" << expm1(-t/tau_rc) << " ";
		if(temp[i] > 1) train[i] = 1/t;
        //spike_mask = temp[i] > 1;
        double eps = pow(10, -9);
        //t_spike = t + tau_rc * log1p(-1*(temp[i]-1)/(train[i]-1)+eps);
        t_spike = t + tau_rc * log1p(-1*(temp[i]-1)/(train[i]-1));
        if(temp[i] < 0) temp[i] =0;

        refractory_time[i] = tau_ref + t_spike;
            
	}
	for(int i=0; i < n; i++){
		nodes[i] = temp[i];
	}
        
}

int main(){
	int layers[3] = {784, 64, 10};
	double pot_rest = 0, pot_thr = 1, tau_ref = 0.002;;


	string wt_files[2] = {"wt1.txt", "wt2.txt"};

	double weights1[64][784];
	double weights2[10][64];

	learned_weights(wt_files[0], wt_files[1], weights1, weights2);
	//cout<<(weights.size())<<"weightsss"<<endl;
	//cout<<(weights[0].size())<<" "<<weights[1].size()<<" "<<weights[2].size()<<" "<<weights[3].size()<<" "<<endl;
	//cout<<(weights[0][1].size())<<endl;
	int T=30;
	double dt[T];
	for(int ind1=0;ind1<30;ind1++)
		dt[ind1] = 0.001;

	for(int t=1; t < 30; t++) dt[t] += dt[t-1];
	
	int epochs = 1, n_images = 3;
	for(int i=0; i < epochs; i++){
		for(int j=2; j < n_images; j++)
		{
			vector< vector<double> > nodes;
			double nodes0[784] = {0};
			double nodes1[64] = {0};
			double nodes2[10] = {0};

			double refractory_time0[784] = {0.002};
			double refractory_time1[64] = {0.002};
			double refractory_time2[10] = {0.002};

			// for(int ind2=0;ind2<784;ind2++)
			// {
			// 	nodes0[ind2]=pot_rest;
			// 	refractory_time0[ind2]=tau_ref;
			// }

			// for(int ind3=0;ind3<64;ind3++)
			// {
			// 	nodes1[ind3]=pot_rest;
			// 	refractory_time1[ind3]=tau_ref;
			// }

			// for(int ind4=0;ind4<10;ind4++)
			// {
			// 	nodes2[ind4]=pot_rest;
			// 	refractory_time2[ind4]=tau_ref;
			// }

			int spikes[30][10];

			string fname = "st"+ to_string(j) + ".txt";
			double inp[30][784];
			read_spike_train(fname, inp);
			
			for(int t=0; t < T; t++)
			{
				cout<<"t"<<" "<<"value is"<<" "<<t<<endl<<endl;
				double train0[784];
				double train1[64];
				double train2[10];

				for(int ind5=0;ind5<784;ind5++)
				{
					train0[ind5]=inp[t][ind5];
					cout<<nodes0[ind5]<<" ";
				}
				cout<<endl<<"node0"<<endl;
				
				update_voltage(train0, nodes0, 784, dt[t], refractory_time0);

				for(int sp=0; sp < 784; sp++)
					if(nodes0[sp] > pot_thr)
						nodes0[sp] = pot_rest;

				m_multiply1(train1, train0, weights1);
				update_voltage(train1, nodes1, 64, dt[t], refractory_time1);
				for(int j=0; j < 64; j++)
					if(nodes1[j] > pot_thr)
						nodes1[j] = pot_rest;

				for(int ind6=0;ind6<64;ind6++)
					cout<<nodes1[ind6]<<" ";
				cout<<endl<<" node1"<<endl;
				m_multiply2(train2, train1, weights2);
				update_voltage(train2, nodes2, 10, dt[t], refractory_time2);

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
				for(int ind7=0;ind7<10;ind7++)
					cout<<nodes2[ind7]<<" ";
				cout<<endl<<" nodes2"<<endl;


			}	
                      
            int spike_sum[10];
			

			for(int sn=0; sn < 10; sn++)
			{
				int tot = 0;
				for(int sm=0; sm < 30; sm++)
					tot += spikes[sm][sn];
				
				spike_sum[sn] = tot;	
			}

			cout<<endl<<endl<<endl;
			cout<<"yayyayyuy";

			for(int num=0;num<10;num++)
				cout << spike_sum[num] << " ";
			cout << endl;

		}
	}
	return 0;
}
