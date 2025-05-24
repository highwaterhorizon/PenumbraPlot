#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Core>
#include "eigen-3.4.0/Eigen/Dense"
#include <iomanip>
#include <ctime>
#include <sstream>

#include "flow.h"
#include "ray_methods.h"

#define G 6.6743015e-11
#define M_e 5.9722e24
#define PI 3.1415926535897
#define r_e 6378e3
#define R_e 149.60e9
#define M_s 1.98892e30
#define R_m 384.1e6
#define M_m 7.346e22
#define M_j 1.89813e27
#define R_j 778.479e9
#define r_s 6.95700e8
#define r_m 1737.4e3


using namespace std;
using namespace std::chrono;
using namespace Eigen;

double r_pow_3(double x,double y) {
	double output,inter;
	inter = pow(x,2)+pow(y,2);
	output = pow(inter,3./2.);
		
	return output;
}


int main(void){
	
	//get config values from config file:
	ifstream conf_ifstream("./main.conf");
	unordered_map<string, string> config;
	string line;
	while (getline(conf_ifstream, line)) {
        // Trim leading whitespaces
        line.erase(0, line.find_first_not_of(" \t"));

        // Skip empty lines or comment lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        string key, value;
        
        // Find the delimiter '=' and split
        if (getline(ss, key, '=') && getline(ss, value)) {
            config[key] = value;
        }
    }
	init_year = stoi(config["init_year"]);
	init_month = stoi(config["init_month"]);
	init_day = stoi(config["init_day"]);
	init_hour = stoi(config["init_hour"]);
	init_min = stoi(config["init_min"]);
	init_second = stoi(config["init_second"]);
	double phi_c_multiplier = stod(config["phi_c_multiplier"]);
	int umbra_res = stoi(config["umbra_res"]);
	int penumbra_res = stoi(config["penumbra_res"]);
	
    Eigen::IOFormat TabFormat(15, Eigen::DontAlignCols, "\t", "\n");
	Eigen::IOFormat CSVFormat(5, Eigen::DontAlignCols, ",", "\n");

	double t_grid;
	int linecount;
	
	int days;
	string input;

	Vector3d omega;
	omega << 1.54277026861935e1, 2.528775803839025e3, 5.832096956033017e3;

	cout << endl << "Reference date is "<< init_day << "." << init_month << "." << init_year << " " << setw(2) << setfill('0') << init_hour << ":" << setw(2) << setfill('0') << init_min << ":" << setw(2) << setfill('0') << init_second << "UTC.\nNumber of days to calculate (integer): [90] ";
	getline(cin, input);
	if (input.empty()) {
		days = 90;
	}
	else {
		days = stoi(input);
	}


	double t_end = days*24*3600;
	//Benchamrk int N = 3e6;	
	long int N; // number of steps
	char input_N;

	cout << "Number of Steps to use (integer): [1e7] ";
	getline(cin, input);

    if (input.empty()) {
        N = 999999;  // Set default value if nothing is entered
    }
	else {
		N = stoi(input);
	}


	/*------------------------------------------------------------------------------------------------------------*/
	/*------------------------------------------------------------------------------------------------------------*/
	/*--------------------------------------------CALCULATION-BEGINS----------------------------------------------*/
	/*------------------------------------------------------------------------------------------------------------*/
	/*------------------------------------------------------------------------------------------------------------*/

	double t_delta = t_end/N;

	//Arrays for velocity an position with initial conditions
	VectorXd r_init(12+1);
	VectorXd v_init(12+1);
	ArrayXXd r(2,12+1);
	ArrayXXd v(2,12+1);

    r_init(0) = 0;
    v_init(0) = 0;
	
	//Sun
	//v_init.segment(1,3) << 1.082108190818754E-02, -9.364414778449289E-03, -1.544564609748103E-04;
	//r_init.segment(1,3) << -1.031800950813197E+06, -6.152204568066536E+05,  2.950524328872960E+04;
	
	v_init.segment(1,3)(0) = stod(config["vs_x"])*1E3;
	v_init.segment(1,3)(1) = stod(config["vs_y"])*1E3;
	v_init.segment(1,3)(2) = stod(config["vs_z"])*1E3;

	r_init.segment(1,3)(0) = stod(config["rs_x"])*1E3;
	r_init.segment(1,3)(1) = stod(config["rs_y"])*1E3;  
	r_init.segment(1,3)(2) = stod(config["rs_z"])*1E3;    
	


	//Earth
	/*r_init.segment(4,3) << 5.052379418322767e10, -1.436902867759990e11,3.715806875360012e7;
	v_init.segment(4,3) << 2.755195392560378e4, 9.987354520126058e3,1.462542185284299e-01;*/
	
	v_init.segment(4,3)(0) = stod(config["ve_x"])*1E3;
	v_init.segment(4,3)(1) = stod(config["ve_y"])*1E3;
	v_init.segment(4,3)(2) = stod(config["ve_z"])*1E3;

	r_init.segment(4,3)(0) = stod(config["re_x"])*1E3;
	r_init.segment(4,3)(1) = stod(config["re_y"])*1E3;  
	r_init.segment(4,3)(2) = stod(config["re_z"])*1E3;


	//Moon
	/*
    r_init.segment(7,3) << 5.011968031234255e10,-1.436804984379607e11,4.401479551874101e7;
	v_init.segment(7,3) << 2.752236630992519e4,9.023592869797513e3,-8.341260037563725e1;*/
	
	v_init.segment(7,3)(0) = stod(config["vm_x"])*1E3;
	v_init.segment(7,3)(1) = stod(config["vm_y"])*1E3;
	v_init.segment(7,3)(2) = stod(config["vm_z"])*1E3;

	r_init.segment(7,3)(0) = stod(config["rm_x"])*1E3;
	r_init.segment(7,3)(1) = stod(config["rm_y"])*1E3;  
	r_init.segment(7,3)(2) = stod(config["rm_z"])*1E3;


	//Jupiter
	/*
	v_init.segment(10,3) << -1.177372976659906e4,6.582793582172885e3,2.360518944456147e2;	
	r_init.segment(10,3) << 3.428615675815491e11,6.686251492898077e11,-1.044446430752873e10;*/ 
	
	v_init.segment(10,3)(0) = stod(config["vj_x"])*1E3;
	v_init.segment(10,3)(1) = stod(config["vj_y"])*1E3;
	v_init.segment(10,3)(2) = stod(config["vj_z"])*1E3;

	r_init.segment(10,3)(0) = stod(config["rj_x"])*1E3;
	r_init.segment(10,3)(1) = stod(config["rj_y"])*1E3; 
	r_init.segment(10,3)(2) = stod(config["rj_z"])*1E3;

    r.row(0) = r_init;
    v.row(0) = v_init;

    ofstream pos_ems_ofstream;
    pos_ems_ofstream.open("./data/pos_ems.dat");


	//Initialising progress bar
	float progress = 0.0;
	int barWidth = 20;
	int pos, spacebar;
	cout << endl;
	
	
	//flow array
	Array<double,1,26> flow;

	auto start = high_resolution_clock::now();
	

	//Initialising Runge-Kutta-4
	Array<double,1,26> k1;
	Array<double,1,26> k2;
	Array<double,1,26> k3;
	Array<double,1,26> k4;


	//Intialising recording parameters

	Vector3d r_ms,r_es;
	double d_ms,d_es,phi,phi_c;
	linecount = 0;
	for (int i = 0; i<N; i++) {

        k1 = full_flow_f(r.row(0),v.row(0));
		k2 = full_flow_f(r.row(0)+0.5*t_delta*k1.segment(13,13),v.row(0)+0.5*t_delta*k1.segment(0,13));
		k3 = full_flow_f(r.row(0)+0.5*t_delta*k2.segment(13,13),v.row(0)+0.5*t_delta*k2.segment(0,13));
		k4 = full_flow_f(r.row(0)+    t_delta*k3.segment(13,13),v.row(0)+    t_delta*k3.segment(0,13));

		//TODO optimization: v is added six times and sum is divided by six this can be ignored and only the 			relevant part should be added, i.e. the acceleration components
		flow = (1/6.)*(k1+2*k2+2*k3+k4);

		v.row(1) = v.row(0) + t_delta * flow.segment(0,13);
		r.row(1) = r.row(0) + t_delta * flow.segment(13,13);

		r.row(1)(0) = (i+1)*t_delta;
		v.row(1)(0) = (i+1)*t_delta;

        //pos_ems_file<<i*t_delta<<"\t"<<r.row(1).segment(4,3)<<"\t"<<r.row(1).segment(7,3)<<"\t"<<r.row(1).segment(1,3)<<"\n";

		

        v.row(0) = v.row(1);
        r.row(0) = r.row(1);


		//Start recording when critical angle is reached:

		r_ms = r.row(1).segment(7,3)-r.row(1).segment(1,3);
		r_es = r.row(1).segment(4,3)-r.row(1).segment(1,3);
		d_ms = r_ms.norm();
		d_es = r_es.norm();
		phi_c = phi_c_multiplier * atan(r_e/d_es);
		phi = acos((r_ms.dot(r_es)) / (r_ms.norm()*r_es.norm()) );
		if (d_ms < d_es && phi < phi_c) {
			
			pos_ems_ofstream << r.row(1).format(TabFormat) << endl;
			linecount +=1;

		}

		




		//Progress bar
		spacebar = 0;
		if (i % (int(N)/20) == 0) 
		{
		progress += 0.05;
		if (progress > 1)
		{
			progress = 1;
		}
		if (progress < 0.1) 
		{
			spacebar = 2; 
		}
		else if (progress > 0.05 && progress < 1)
		{
			spacebar = 1;
		}
		else 
		{
			spacebar = 0;
		}
		cout << "Integrating differential equation" << setw(6+spacebar) << setfill(' ') << " \033[32m" << progress*100.0 << "%" << setw(spacebar) << "\033[0m" << " [";
		pos = barWidth * progress;
		for (int j = 0; j < barWidth;  j++) {
			if (j < pos) cout << "=";
			else if (j == pos) cout << ">";
			else cout << " ";
		}
		cout << "] " << " (" << int(progress * days) << "/" << days << " days)\r";
		cout.flush();
		}
		

	}

	auto stop = high_resolution_clock::now();
	
	//Resetting progress bar
	progress = 0.0;
	cout << endl;

	pos_ems_ofstream.close();
	//Creating latlong file for recorded values:

	//1st Step: Calculating Tangent Parameters for every recorded set of positions
	//ifstream pos_ems_ifstream;


	ifstream pos_ems_ifstream("./data/pos_ems.dat");
	ofstream param_ems_ofstream("./data/param_ems.dat");
	ofstream umbraisx_ems_ofstream("./data/umbraisx_ems.dat");
	ofstream latlong_ofstream("./data/latlong.dat");
	ofstream param_ems_in_ofstream("./data/param_ems_in.dat");
	ofstream umbraisx_ems_in_ofstream("./data/umbra_ems_in.dat");
	ofstream latlong_in_ofstream("./data/latlong_in.dat");

	//Initialising progress bar
	progress = 0.0;




	stringstream ss;
	int while_idx = 0, vl_idx;
	ArrayXXd ps(3,umbra_res), us(3,umbra_res), earth_isx(3,umbra_res);
	ArrayXXd ps_in(3,penumbra_res), us_in(3,penumbra_res), earth_isx_in(3,penumbra_res);
	string value;

	while (getline(pos_ems_ifstream,line)) {
		while_idx +=1;
		vl_idx = 0;
		ss.clear();
		ss.str(line);
		//cout << ss.str() << endl;
        while(getline(ss,value,'\t')) {
			//reusing r array
			
            r.row(0)(vl_idx) = stod(value);
			//cout << value << endl;
			vl_idx += 1;               
            }
		//cout << r.row(0).segment(7,3) << endl;
		ps = umbrapar(r.row(0).segment(7,3),r.row(0).segment(1,3),r_m,r_s,umbra_res,'p','o');
		us = umbrapar(r.row(0).segment(7,3),r.row(0).segment(1,3),r_m,r_s,umbra_res,'u','o');
		ps_in = umbrapar(r.row(0).segment(7,3),r.row(0).segment(1,3),r_m,r_s,penumbra_res,'p','i');
		us_in = umbrapar(r.row(0).segment(7,3),r.row(0).segment(1,3),r_m,r_s,penumbra_res,'u','i');
		earth_isx = umbraisx(ps,us,r.row(0).segment(4,3),r_e);
		earth_isx_in = umbraisx(ps_in,us_in,r.row(0).segment(4,3),r_e);
		for (size_t i = 0; i < umbra_res; i++) {
			//Here we transpose the ps,us, ands earth_isx because printing works easier when printing rows instead of columns
			//of an eigen array
			param_ems_ofstream << setprecision(10) << r.row(0)(0) << "\t" << ps.transpose().row(i).format(TabFormat) <<"\t" << us.transpose().row(i).format(TabFormat) << "\n";
			umbraisx_ems_ofstream << setprecision(10) << r.row(0)(0) << "\t" << earth_isx.transpose().row(i).format(TabFormat) << "\n";
			latlong_ofstream << longlat(earth_isx.transpose().row(i),r.row(0).segment(4,3),r.row(0).segment(1,3), omega,int(r.row(0)(0))).transpose().format(CSVFormat) << "," << int(r.row(0)(0)) << "\n";
		}
		for (size_t i = 0; i < penumbra_res; i++) {
			param_ems_in_ofstream << setprecision(10) << r.row(0)(0) << "\t" << ps_in.transpose().row(i).format(TabFormat) <<"\t" << us_in.transpose().row(i).format(TabFormat) << "\n";
			umbraisx_ems_in_ofstream << setprecision(10) << r.row(0)(0) << "\t" << earth_isx_in.transpose().row(i).format(TabFormat) << "\n";
			latlong_in_ofstream << longlat(earth_isx_in.transpose().row(i),r.row(0).segment(4,3),r.row(0).segment(1,3), omega,int(r.row(0)(0))).transpose().format(CSVFormat) << "," << int(r.row(0)(0)) << "\n";
		}


		//Progress bar
		if (while_idx % (int(linecount)/20) == 0) {
			progress += 0.05;
			if (progress < 0.1) //condition never met when condition is defined as (progress == 0.05) !! Why ?
			{

				spacebar = 2; 
			}
			else if (progress > 0.05 && progress < 1)
			{
				spacebar = 1;
			}
			else {
				spacebar = 0;
			}

			cout << "Calculating shadow " << setw(20+spacebar) << "\033[32m" << int(progress*100.0) << "%" << setw(spacebar) << "\033[0m" << " [";
			pos = barWidth * progress;
			for (int j = 0; j < barWidth;  j++) {
				if (j < pos) cout << "=";
				else if (j == pos) cout << ">";
				else cout << " ";
		}
		
		cout << "]\r";
		cout.flush();
		}
		

	}
	cout << endl;
	pos_ems_ifstream.close();
	umbraisx_ems_ofstream.close();	
	param_ems_ofstream.close();
	latlong_ofstream.close();
	param_ems_in_ofstream.close();
	umbraisx_ems_in_ofstream.close();
	latlong_in_ofstream.close();
	


	

	







	
	auto duration = duration_cast<milliseconds>(stop-start);

	double r_end = sqrt(pow(r.row(1)(4),2)+pow(r.row(1)(5),2)+pow(r.row(1)(6),2));

	double r_m_end = sqrt(pow(r.row(1)(7)-r.row(1)(4),2)+pow(r.row(1)(8)-r.row(1)(5),2)+pow(r.row(1)(9)-r.row(1)(6),2));	

	cout << endl << "=> END OF INTEGRATION: " << duration.count() << " ms" << endl << endl;
	
	cout << "Distance Sun-Earth after " << days << " days: " <<  r_end/1e9 << " million km" << endl;
	
	cout << "Distance Earth-Moon after " << days << " days: " << r_m_end/1e6 << " thousand km\n" << endl;
	











	

	

	//Asking for eclipse calculation + eclipse calculation
	askTheUser:

	cout << "Do you want to look for solar eclipses in that timeframe ? ([Y],N):";
	string answer;
	getline(cin,input);
	if (input.empty()) {
		answer = "yes";
	}
	else {
		answer = input;
	}

	int eclipse_seconds, i_eclipse;
	//Actual execution depending on the answer:
	if (answer == string("no") || answer == string("n") || answer == string("No") || answer == string("N") || 		answer == string("NO")) {
	return 0;
	}
	else if (answer == string("yes") || answer == string("y") || answer == string("Yes") || answer == string("Y") ||	answer == string("YES") || answer == string("")) {
		ifstream pos_ems_ifstream("./data/pos_ems.dat");
        string line;
        //TODO Boolean variable type
		int flag = 0;
		//Looking for total eclpise condition
        while (getline(pos_ems_ifstream,line)) {
            stringstream ss(line);
            string value;
            int vl_idx = 0;
            while(getline(ss,value,'\t')) {
				//reusing r array
                r.row(0)(vl_idx) = stod(value);
                vl_idx += 1;                
            }
            Vector3d r_ms = r.row(0).segment(7,3)-r.row(0).segment(1,3);
			Vector3d r_es = r.row(0).segment(4,3)-r.row(0).segment(1,3);
			double d_ms = r_ms.norm();
			double d_es = r_es.norm();
			double phi_c = atan(r_e/d_es);
			double phi = acos((r_ms.dot(r_es)) / (r_ms.norm()*r_es.norm()) );
			if (d_ms < d_es && phi < phi_c) {
				eclipse_seconds = r.row(0)(0);
				flag = 1;
				break;
				//print_time("Jul 12 2024 00:00:00", eclipse_seconds);
			}
        }
		/*
		for (int s; s < N; s++) {
			Vector3d r_ms = r.row(s).segment(7,3)-r.row(s).segment(1,3);
			Vector3d r_es = r.row(s).segment(4,3)-r.row(s).segment(1,3);
			double d_ms = r_ms.norm();
			double d_es = r_es.norm();
			double phi_c = atan(r_e/d_es);
			//double phi_c = 0.1;
			double phi = acos((r_ms.dot(r_es)) / (r_ms.norm()*r_es.norm()) );
			//if ( s % int(N)/100 ) {cout << "phi: " << phi << "phi_c: " << phi_c << endl;}
			if (d_ms < d_es && phi < phi_c) {
				eclipse_seconds = s*t_delta;
				i_eclipse = s;
				flag = 1;
				break;
				//print_time("Jul 12 2024 00:00:00", eclipse_seconds);
			}
		}*/
		if (flag == 1)
		{
			cout << endl;
			cout << "NEXT ECLIPSE ON: "<< endl;
			//print_time("Jul 11 2024 12:00:00", eclipse_seconds);
			cout << "Condition met!" << endl;
		}
		else {
			cout << "\nNO ECLIPSE FOUND!\n(Eclipse condition not met in specified timeframe.)" << endl;
		}

		return 0;
	}
	else {
		goto askTheUser;
	}



	

}
