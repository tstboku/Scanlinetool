#include <fstream>
#include <iostream>
#include <vector>
#include <list>

#include <set>
#include <algorithm>
#include <stdexcept>

#include <stdint.h>
#include <stdio.h>
#include <math.h>

using namespace std;

double LX=3,LY=3,LZ=3;
double xmin=0, ymin=0, zmin=0;
uint64_t Nx= 50, Ny=50, Nz=50, N_planes;
uint64_t Nzy=Nz*Ny;
double DX=LX/Nx, DY=LY/Ny,DZ=LZ/Nz;
double vol_min=DX*DY*DZ;
typedef uint32_t  SET_NUMBER_TYPE;
//typedef uint64_t  SET_NUMBER_TYPE;

int body_next_num;
int N_CPU = 2;

class Cell{
public:
	//	double xc,yc,zc,vol;
	SET_NUMBER_TYPE  number;
	double zc(){
		return zmin+(number%Nz+0.5)*DZ;
	}
	double yc(){
		return ymin+((number/Nz)%Ny+0.5)*DY;
	}
	double xc(){
		return xmin+((number/Nzy)%Nx + 0.5)*DX;
	}
	double vol(){
		return vol_min;
	}
};

double zc_n(SET_NUMBER_TYPE number){
	return zmin+(number%Nz+0.5)*DZ;
}
double yc_n(SET_NUMBER_TYPE number){
	return ymin+((number/Nz)%Ny+0.5)*DY;
}
double xc_n(SET_NUMBER_TYPE number){
	return xmin+((number/Nzy)%Nx + 0.5)*DX;
}
double vol_n(){
	return vol_min;
}


typedef vector<SET_NUMBER_TYPE>::iterator CELL_ITERATOR;
typedef vector<SET_NUMBER_TYPE> CELL_CONTAINER;

//typedef set<SET_NUMBER_TYPE>::iterator CELL_ITERATOR;
//typedef set<SET_NUMBER_TYPE> CELL_CONTAINER;

//typedef list<SET_NUMBER_TYPE>::iterator CELL_ITERATOR;
//typedef list<SET_NUMBER_TYPE> CELL_CONTAINER;


class Body{
public:
	int body_id;
	int parent_id;
	CELL_CONTAINER body_cells;
	set<SET_NUMBER_TYPE> body_planes;
};

vector<Body> bodies;

class Disk
{
public:
	double a,b,c,d,R;
	double x,y,z;
	Disk(double kx, double ky, double kz,double xp,double yp,double zp,double r):a(kx),b(ky),c(kz),d(kx*xp+ky*yp+kz*zp),R(r),x(xp),y(yp),z(zp)
	{}

};
vector<Disk> planes;

bool larger(const Disk& lhs, const Disk& rhs)
{
	return lhs.R > rhs.R;
}




struct CELL_PARAMS{
     CELL_ITERATOR i_cell_0;
	 CELL_ITERATOR i_cell_1;

	 int N_active;
	 int N_cut; 
	 int flag_intersect; 
	 Body *bm;
	 Body *b2;
	 CELL_CONTAINER &b_parent;
	 Disk *D;
	 CELL_PARAMS(CELL_ITERATOR i_0, CELL_ITERATOR i_1, Body *_bm, Body *_b2, CELL_CONTAINER &b_p, Disk *_D):
		 i_cell_0(i_0), i_cell_1(i_1), bm(_bm), b2(_b2), b_parent(b_p), D(_D), flag_intersect(0),N_active(0),N_cut(0)
	 {}
};


void proc_split(void *params);
void proc_split2(Body *bm, Body *bcur, Disk *D);

int main ()
{
	uint64_t ix,iy,iz;
	double xc,yc;
	std::cout << "Plane intersections" << endl;


	try {
		uint64_t  Ncells = Nx*Ny*Nz;

		int num = 0;
		Body  b0;
		b0.body_id = 0;
		b0.parent_id = -1;  //First body
		 
		ifstream myfile ("geom.txt");
		myfile >> LX>>LY>>LZ;
		myfile >> xmin>>ymin>>zmin;
		myfile >>Nx>>Ny>>Nz;
		myfile >>N_planes;
		myfile.close();

		Nzy=Nz*Ny;
		DX=LX/Nx; DY=LY/Ny; DZ=LZ/Nz;
		vol_min=DX*DY*DZ;
			
		vector<double> xpl(N_planes,0);
		vector<double> ypl(N_planes,0);
		vector<double> zpl(N_planes,0);
		vector<double> kxpl(N_planes,0);
		vector<double> kypl(N_planes,0);
		vector<double> kzpl(N_planes,0);
		vector<double> rpl(N_planes,0);

		myfile.open ("c.txt");
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror1");
			myfile >> xpl[i];
		}
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror2");
			myfile >> ypl[i];
		}
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror3");
			myfile >> zpl[i];
		}
       myfile.close();

	   myfile.open("n.txt");
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror4");
			myfile >> kxpl[i];
		}
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror5");
			myfile >> kypl[i];
		}
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror6");
			myfile >> kzpl[i];
		}
		myfile.close();

		myfile.open ("r.txt");
		for (int i=0;i<N_planes;i++){
			if (!myfile.good()) throw std::invalid_argument ("Input data errror7");
			myfile >> rpl[i];
		}
		myfile.close();

		for (int i=0;i<N_planes;i++){
			Disk  D(kxpl[i],kypl[i],kzpl[i],xpl[i],ypl[i],zpl[i],rpl[i]);
			planes.push_back(D);
		}

		std::sort(planes.begin(),planes.end(), larger);

		bodies.reserve(10000);
		bodies.push_back(b0);
		DX=LX/Nx;
		DY=LY/Ny;
		DZ=LZ/Nz;
		printf("Disks:");
		for( ix=0; ix<Nx;ix++){
			xc = xmin+ (ix + 0.5)*DX;
			for( iy=0; iy<Ny;iy++){
				yc = ymin+ (iy+0.5)*DY;
				for( iz=0; iz<Nz;iz++){
					Cell c;
					c.number=num;
					bodies[0].body_cells.push_back(num);

					num++;
				}
			}
		}

		// Step 1
		int i_planes = 0;
		try {
			body_next_num = bodies.size();
			for(vector<Disk>::iterator iter_d = planes.begin(); iter_d != planes.end(); iter_d++)
			{
				Disk  D = *iter_d;
				printf(" %d",i_planes);
				std::vector<Body>::size_type size = bodies.size();
				for (std::vector<Body>::size_type i = 0; i < size; ++i)
				{
					int flag_intersect = 0;
					int N_active = 0;
					int N_cut = 0;
					CELL_CONTAINER &cur_cells = bodies[i].body_cells;

					Body bm;
					Body b2;

					CELL_PARAMS  data(cur_cells.begin(), cur_cells.end(), &bm, &b2, bodies[i].body_cells,  &D);
					proc_split(&data);
					flag_intersect = data.flag_intersect;
					N_active = data.N_active;
					N_cut = data.N_cut;

					if (flag_intersect && N_active==bodies[i].body_cells.size() && N_cut < bodies[i].body_cells.size()){


						//for(CELL_ITERATOR iter = bm.body_cells.begin(); iter != bm.body_cells.end();)
						//{
						//	bodies[i].body_cells.erase(*iter);

						//	iter++;
						//}

						//CELL_CONTAINER &  b_cur =bodies[i].body_cells;
						//for(CELL_ITERATOR iter = bm.body_cells.begin(); iter != bm.body_cells.end();)
						//{
						//	b_cur.erase(std::remove(b_cur.begin(), b_cur.end(), *iter), b_cur.end());
						//	iter++;
						//}


						bodies[i].body_cells = b2.body_cells;


						bodies[i].body_planes.insert(i_planes);
						bm.body_id = body_next_num;
						bm.parent_id = bodies[i].body_id;
						bm.body_planes.insert(bodies[i].body_planes.begin(), bodies[i].body_planes.end());
						bodies.push_back(bm);

						body_next_num++;
					}
				}

				i_planes++;
				printf("(%d)", body_next_num);
			}
		}
		catch (std::system_error& e) {
			std::cerr << "Error found" << std::endl;
		}
		// Step 2
		printf("\n Checks Rcell < R disks for all bodies \n");
		for (unsigned i = bodies.size(); i-- > 1; ){
			// 0  initial body don't checked
			CELL_CONTAINER &cur_cells = bodies[i].body_cells;
			set<SET_NUMBER_TYPE> &cur_planes = bodies[i].body_planes;
			printf(" Body %d: ", i);
			double vol = 0;
			int planes_ok=-1;

			
			for(set<SET_NUMBER_TYPE>::iterator iter_d = cur_planes.begin(); iter_d != cur_planes.end(); iter_d++){
				Disk  *D = &planes[*iter_d];
				int N_active=0;
				for(CELL_ITERATOR cl = cur_cells.begin(); cl != cur_cells.end(); cl++){

					int i_cell = *cl;
					//double t_xc=cells[i_cell].xc();
					double t_xc=xc_n(i_cell);

					//double t_yc=cells[i_cell].yc();
					double t_yc=yc_n(i_cell);

					//double t_zc=cells[i_cell].zc();
					double t_zc=zc_n(i_cell);

					double d1 = t_xc - D->x;
					double d2 =t_yc -D->y;
					double d3 = t_zc- D->z;

					double Lpc2 = (d1*d1 + d2*d2 +d3*d3);
					double Lp2   = (d1* D->a + d2*D->b +d3* D->c);
					Lp2 = Lp2*Lp2/(D->a*D->a + D->b*D->b +D->c*D->c);

					double D_axis;
					if (Lp2 > Lpc2  ){
						D_axis = 0;
						if ((Lp2 - Lpc2)> 1e-10) throw std::invalid_argument ("SQRT negative");
					}
					else 
						D_axis = sqrt(Lpc2 - Lp2);

					if (D_axis <= D->R){   //
						N_active ++;
					}
				}
				if (N_active < cur_cells.size()){
					// merge to parent body
					planes_ok = *iter_d;
					break;
				}
			}
			if (planes_ok!=-1){
				// merge to parent body
				Body &b_parent = bodies[bodies[i].parent_id];
				b_parent.body_planes.erase(planes_ok);
				//b_parent.body_cells.insert(cur_cells.begin(),cur_cells.end());
				b_parent.body_cells.insert(b_parent.body_cells.end(), cur_cells.begin(),cur_cells.end());

				bodies.erase(bodies.begin()+i);
				printf(" deleted \n", i);
			}
			else{
				printf(" OK \n", i);
			}
		}
		printf(" Body 0 OK\n");

		double tot = 0;
		int N_tot=0, N_bd=0;
		FILE *f_out; 
		f_out = fopen("result.txt","w");

		printf("\n");
		for(vector<Body>::iterator bd = bodies.begin(); bd != bodies.end(); bd++){
			CELL_CONTAINER &cur_cells = bd->body_cells;
			printf(" Body %d: Cells: ", bd->body_id);
			double vol = 0;
			int N_c=0;
			for(CELL_ITERATOR cl = cur_cells.begin(); cl != cur_cells.end(); cl++){
				//			printf(" %d ", *cl);
				// fprintf (f_out, " %d", *cl);
//				vol +=  cells[*cl].vol(); 
				vol +=  vol_n(); 

				N_c++;
			}
			printf(" %d Volume=%lf \n",N_c, vol);
			fprintf(f_out, "%lf\n",vol);
			tot += vol;
			N_tot += N_c;
			N_bd++;
		}
		printf(" Total Bodies %d Cells %d Volume %lf ",N_bd, N_tot, tot);
		//fprintf (f_out, " Total Cells %d Volume %lf ", N_tot, tot);
		fclose(f_out);

	}
	catch (std::bad_alloc& ba){
		cerr << "bad_alloc caught: " << ba.what() << endl;
		cerr << "Decrease Nx Ny Nz "<< endl;
	}
	catch (std::exception& e) {
		cerr << "Error: " << e.what() << endl;
	}
	catch (...) {
		cerr << "Error2" << endl;
	}

}


void proc_split(void *params)
{ 
	CELL_PARAMS *pa = (CELL_PARAMS*)params;

	for(CELL_ITERATOR cl = pa->i_cell_0; cl != pa->i_cell_1; ++cl)
	{
		int i_cell = *cl;
		double t_xc=xc_n(i_cell);
		double t_yc=yc_n(i_cell);
		double t_zc=zc_n(i_cell);


		double d1 = t_xc - pa->D->x;
		double d2 =t_yc -pa->D->y;
		double d3 = t_zc- pa->D->z;
		double D_plane = t_xc * pa->D->a + t_yc * pa->D->b + t_zc * pa->D->c - pa->D->d;

		pa->N_active ++;
		if( D_plane < 0 ){
			// plane split the body
			pa->N_cut++;
			if(!pa->flag_intersect) {
				pa->flag_intersect = 1;
			}
			//pa->bm->body_cells.insert(i_cell);
			pa->bm->body_cells.push_back(i_cell);
		} else {
			pa->b2->body_cells.push_back(i_cell);
		}
	}
}
