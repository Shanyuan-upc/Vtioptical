/* Copyright (c) China university of Petroelum, 2016.*/
/* All rights reserved.                       */

/* FDMOD2_PML_staggered: $Revision: 1.1 $ ; $Date: 2016/03/31 08:58:00 $	*/


#include "main.h"
#include"mpi.h"



/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" FDMOD2_PML_STAGGERED - Finite-Difference MODeling (high order) ",
" for acoustic wave equation with PML absorbing boundary conditions.	",
" 									",
" 									",
" Required Parameters:							",
" <vfile		file containing velocity[nx][nz]		",
" >wfile		file containing waves[nx][nz] for time steps	",
" nx=			number of x samples (2nd dimension)		",
" nz=			number of z samples (1st dimension)		",
" xs=			x coordinates of source				",
" zs=			z coordinates of source				",
" tmax=			maximum time					",
" 									",
" Optional Parameters:							",
" nt=1+tmax/dt		number of time samples (dt determined for stability)",
" mt=1			number of time steps (dt) per output time step	",
" 									",
" dx=1.0		x sampling interval				",
" fx=0.0		first x sample					",
" dz=1.0		z sampling interval				",
" fz=0.0		first z sample					",
" 									",
" fmax = vmin/(10.0*h)	maximum frequency in source wavelet		",
" fpeak=0.5*fmax	peak frequency in ricker wavelet		",
" 									",
" name_density=		input file containing density[nx][nz]		",
" vsx=			x coordinate of vertical line of seismograms	",
" hsz=			z coordinate of horizontal line of seismograms	",
" vsfile=		output file for vertical line of seismograms[nz][nt]",
" name_shot=		output file for horizontal line of seismograms[nx][nt]",
" ssfile=		output file for source point seismograms[nt]	",
" verbose=0		=1 for diagnostic messages, =2 for more		",
" 									",
" abs=1,1,1,1		Absorbing boundary conditions on top,left,bottom,right",
" 			sides of the model. 				",
" 		=0,1,1,1 for free surface condition on the top		",
"                                                                       ",
" ...PML parameters....                                                 ",
" pml_max=1000.0        PML absorption parameter                        ",
" pml_thick=0           half-thickness of pml layer (0 = do not use PML)",
" 									",
" Notes:								",
" This program uses the traditional explicit second order differencing	",
" method. 								",
" 									",
" Two different absorbing boundary condition schemes are available. The ",
" first is a traditional absorbing boundary condition scheme created by ",
" Hale, 1990. The second is based on the perfectly matched layer (PML)	",
" method of Berenger, 1995.						",
" 									",
NULL};

/*
 * Authors:  CPU:Jidong  Yang, SIIG SWPI Research Group
 *
 */
/**************** end self doc ********************************/


/* PML related global variables */
float pml_max=0;
int pml_thick=0;
int pml_thickness=0;

/* CPML related global variables */
float *d_x,*K_x,*alpha_x,*a_x,*b_x,*d_x_half,*K_x_half,*alpha_x_half,*a_x_half,*b_x_half; 
float *d_z,*K_z,*alpha_z,*a_z,*b_z,*d_z_half,*K_z_half,*alpha_z_half,*a_z_half,*b_z_half; 

float sigma, sigma_ex, sigma_ez, sigma_mx, sigma_mz;


/* functions used internally*/
void read_par(FILE*fp_par,char *name_velo,char *name_epsilon,char *name_delta,char *name_density,char *name_ele,
	char *name_shot,char *name_snap,int *flag_snap,char *name_snapu,int *flag_snapu,char *name_snapv,int *flag_snapv,int *nx,int *nz,int*nxc,
	int *nzc,float *dx,float *fx,float *dz,float *fz,int *ns,int *is_beg,
	int *is_end,int *xsn,int *zsn,int *ds,float *tmax,float *dt,
	float *dtout,int*nmt,int *ndtr,float *fpeak,
	int *order,int *ssize,float *tol,float *lambda,float *gama,int *l,int *iteration,float *pml_max,int *pml_thick,int *verbose,int *itype_acq);
void load_coeff(int order,float **coeff);
static void pml_init (int nx, int nz, float dx, float dz, float dt,
	float **dvv, float **od, int verbose,int order);
void load_dvvc_odc(int is,int ds,float **dvv,float **epsilon,float **delta,float **od,float **dvvc,float **epsilonc,float **deltac,
	float **odc,int nx,int nz,int nxc,int nzc,int *hele,float *ele,
	int pml_thick,int pml_thickness,int itype_acq);
void ptsrc (int ixs, int izs,int nx, float dx, float fx,int nz, float dz,
	 float fz,float dt, float t, float fmax, float fpeak, float *tdelay,
 	 float **s,int pml_thick,int pml_thickness);
void tstep2 (int nx, float dx, int nz, float dz, float dt,float **dvv, float **epsilon,float **delta,float **od, float **s,float **p,float **Ep,float**ud,float **E,float **u,float **v,float **pp,float **pm,float **dpdx,
float **dpdz,float **vx,float **vz,float **dpdxx,float **dpdzz,int order,int ssize,float tol,float lambda,float gama,int iteration,float**coeff,int pml_thick,int pml_thickness);

static void write_su(int is,int ds,int xsn,int nx,int ndtr,float fx, 
	int nt,float dt,float dtout,float dx,float **hs,FILE*fp_out);



int main(int argc,char**argv)
{
	register int ix,iz,it,is;	/* counters 			*/
	int nx,nz,nt,nmt;		/* x,z,t,tsizes 		*/
	int nxc,nzc;			/* calculation grid size	*/

	/* finite difference parameter*/
	int order;		/* finite difference order		*/
	int ssize;
	float tol;
	float lambda;		/* finite difference order		*/
	float gama;
	int l;		/* finite difference order		*/
	int iteration;		/* finite difference order		*/
	float **coeff;		/* diference coefficient 		*/

	/* shot paremter */	
	int ns;			/* shot number				*/
	int is_beg,is_end;	/* shot number of beginning and ending 	*/
	int xsn,zsn;		/* first shot position			*/		
	int xsn1,zsn1;		/* first shot position			*/		
	int ds;			/* grid spacing of shot			*/
	int ndtr;		/* trace spacing for output trace	*/

	/* grid parameter */
	float dx;		/* horizontal grid spacing		*/
	float dz;		/* vertical grid spacing		*/
	float h;		/* minumum spatial sample interval	*/
	float fx;		/* first x value of model		*/
	float fz;		/* first z value of model		*/

	/* wavelet parameter*/
	float dt,dtout;		/* time sample interval */
	float fmax;		/* maximum temporal frequency allowable */
	float fpeak;		/* peak frequency of ricker wavelet 	*/
	float tdelay=0.;	/* time delay of source beginning 	*/
	float tmax;		/* maximum time to compute */
	float t;		/* time parameter */
	
	/* maximum and minimum velocity and density */
	float vmin;		/* minimum wavespeed in vfile 		*/
	float vmax;		/* maximum wavespeed in vfile 		*/
	float dmin=0.0;		/* minimum density in name_density 	*/
	float dmax=0.0;		/* maximum density in name_density 	*/

	/* parameter and wavefield array*/
	float **s;		/* array of source pressure values */
	float **dvv;		/* array of velocity values from vfile */
	float **epsilon;		/* array of velocity values from vfile */
	float **delta;		/* array of velocity values from vfile */
	float **od;		/* array of density values from name_density */
	float **dvvc;		/* array of calculatiing velocity values  */
	float **epsilonc;		/* array of calculatiing velocity values  */
	float **deltac;		/* array of calculatiing velocity values  */
	float **odc;		/* array of calculatiing density values  */

	/* pressure field arrays */
	float **p,**pp,**pm,**dpdx,**dpdz;	/* pressure field at t-0.5dt */
	float **u,**v,**E,**Ep,**ud;	/* pressure field at t-0.5dt */

        float **vx,**dpdxx;             /* velocity x field at t */
        float **vz,**dpdzz;             /* velocity z field at t */

	float *ws;		/* record near the source*/

	/* elevation parameter*/
	float *ele;	/* elevation array			*/
	int *hele;	/* vertical sample of horiz rec line 	*/
	
	/* output data arrays */
	float **hs;		/* seismograms from horiz receiver line */

	/* flag of wavefield output */
	int flag_snap;
	int flag_snapu;
	int flag_snapv;
	
	/* file names */
	char name_density[256]="";	/* density file name 				*/
	char name_shot[256]="";		/* horiz receiver seismogram line file name 	*/
	char name_snap[256]="";		/* wavefield snap slice at different time	*/
	char name_snapu[256]="";		/* wavefield snap slice at different time	*/
	char name_snapv[256]="";		/* wavefield snap slice at different time	*/
	char name_velo[256]="";		/* p-wave velocity file 			*/
	char name_epsilon[256]="";		/* p-wave velocity file 			*/
	char name_delta[256]="";		/* p-wave velocity file 			*/
	char name_ele[256]="";		/* elevation file for irregular topography	*/
	
	/* input file pointers */
	FILE *fp_velo; 			/* pointer to input velocity data 	*/
	FILE *fp_epsilon; 			/* pointer to input velocity data 	*/
	FILE *fp_delta; 			/* pointer to input velocity data 	*/
	FILE *fp_density;		/* pointer to input density data file 	*/
	FILE *fp_ele;			/* pointer to input elavation file	*/
	FILE *fp_shot;			/* pointer to output shot record file	*/
	FILE *fp_snap;			/* pointer to output wavefield 		*/
	FILE *fp_snapu;			/* pointer to output wavefield 		*/
	FILE *fp_snapv;			/* pointer to output wavefield 		*/
	FILE *fp_snapud;			/* pointer to output wavefield 		*/
	
	/* SEGY fields */
	long tracl=0;		/* trace number within a line */
	long tracr=0;		/* trace number within a reel */


	/* MPI parameter */
	int myid;	/* current id of processor 	*/
	int np;		/* number of processors		*/

	int verbose;
	int itype_acq;

	/* parameter file*/
	char name_par[256],*name_siig;
	FILE* fp_par,*fp_siig;
	char temp_char[256];
	/* time tick parameter*/
	double time_beg,time_end;

	
	/* open parameter file*/
	name_siig = "par_modelling.siig";
	if((fp_siig=fopen(name_siig,"rb"))==NULL)
		printf("Open par file error!\n"),exit(0);
	fgets(temp_char,240,fp_siig);
	printf("%s",temp_char);
	fscanf(fp_siig,"%s\n",name_par);
	printf("%s\n",name_par);
	fclose(fp_siig);

	/* read parameter from par file*/
	if((fp_par=fopen(name_par,"rb"))==NULL)
		printf("Open par file error!\n"),exit(0);
	read_par(fp_par,name_velo,name_epsilon,name_delta,name_density,name_ele,name_shot,name_snap,
		&flag_snap,name_snapu,&flag_snapu,name_snapv,&flag_snapv,&nx,&nz,&nxc,&nzc,&dx,&fx,&dz,&fz,&ns,&is_beg,&is_end,
		&xsn,&zsn,&ds,&tmax,&dt,&dtout,&nmt,&ndtr,&fpeak,&order,&ssize,&tol,&lambda,&gama,&l,&iteration,&pml_max,
		&pml_thick,&verbose,&itype_acq);
	fclose(fp_par);

	/* set order=4*/
	order = 4;
	
	/* open MPI */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	time_beg = MPI_Wtime();

	/* get the beggining time*/
	if(myid==0)
	{
		time_beg = MPI_Wtime();
	}

	//MPI_Barrier(MPI_COMM_WORLD);

	/* judge shot number  */
	if (ns<=0) printf(" shot number must be larger than 1!\n"),exit(0);
                printf("The source xsn %d zsn %d nxc %d nzc %d \n",xsn,zsn,nxc,nzc);
        if (xsn<1 || xsn>nxc || zsn <1 || zsn >nzc) 
                printf("The source must be in the model!\n"),exit(0);

	/* boundary thichness */
	pml_thickness = 2 * pml_thick;

	/* determine number of time steps required to reach maximum time */
	nt = 1+tmax/dt;

	/* open velocity file and read velocity */
	if((fp_velo = fopen(name_velo,"rb"))==NULL)
	{
		printf("cann't open velo file\n"),exit(0);
	}
	else
	{
		/* read velocities */
		dvv = alloc2float(nz,nx);
		if(fread(dvv[0],sizeof(float),nx*nz,fp_velo)!=nx*nz)
			printf("read velocity error!\n"),exit(0);
	}

	/* open epsilon file and read epsilon */
	if((fp_epsilon = fopen(name_epsilon,"rb"))==NULL)
	{
		printf("cann't open epsilon file\n"),exit(0);
	}
	else
	{
		/* read epsilon */
		epsilon = alloc2float(nz,nx);
		if(fread(epsilon[0],sizeof(float),nx*nz,fp_epsilon)!=nx*nz)
			printf("read epsilon error!\n"),exit(0);
	}
		printf("epsilon=%f file\n",epsilon[100][110]);

	/* open delta file and read delta */
	if((fp_delta = fopen(name_delta,"rb"))==NULL)
	{
		printf("cann't open delta file\n"),exit(0);
	}
	else
	{
		/* read delta */
		delta = alloc2float(nz,nx);
		if(fread(delta[0],sizeof(float),nx*nz,fp_delta)!=nx*nz)
			printf("read delta error!\n"),exit(0);
	}

	/* if specified, read densities */
	od = alloc2float(nz,nx);
	if (strcmp(name_density,"NULL")==0) 
	{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				od[ix][iz] = 1.0;
		dmin = dmax = 1.0;
	}
	else 
	{
		if((fp_density=fopen(name_density,"r"))==NULL)
			printf("cannot open name_density=%s\n",name_density),exit(0);
		if (fread(od[0],sizeof(float),nx*nz,fp_density)!=nx*nz)
			printf("error reading name_density=%s\n",name_density),exit(0);
		fclose(fp_density);
		dmin = dmax = od[0][0];
		for (ix=0; ix<nx; ++ix) {
			for (iz=0; iz<nz; ++iz) {
				dmin = MIN(dmin,od[ix][iz]);
				dmax = MAX(dmax,od[ix][iz]);
			}
		}
	}

	/* read or define receiver line position */
	ele = alloc1float(nx);
	hele = alloc1int(nxc);
	if(strcmp(name_ele,"NULL")==0)
	{
		for (ix=0; ix<nx; ++ix)
			ele[ix] = 0.0;
	}
	else
	{
		if((fp_ele=fopen(name_ele,"r"))==NULL)
			printf("cannot open name_density=%s\n",name_density),exit(0);
		if (fread(ele,sizeof(float),nx,fp_ele)!=nx)
			printf("error reading name_ele=%s\n",name_ele),exit(0);
		fclose(fp_ele);
	}

	/* if requested, open file and allocate space for seismograms */
	/* ... horizontal line of seismograms */
	if (*name_shot!='\0') {
		if((fp_shot=fopen(name_shot,"wb"))==NULL)
			printf("cannot open name_shot=%s\n",name_shot),exit(0);
		hs = alloc2float(nt,nxc);
	} else {
		hs = NULL;
	}

	/* open wavesnap file */
	if(flag_snap)
	{
		if((fp_snap=fopen(name_snap,"wb"))==NULL)	
			printf("cannot open wavesnapfile!\n"),exit(0);	
	}
	else
	{
		fp_snap=NULL;
	}	

	if(flag_snapu)
	{
		if((fp_snapu=fopen(name_snapu,"wb"))==NULL)	
			printf("cannot open wavesnapufile!\n"),exit(0);	
	}
	else
	{
		fp_snapu=NULL;
	}	
	if(flag_snapv)
	{
		if((fp_snapv=fopen(name_snapv,"wb"))==NULL)	
			printf("cannot open wavesnapvfile!\n"),exit(0);	
	}
	else
	{
		fp_snapv=NULL;
	}	
		if((fp_snapud=fopen("snapshotud.dat","wb"))==NULL)	
			printf("cannot open wavesnapfile!\n"),exit(0);	
	/* allocate workplace*/
	coeff = alloc2float(order,order);
	dvvc = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	epsilonc = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	deltac = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	odc = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	s = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	p = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	pm = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	pp = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        dpdx = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        dpdz = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	E = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	Ep = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	ud = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	u = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
	v = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        vx = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        vz = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        dpdxx = alloc2float(nzc+pml_thickness,nxc+pml_thickness);
        dpdzz = alloc2float(nzc+pml_thickness,nxc+pml_thickness);

	ws = alloc1float(nt);

	/* obtain the difference coefficient */
	load_coeff(order,coeff);

	/* determine minimum and maximum velocities */
	vmin = vmax = dvv[0][0];
	for (ix=0; ix<nx; ++ix) {
		for (iz=0; iz<nz; ++iz) {
			vmin = MIN(vmin,dvv[ix][iz]);
			vmax = MAX(vmax,dvv[ix][iz]);
		}
	}

        /* determine mininum spatial sampling interval */
	h = MIN(fabs(dx),fabs(dz));

	/* determine time sampling interval to ensure stability */
	if(vmax*dt/h>0.5)
	{
		printf("Input dt is too larger or space interval is too small!"),exit(0);
		//dt = h/(2.0*vmax);		
	}

	/* determine maximum temporal frequency to avoid dispersion */
	fmax = vmin/(5.0*h);
	
	/* compute or set peak frequency for ricker wavelet */
//	if(fpeak>0.5*fmax){
//		printf("fpeak is too large and is reset as %f Hz\n",0.5*fmax);	
//		fpeak = 0.5*fmax;
//	}	

	/* compute density*velocity^2 and 1/density and zero time slices */	
	for (ix=0; ix<nx; ++ix) {
		for (iz=0; iz<nz; ++iz) {
			dvv[ix][iz] = dvv[ix][iz]*dvv[ix][iz];
			od[ix][iz] = 1.0/od[ix][iz];
		}
	}

	/* compute elevation/dz*/
	for (ix=0; ix<nx; ++ix) 
	{
		ele[ix] /= dz;
	}	
	
	/* if verbose, print parameters */
	if (verbose) {
		printf("nx = %d\n",nx);
		printf("dx = %g\n",dx);
		printf("nz = %d\n",nz);
		printf("dz = %g\n",dz);
		printf("nt = %d\n",nt);
		printf("dt = %g\n",dt);
		printf("tmax = %g\n",tmax);
		printf("fmax = %g\n",fmax);
		printf("fpeak = %g\n",fpeak);
		printf("vmin = %g\n",vmin);
		printf("vmax = %g\n",vmax);
		printf("nmt = %d\n",nmt);
                printf("pml_max = %g\n",pml_max);
                printf("pml_half = %d\n",pml_thick);
                printf("pml_thickness = %d\n",pml_thickness);
		if (dmin==dmax) {
			printf("constant density\n");
		} else {
			printf("name_density=%s\n",name_density);
			printf("dmin = %g\n",dmin);
			printf("dmax = %g\n",dmax);
		}
	}

	/* initialize MPL boundary coefficients*/
	if (pml_thick != 0)
		pml_init (nxc, nzc, dx, dz, dt, dvv, od, verbose,order);


		zsn1=zsn;
		xsn1=xsn;	
	/* shot loop for seismic modelling*/
	for(is=is_beg+myid;is<=is_end;is += np)
	{
	        /* set shot vertical position */
		if(itype_acq==0){
		  xsn = (is-1)*ds+xsn1; 
	          zsn = ceil(ele[(is-1)*ds+xsn]+zsn1);
        	}else{
	          zsn = ceil(ele[(is-1)*ds+xsn]+zsn1);
		}	
	
		/* zeros p,vx,vz */
		memset(p[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(pp[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(pm[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(Ep[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(ud[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(E[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(u[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
		memset(v[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(vx[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(vz[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(dpdxx[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(dpdzz[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(dpdx[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));
                memset(dpdz[0],0,sizeof(float)*(nxc+pml_thickness)*(nzc+pml_thickness));

		/* load velocity and density into calculating grid */
		load_dvvc_odc(is,ds,dvv,epsilon,delta,od,dvvc,epsilonc,deltac,odc,nx,nz,nxc,nzc,hele,ele,
			pml_thick,pml_thickness,itype_acq);



		/* loop  over time steps */
		for (it=0,t=0.0; it<nt; ++it,t+=dt) 
		{
			/* update source function */
			ptsrc(xsn,zsn,nxc,dx,fx,nzc,dz,fz,dt,t,
			      fmax,fpeak,&tdelay,s,pml_thick,pml_thickness);


#if 0
/* for test */
     for (ix=0;ix<nxc+2*pml_thick;ix++){
       printf("Hello: d=%f K=%f alpha=%f a=%f b=%f %f %f %f %f %f\n",d_x[ix],K_x[ix],alpha_x[ix],a_x[ix],b_x[ix],d_x_half[ix],
            K_x_half[ix],alpha_x_half[ix],a_x_half[ix],b_x_half[ix]);
     }



for (ix=0;ix<nzc+2*pml_thick;ix++){
       printf("Hello: iz=%d d=%f K=%f alpha=%f a=%f b=%f %f %f %f %f %f\n",ix,d_z[ix],K_z[ix],alpha_z[ix],a_z[ix],b_z[ix],d_z_half[ix],
            K_z_half[ix],alpha_z_half[ix],a_z_half[ix],b_z_half[ix]);
     }

exit(0);
#endif

			/* do one time step */
			tstep2(nxc,dx,nzc,dz,dt,dvvc,epsilonc,deltac,odc,s,p,Ep,ud,E,u,v,pp,pm,dpdx,dpdz,vx,vz,dpdxx,dpdzz,order,ssize,tol,lambda,gama,iteration,coeff,pml_thick,pml_thickness);

#if 0			
			/* save source wavelet*/
			ws[it] = p[xsn+pml_thick][zsn+pml_thick];

			/* set the matched laryer boundary as 0*/
			for(ix=0;ix<nxc+pml_thickness;ix++)
			{
				vx[ix][0] = 0.;
				vz[ix][0] = 0.;
				vx[ix][nzc+pml_thickness-1] = 0.;
				vz[ix][nzc+pml_thickness-1] = 0.;
			}
			for(iz=0;iz<nzc+pml_thickness;iz++)
			{
				vx[0][iz] = 0.;
				vz[0][iz] = 0.;
				vx[nxc+pml_thickness-1][iz] = 0.;
				vz[nxc+pml_thickness-1][iz] = 0.;
			}
#endif
			/* output wavefiled snaps */
			if (it%nmt==0) 
			{
				printf("is = %d | it/nt = %d / %d\n  xsn= %d",is,it,nt,xsn);
				
				if(flag_snap )//&& it==1200
				{
					/* set selected trace header fields trace by trace */
					for (ix=0 ; ix < nxc+pml_thick+pml_thick ; ++ix) 
					{
						/* output traces of data cube */
						fwrite(&p[0+ix][0],
							sizeof(float),nzc+pml_thick+pml_thick,fp_snap);
					}
				}	
				if(flag_snapu )//&& it==1200
				{
					/* set selected trace header fields trace by trace */
					for (ix=0 ; ix < nxc+pml_thick+pml_thick ; ++ix) 
					{
						/* output traces of data cube */
						fwrite(&u[0+ix][0],
							sizeof(float),nzc+pml_thick+pml_thick,fp_snapu);
					}
				}	
				if(flag_snapv )//&& it==1200
				{
					/* set selected trace header fields trace by trace */
					for (ix=0 ; ix < nxc+pml_thick+pml_thick ; ++ix) 
					{
						/* output traces of data cube */
						fwrite(&E[0+ix][0],
							sizeof(float),nzc+pml_thick+pml_thick,fp_snapv);
					}
				}	
					for (ix=0 ; ix < nxc+pml_thick+pml_thick ; ++ix) 
					{
						/* output traces of data cube */
						fwrite(&ud[0+ix][0],
							sizeof(float),nzc+pml_thick+pml_thick,fp_snapud);
					}
				//printf("is= %d  p = %.10f\n",is,p[00][300]);
			}
			
			
			/* if requested, save horizontal line of seismograms */
			if (hs!=NULL) {
			    #pragma omp parallel for shared(nxc,hs,p,pml_thick,hele,it) private(ix) 
				for (ix=0; ix<nxc; ++ix)
					hs[ix][it] = p[pml_thick+ix][pml_thick+hele[ix]];
			}
		}/* time loop over*/


		/* if requested, write horizontal line of seismograms */
		write_su(is,ds,xsn,nxc,ndtr,fx,nt,dt,dtout,dx,hs,fp_shot);
		printf("is  %d shot is writen \n",is);
//	MPI_Barrier(MPI_COMM_WORLD);
	#if 0
		{
			FILE* fp_ws;
			fp_ws = fopen("wavlet.dat","wb");
			fwrite(ws,sizeof(float),nt,fp_ws);
			fclose(fp_ws);
		}
      #endif
	}/* shot loop over*/




	/* close file*/
	fclose(fp_shot);
	if(fp_snap!=NULL) fclose(fp_snap);
	if(fp_snapu!=NULL) fclose(fp_snapu);
	if(fp_snapv!=NULL) fclose(fp_snapv);
	if(fp_snapud!=NULL) fclose(fp_snapud);

	/* free space before returning */
	free2float(s);
	free2float(dvv);
	free2float(epsilon);
	free2float(delta);
	free2float(p);
        free2float(dpdx);
        free2float(dpdz);
        free2float(vx);
        free2float(vz);
        free2float(dpdxx);
        free2float(dpdzz);
	free2float(pm);
	free2float(pp);
	free2float(Ep);
	free2float(ud);
	free2float(E);
	free2float(u);
	free2float(v);
	free1float(ws);
	if (od!=NULL) free2float(od);
	if (odc!=NULL) free2float(odc);
	if (dvvc!=NULL) free2float(dvvc);
	if (epsilonc!=NULL) free2float(epsilonc);
	if (deltac!=NULL) free2float(deltac);
	if (hs!=NULL) free2float(hs);
	if(coeff!=NULL) free2float(coeff);
	if(ele!=NULL) free1float(ele);
	if(hele!=NULL) free1int(hele);
	
	time_end = MPI_Wtime();

	if(myid==0)
		printf("Time | %f s\n",(double)(time_end-time_beg));
	/* close MPI*/
	MPI_Finalize();

	return(SIIG_Exit());
}


/* read parameter from *.par file*/
void read_par(FILE*fp_par,char *name_velo, char *name_epsilon,char *name_delta,char *name_density,char *name_ele,
	char *name_shot,char *name_snap,int *flag_snap,char *name_snapu,int *flag_snapu,char *name_snapv,int *flag_snapv,int *nx,int *nz,int*nxc,
	int *nzc,float *dx,float *fx,float *dz,float *fz,int *ns,int *is_beg,
	int *is_end,int *xsn,int *zsn,int *ds,float *tmax,float *dt,float *dtout,
	int*nmt,int *ndtr,float *fpeak,int *order,int *ssize,float *tol,float *lambda,float *gama,int *l,int *iteration,float *pml_max,int *pml_thick,
	int *verbose,int *itype_acq)
/*******************************************************************************
*******************************************************************************/
{
	char temp_char[250]={0};	
	
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_velo);
	printf("%s\n",name_velo);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_epsilon);
	printf("%s\n",name_epsilon);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_delta);
	printf("%s\n",name_delta);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_density);
	printf("%s\n",name_density);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_ele);
	printf("%s\n",name_ele);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_shot);
	printf("%s\n",name_shot);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_snap);
	printf("%s\n",name_snap);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",flag_snap);
	printf("%d\n",*flag_snap);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_snapu);
	printf("%s\n",name_snapu);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",flag_snapu);
	printf("%d\n",*flag_snapu);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%s\n",name_snapv);
	printf("%s\n",name_snapv);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",flag_snapv);
	printf("%d\n",*flag_snapv);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",nx);
	printf("%d\n",*nx);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",nz);
	printf("%d\n",*nz);
	
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",nxc);
	printf("%d\n",*nxc);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",nzc);
	printf("%d\n",*nzc);
	
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",dx);
	printf("%f\n",*dx);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",dz);
	printf("%f\n",*dz);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",fx);
	printf("%f\n",*fx);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",fz);
	printf("%f\n",*fz);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",ns);
	printf("%d\n",*ns);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",is_beg);
	printf("%d\n",*is_beg);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",is_end);
	printf("%d\n",*is_end);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",xsn);
	printf("%d\n",*xsn);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",zsn);
	printf("%d\n",*zsn);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",ds);
	printf("%d\n",*ds);
	
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",tmax);
	printf("%f\n",*tmax);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",dt);
	printf("%f\n",*dt);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",dtout);
	printf("%f\n",*dtout);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",nmt);
	printf("%d\n",*nmt);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",ndtr);
	printf("%d\n",*ndtr);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",fpeak);
	printf("%f\n",*fpeak);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",order);
	printf("%d\n",*order);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",ssize);
	printf("%d\n",*ssize);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",tol);
	printf("%f\n",*tol);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",lambda);
	printf("%f\n",*lambda);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",gama);
	printf("%f\n",*gama);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",iteration);
	printf("%d\n",*iteration);
		
	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%f\n",pml_max);
	printf("%f\n",*pml_max);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",pml_thick);
	printf("%d\n",*pml_thick);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",verbose);
	printf("%d\n",*verbose);

	fgets(temp_char,240,fp_par);
	printf("%s",temp_char);
	fscanf(fp_par,"%d\n",itype_acq);
	printf("%d\n",*itype_acq);

}

/* determine 1-order finite difference coefficients of staggered grids*/
void load_coeff(int order,float **coeff)
{
	if(order==1)
	{
		coeff[0][0] = 1.0;
	}
	else if(order==2)
	{
		coeff[0][0] = 1.0;
		coeff[1][0] = 1.12500,coeff[1][1] = -0.04166667;
	}
	else if(order==3)
	{
		coeff[0][0] = 1.0;
		coeff[1][0] = 1.12500,coeff[1][1] = -0.04166667;
		coeff[2][0] = 1.1718750,coeff[2][1] = -0.065104167,
		coeff[2][2] = 0.0046875;
	}
	else if(order==4)
	{
		coeff[0][0] = 1.0;
		coeff[1][0] = 1.12500,coeff[1][1] = -0.04166667;
		coeff[2][0] = 1.1718750,coeff[2][1] = -0.065104167,
		coeff[2][2] = 0.0046875;
		coeff[3][0] = 1.196289,coeff[3][1] = -0.0797526,
		coeff[3][2] = 0.009570313,coeff[3][3] = -0.0006975447;
	}
	else if(order==5)
	{
		coeff[0][0] = 1.0;
		coeff[1][0] = 1.12500,coeff[1][1] = -0.04166667;
		coeff[2][0] = 1.1718750,coeff[2][1] = -0.065104167,
		coeff[2][2] = 0.0046875;
		coeff[3][0] = 1.196289,coeff[3][1] = -0.0797526,
		coeff[3][2] = 0.009570313,coeff[3][3] = -0.0006975447;
		coeff[4][0] = 1.211243,coeff[4][1] = -0.08972168 ,
		coeff[4][2] = 0.01384277,coeff[4][3] = -0.00176566,
		coeff[4][4] = 0.0001186795;
	}
	else if(order==6)
	{
		coeff[0][0] = 1.0;
		coeff[1][0] = 1.12500,coeff[1][1] = -0.04166667;
		coeff[2][0] = 1.1718750,coeff[2][1] = -0.065104167,
		coeff[2][2] = 0.0046875;
		coeff[3][0] = 1.196289,coeff[3][1] = -0.0797526,
		coeff[3][2] = 0.009570313,coeff[3][3] = -0.0006975447;
		coeff[4][0] = 1.211243,coeff[4][1] = -0.08972168 ,
		coeff[4][2] = 0.01384277,coeff[4][3] = -0.00176566,
		coeff[4][4] = 0.0001186795;
		coeff[5][0] = 1.22133636,coeff[5][1] = -0.09693146,
		coeff[5][2] = 0.01744766,coeff[5][3] = -0.00296729,
		coeff[5][4] = 0.00035901;coeff[5][5] = -0.00002185;

	}
	else
	{
		printf("Input fd order is not in the range of 1~6!\n"),exit(0);
	}

}
/* initialize the MPL boundary coefficients */
static void pml_init (int nx, int nz, float dx, float dz, float dt,
	float **dvv, float **od, int verbose,int order)
/*******************************************************************
*********************************************************************/
{
   	int ix, iz;
   	float sigma; 
	float K_MAX_PML;
	float ALPHA_MAX_PML;
	float NPOWER;
	float Rcoef,cp;
	float d0_x,d0_z;
   	float thickness_PML_x,thickness_PML_z;
	float xoriginleft,xoriginright;
	float zorigintop,zoriginbottom;
	float xval,zval;
	float abscissa_in_PML;
	float abscissa_normalized;
    

    	/* thickness of the PML layer in meters */
    	thickness_PML_x = dx*pml_thick;
    	thickness_PML_z = dz*pml_thick;


	/* parameter*/
	K_MAX_PML= 1.0 ;
	//ALPHA_MAX_PML = 2.0*PI*(f0/2.0);
	ALPHA_MAX_PML = 2.0*PI*(20./2.0);
	NPOWER = 2.;
	Rcoef = 0.0001;
	cp = 8000.;

	/* compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf*/
	d0_x = - (NPOWER + 1.) * cp * log(Rcoef) / (2.0 * thickness_PML_x);
  	d0_z = - (NPOWER + 1.) * cp * log(Rcoef) / (2.0 * thickness_PML_z);


	/* Allocate arrays for CPML boundary on front and back */
   	d_x = alloc1float (pml_thickness+nx);
   	K_x = alloc1float (pml_thickness+nx);
	a_x = alloc1float (pml_thickness+nx);
	b_x = alloc1float (pml_thickness+nx);	
	alpha_x = alloc1float (pml_thickness+nx);
	memset(d_x,0,sizeof(float)*(pml_thickness+nx));
	memset(a_x,0,sizeof(float)*(pml_thickness+nx));
	memset(alpha_x,0,sizeof(float)*(pml_thickness+nx));
	for(ix=0;ix<nx+pml_thickness;ix++)
		K_x[ix] = 1.0;	

	d_x_half = alloc1float (pml_thickness+nx);
	K_x_half = alloc1float (pml_thickness+nx);
	a_x_half = alloc1float (pml_thickness+nx);
	b_x_half = alloc1float (pml_thickness+nx);
	alpha_x_half = alloc1float (pml_thickness+nx);
	memset(d_x_half,0,sizeof(float)*(pml_thickness+nx));	
	memset(a_x_half,0,sizeof(float)*(pml_thickness+nx));
	memset(alpha_x_half,0,sizeof(float)*(pml_thickness+nx));
	for(ix=0;ix<nx+pml_thickness;ix++)
		K_x_half[ix] = 1.0;


  	/* Allocate arrays for CPML boundary on upper and bottom */
	d_z = alloc1float (pml_thickness+nz);
   	K_z = alloc1float (pml_thickness+nz);
	a_z = alloc1float (pml_thickness+nz);
	b_z = alloc1float (pml_thickness+nz);	
	alpha_z = alloc1float (pml_thickness+nz);
	memset(d_z,0,sizeof(float)*(pml_thickness+nz));
	memset(a_z,0,sizeof(float)*(pml_thickness+nz));
	memset(alpha_z,0,sizeof(float)*(pml_thickness+nz));
	for(iz=0;iz<nz+pml_thickness;iz++)
		K_z[iz] = 1.0;

	d_z_half = alloc1float (pml_thickness+nz);
	K_z_half = alloc1float (pml_thickness+nz);
	a_z_half = alloc1float (pml_thickness+nz);
	b_z_half = alloc1float (pml_thickness+nz);
	alpha_z_half = alloc1float (pml_thickness+nz);
	memset(d_z_half,0,sizeof(float)*(pml_thickness+nz));
	memset(a_z_half,0,sizeof(float)*(pml_thickness+nz));
	memset(alpha_z_half,0,sizeof(float)*(pml_thickness+nz));
	for(iz=0;iz<nz+pml_thickness;iz++)
		K_z_half[iz] = 1.0;	



	/*************************** damping in the X direction *********************************/
	/*origin of the PML layer (position of front and back edge minus thickness, in meters)*/
  	xoriginleft = thickness_PML_x;
  	xoriginright = (pml_thick+nx-1)*dx;
 
	for(ix=0;ix<pml_thickness+nx;ix++)
	{

		xval = dx*ix;
		
		/* for front*/
		abscissa_in_PML = xoriginleft - xval;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at the grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_x;			
			d_x[ix] = d0_x * pow(abscissa_normalized,NPOWER);
			K_x[ix] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
       	 	alpha_x[ix] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;

			
		}
		abscissa_in_PML = xoriginleft - (xval + 0.5*dx);
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at half grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_x;
			d_x_half[ix] = d0_x * pow(abscissa_normalized,NPOWER);
			K_x_half[ix] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
			alpha_x_half[ix] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}
		
		
		/* for back */
		abscissa_in_PML = xval - xoriginright;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at the grid points*/
			abscissa_in_PML = xval - xoriginright;
			abscissa_normalized = abscissa_in_PML / thickness_PML_x;
			d_x[ix] = d0_x * pow(abscissa_normalized,NPOWER);
			K_x[ix] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
			alpha_x[ix] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}
		abscissa_in_PML = (xval + 0.5*dx)-xoriginright;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at half grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_x;
			d_x_half[ix] = d0_x * pow(abscissa_normalized,NPOWER);
			K_x_half[ix] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
			alpha_x_half[ix] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}


		/* define b and a*/
		if(alpha_x[ix]<0.)  alpha_x[ix] = 0.;
		if(alpha_x_half[ix]<0.)  alpha_x_half[ix] = 0.;
		b_x[ix] = exp(- (d_x[ix] / K_x[ix] + alpha_x[ix]) * dt);
   		b_x_half[ix] = exp(- (d_x_half[ix] / K_x_half[ix] + alpha_x_half[ix]) * dt);

		if(fabs(d_x[ix]) > 1.0e-6) 
				a_x[ix] = d_x[ix] * (b_x[ix] - 1.0) / 
					(K_x[ix] * (d_x[ix] + K_x[ix] * alpha_x[ix]));
    		if(fabs(d_x_half[ix]) > 1.0e-6) 
				a_x_half[ix] = d_x_half[ix] *
						(b_x_half[ix] - 1.0) / 
						(K_x_half[ix] * (d_x_half[ix] + 
						K_x_half[ix] * alpha_x_half[ix]));
	}


	/* damping in the Z direction */
	/*origin of the PML layer (position of front and back edge minus thickness, in meters)*/
	zorigintop = thickness_PML_z;
  	zoriginbottom = (pml_thick+nz-1)*dz;
	for(iz=0;iz<nz+pml_thickness;iz++)
	{

		zval = dz*iz;

		/* for top*/
		abscissa_in_PML = zorigintop - zval;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at the grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_z;
			d_z[iz] = d0_z * pow(abscissa_normalized,NPOWER);
			K_z[iz] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
       	 	alpha_z[iz] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}
		abscissa_in_PML = zorigintop - (zval + 0.5*dz);
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at half grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_z;
			d_z_half[iz] = d0_z * pow(abscissa_normalized,NPOWER);
			K_z_half[iz] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
			alpha_z_half[iz] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;	
		}

		/* for bottom*/
		abscissa_in_PML = zval - zoriginbottom;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at the grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_z;
			d_z[iz] = d0_z * pow(abscissa_normalized,NPOWER);
			K_z[iz] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
			alpha_z[iz] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}
		abscissa_in_PML = zval + 0.5*dz - zoriginbottom;
		if(abscissa_in_PML>=0.)
		{
			/* define damping profile at half the grid points*/
			abscissa_normalized = abscissa_in_PML / thickness_PML_z;
			d_z_half[iz] = d0_z * pow(abscissa_normalized,NPOWER);
			K_z_half[iz] = 1.0 + (K_MAX_PML - 1.0) * pow(abscissa_normalized,NPOWER);
        		alpha_z_half[iz] = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML;
		}

		/* define b and a */
		if(alpha_z[iz]<0.)  alpha_z[iz] = 0.;
		if(alpha_z_half[iz]<0.)  alpha_z_half[iz] = 0.;
		b_z[iz] = exp(- (d_z[iz] / K_z[iz] + alpha_z[iz]) * dt);
    		b_z_half[iz] = exp(- (d_z_half[iz] / K_z_half[iz] + alpha_z_half[iz]) * dt);

		if(fabs(d_z[iz]) > 1.0e-6) 
			a_z[iz] = d_z[iz] * (b_z[iz] - 1.0) / 
					(K_z[iz] * (d_z[iz] + K_z[iz] * alpha_z[iz]));
    		if(fabs(d_z_half[iz]) > 1.0e-6) 
			a_z_half[iz] = d_z_half[iz] *
					(b_z_half[iz] - 1.0) / 
					(K_z_half[iz] * (d_z_half[iz] + 
					K_z_half[iz] * alpha_z_half[iz]));
	}

}

/* load the velocity and density froom the inout model to the current calculation mocel*/
void load_dvvc_odc(int is,int ds,float **dvv,float **epsilon,float **delta,float **od,float **dvvc,float **epsilonc,float **deltac,
	float **odc,int nx,int nz,int nxc,int nzc,int *hele,float *ele,
	int pml_thick,int pml_thickness,int itype_acq)
/***************************************************************************************
****************************************************************************************/
{
	register int ix,iz;
	int ix_beg,ix_end,iz_beg,iz_end;
	
	if(itype_acq==0){
	  ix_beg = 0;
	  ix_end = nx;
	  iz_beg = 0;
	  iz_end = nz;
        }else{
	  ix_beg = (is-1)*ds;
	  ix_end = nxc+ix_beg;
	  iz_beg = 0;
	  iz_end = nzc;
        }	
	printf("is=%d ix_end=%d \n",is,ix_end);

	if(ix_beg<0 || ix_end>nx || iz_beg<0 || iz_beg>nz)
		err("calculating break bounds");
	
	for(ix=ix_beg;ix<ix_end;++ix)
		for(iz=iz_beg;iz<iz_end;++iz)
		{
			dvvc[pml_thick+ix-ix_beg][pml_thick+iz-iz_beg] = dvv[ix][iz];
			epsilonc[pml_thick+ix-ix_beg][pml_thick+iz-iz_beg] = epsilon[ix][iz];
			deltac[pml_thick+ix-ix_beg][pml_thick+iz-iz_beg] = delta[ix][iz];
			odc[pml_thick+ix-ix_beg][pml_thick+iz-iz_beg] = od[ix][iz];
		}
	
	for(ix=ix_beg;ix<ix_end;++ix)
		hele[ix-ix_beg] = ceil(ele[ix])+2;

	/* pad left side */
	for(ix=0;ix<pml_thick;++ix)
		for(iz=pml_thick;iz<pml_thick+nzc;++iz)
		{
			dvvc[ix][iz] = dvvc[pml_thick][iz];
			epsilonc[ix][iz] = epsilonc[pml_thick][iz];
			deltac[ix][iz] = deltac[pml_thick][iz];
			odc[ix][iz] = odc[pml_thick][iz];
		}
	/* pad right side */
	for(ix=pml_thick+nxc;ix<nxc+pml_thickness;++ix)
		for(iz=pml_thick;iz<pml_thick+nzc;++iz)
		{
			dvvc[ix][iz] = dvvc[pml_thick+nxc-1][iz];
			epsilonc[ix][iz] = epsilonc[pml_thick+nxc-1][iz];
			deltac[ix][iz] = deltac[pml_thick+nxc-1][iz];
			odc[ix][iz] = odc[pml_thick+nxc-1][iz];
		}
	/* pad upp side */
	for(ix=0;ix<nxc+pml_thickness;++ix)
		for(iz=0;iz<pml_thick;++iz)
		{
			dvvc[ix][iz] = dvvc[ix][pml_thick];
			epsilonc[ix][iz] = epsilonc[ix][pml_thick];
			deltac[ix][iz] = deltac[ix][pml_thick];
			odc[ix][iz] = odc[ix][pml_thick];
		}
	/* pad down side */
	for(ix=0;ix<nxc+pml_thickness;++ix)
		for(iz=pml_thick+nzc;iz<nzc+pml_thickness;++iz)
		{
			dvvc[ix][iz] = dvvc[ix][pml_thick+nzc-1];
			epsilonc[ix][iz] = epsilonc[ix][pml_thick+nzc-1];
			deltac[ix][iz] = deltac[ix][pml_thick+nzc-1];
			odc[ix][iz] = odc[ix][pml_thick+nzc-1];
		}
}

/* prototype of subroutine used internally */
static float ricker (float t, float fpeak);

void ptsrc (int ixs, int izs,int nx, float dx, float fx,int nz, float dz,
	 float fz,float dt, float t, float fmax, float fpeak, float *tdelay,
 	 float **s,int pml_thick,int pml_thickness)
/*****************************************************************************
ptsrc - update source pressure function for a point source
******************************************************************************
Input:
xs		x coordinate of point source
zs		z coordinate of point source
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
dt		time step (ignored)
t		time at which to compute source function
fmax		maximum frequency (ignored)
fpeak		peak frequency

Output:
tdelay		time delay of beginning of source function
s		array[nx][nz] of source pressure at time t+dt
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/01/90
******************************************************************************/
{
	int ix,iz;
	float ts,xn,zn,td;
	
	/* zero source array */
	for (ix=0; ix<nx+pml_thickness; ++ix)
		for (iz=0; iz<nz+pml_thickness; ++iz)
			s[ix][iz] = 0.0 * dt*fmax;
	
	/* compute time-dependent part of source function */
	/* fpeak = 0.5*fmax;  this is now getparred */
	td = 1.0/fpeak;
	if (t>2.0*(td)) return;
	ts = ricker(t-td,fpeak);


	/* let source contribute within limited distance */
	#pragma omp parallel for shared(pml_thick,ixs,izs,nx,nz,s,ts,dx,dz) private(ix,iz,xn,zn) 
	for (ix=MAX(pml_thick,pml_thick+ixs-3); ix<=MIN(pml_thick+nx-1,pml_thick+ixs+3); ++ix) 
	{
	   for (iz=MAX(pml_thick,pml_thick+izs-3);iz<=MIN(pml_thick+nz-1,pml_thick+izs+3); ++iz) 
	   {
		xn = (ix-ixs-pml_thick)*dx/sqrt(dx*dx+dz*dz);
		zn = (iz-izs-pml_thick)*dz/sqrt(dx*dx+dz*dz);
		s[ix][iz] = ts*exp(-xn*xn-zn*zn);
	   }
	}


	
	//s[ixs+pml_thick][izs+pml_thick] = ts;

	*tdelay = td;
}

static float ricker (float t, float fpeak)
/*****************************************************************************
ricker - Compute Ricker wavelet as a function of time
******************************************************************************
Input:
t		time at which to evaluate Ricker wavelet
fpeak		peak (dominant) frequency of wavelet
******************************************************************************
Notes:
The amplitude of the Ricker wavelet at a frequency of 2.5*fpeak is 
approximately 4 percent of that at the dominant frequency fpeak.
The Ricker wavelet effectively begins at time t = -1.0/fpeak.  Therefore,
for practical purposes, a causal wavelet may be obtained by a time delay
of 1.0/fpeak.
The Ricker wavelet has the shape of the second derivative of a Gaussian.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/29/90
******************************************************************************/
{
	float x,xx;
	
	x = PI*fpeak*t;
	xx = x*x;
	/* return (-6.0+24.0*xx-8.0*xx*xx)*exp(-xx); */
	/* return PI*fpeak*(4.0*xx*x-6.0*x)*exp(-xx); */
	return exp(-xx)*(1.0-2.0*xx);
}


/* 2D finite differencing subroutine */
void ptvel(int nx,int nz,float dx,float dz,float dt,float **p,float **vx,
	float**vz,float **dpdx,float**dpdz,float**od,int order,float**coeff,
	int pml_thick,int pml_thickness);


void ptprs(int nx,int nz,float dx,float dz,float dt,float**p,float**pp,float**pm,float**Ep,float**ud,float **E,float **u,float **v,float **dpdx,float**dpdz,float**vx,float**vz,float**dpdxx,float**dpdzz,float**dvv,float**epsilon,float**delta,float**s,int order,int ssize,float tol,float lambda,float gama,int iteration,float**coeff,int pml_thick,int pml_thickness);

void gausssmooth_2d(float**in,int n1,int n2, float r1);

void tstep2 (int nx, float dx, int nz, float dz, float dt,float **dvv,float **epsilon,float **delta,float **od, float **s,float **p,float **Ep,float **ud,float **E,float **u,float **v,float **pp,float **pm,float **dpdx,float **dpdz,float **vx,float **vz,float **dpdxx,
float **dpdzz,int order,int ssize,float tol,float lambda,float gama,int iteration,float **coeff,int pml_thick,int pml_thickness)
/*****************************************************************************
One time step of FD solution (2nd order in space) to acoustic wave equation
tstep2(nxc,dx,nzc,dz,dt,dvvc,odc,s,p,dpdx,dpdz,vx,vz,dpdxx,dpdzz,order,
				coeff,pml_thick,pml_thickness)
******************************************************************************
Input:
nx		number of x samples
dx		x sampling interval
nz		number of z samples
dz		z sampling interval
dt		time step
dvv		array[nx][nz] of density*velocity^2
od		array[nx][nz] of 1/density (NULL for constant density=1.0)
s		array[nx][nz] of source pressure at time t+dt
pm		array[2][nx][nz] of pressure at time t-dt
p		array[2][nx][nz] of pressure at time t

Output:
pp		array[nx][nz] of pressure at time t+dt
******************************************************************************
Notes:
This function is optimized for special cases of constant density=1 and/or
equal spatial sampling intervals dx=dz.  The slowest case is variable
density and dx!=dz.  The fastest case is density=1.0 (od==NULL) and dx==dz.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/13/90
******************************************************************************/
{
	/* convolve with finite-difference star (special cases for speed) */
//        ptvel(nx,nz,dx,dz,dt,p,vx,vz,dpdx,dpdz,od,order,coeff,pml_thick,pml_thickness);
	ptprs(nx,nz,dx,dz,dt,p,pp,pm,Ep,ud,E,u,v,dpdx,dpdz,vx,vz,dpdxx,dpdzz,dvv,epsilon,delta,s,order,ssize,tol,lambda,gama,iteration,coeff,pml_thick,pml_thickness);

}


#if 0
void ptvel(int nx,int nz,float dx,float dz,float dt,float **p,float **vx,
	float**vz,float **dpdx,float**dpdz,float**od,int order,float**coeff,
	int pml_thick,int pml_thickness)
{
	register int ix,iz,i,ir;
	float txscl,tzscl,epsilon;
	float ux,uz;
	
	txscl = 1/dx;
	tzscl = 1/dz;
	/* x*/
	#pragma omp parallel for default(shared) private(ix,iz,ux)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=0;iz<nz+pml_thickness;iz++)
	   {
		ux = 0.;
		ux =  coeff[3][0]*(p[ix+1][iz]-p[ix][iz])+
		      coeff[3][1]*(p[ix+2][iz]-p[ix-1][iz])+
			coeff[3][2]*(p[ix+3][iz]-p[ix-2][iz])+
			coeff[3][3]*(p[ix+4][iz]-p[ix-3][iz]);
		
                dpdx[ix][iz] = b_x[ix] * dpdx[ix][iz] + a_x[ix] * ux;
                ux = ux / K_x[ix] + dpdx[ix][iz];
		vx[ix][iz] = txscl*ux;

	   }	

	/* z*/
	#pragma omp parallel for default(shared) private(ix,iz,uz)
	for(ix=0;ix<nx+pml_thickness;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	  {
		uz = 0.;
		uz =  coeff[order-1][0]*(p[ix][iz+1]-p[ix][iz])+
			coeff[order-1][1]*(p[ix][iz+2]-p[ix][iz-1])+
			coeff[order-1][2]*(p[ix][iz+3]-p[ix][iz-2])+
			coeff[order-1][3]*(p[ix][iz+4]-p[ix][iz-3]);
		
                dpdz[ix][iz] = b_z[iz] * dpdz[ix][iz] + a_z[iz] * uz;
                uz = uz / K_z[iz] + dpdz[ix][iz];
		vz[ix][iz] = tzscl*uz;
	  }

}

#endif

void ptprs(int nx,int nz,float dx,float dz,float dt,float**p,float**pp,float**pm,float**Ep,float**ud,float **E,float **u,float **v,float**dpdx,float**dpdz,float**vx,float**vz,float **dpdxx,float**dpdzz,
float**dvv,float**epsilon,float**delta,float**s,int order,int ssize,float tol,float lambda,float gama,int iteration,float**coeff,int pml_thick,int pml_thickness)
{
	register int i,j,ix,iz,m,n;
	float tscl2,txscl,tzscl;
	float pc,sum;
	float maxl,Emax,abn;
	float ux,uz,ud1,uvm,uit,vit;
	float **Ex,**Ez,**Et,**um,**vm,**vxx,**vzz;
	float *uinline,*udepth;
	float *vinline,*vdepth;

	Ex = alloc2float(nz+pml_thickness, nx+pml_thickness);
	Ez = alloc2float(nz+pml_thickness, nx+pml_thickness);
	Et = alloc2float(nz+pml_thickness, nx+pml_thickness);
	um = alloc2float(nz+pml_thickness, nx+pml_thickness);
	vm = alloc2float(nz+pml_thickness, nx+pml_thickness);
	vxx = alloc2float(nz+pml_thickness, nx+pml_thickness);
	vzz = alloc2float(nz+pml_thickness, nx+pml_thickness);
	uinline = alloc1float(nz+pml_thickness);
	udepth = alloc1float(nx+pml_thickness);
	vinline = alloc1float(nz+pml_thickness);
	vdepth = alloc1float(nx+pml_thickness);
        memset(Ex[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(Ez[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(Et[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(um[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(vm[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(u[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(v[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(vxx[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));
        memset(vzz[0],0,sizeof(float)*(nx+pml_thickness)*(nz+pml_thickness));

	tscl2 = dt*dt;
	txscl = 1/dx;
	tzscl = 1/dz;
//velocity
      #pragma omp parallel for default(shared) private(ix,iz,ud1) 
	for(ix=0;ix<nx+pml_thickness;ix++)
	   for(iz=0;iz<nz+pml_thickness;iz++)
	   {
		ud1=(1 + 2*epsilon[ix][iz]) + 1;
          	ud[ix][iz]=1+sqrt(1-8*(epsilon[ix][iz]-delta[ix][iz])/(ud1*ud1));

	   }

	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	   {
		ux =  coeff[3][0]*(p[ix+1][iz]-p[ix][iz])+
		      coeff[3][1]*(p[ix+2][iz]-p[ix-1][iz])+
			coeff[3][2]*(p[ix+3][iz]-p[ix-2][iz])+
			coeff[3][3]*(p[ix+4][iz]-p[ix-3][iz]);
		
		uz =  coeff[order-1][0]*(p[ix][iz+1]-p[ix][iz])+
			coeff[order-1][1]*(p[ix][iz+2]-p[ix][iz-1])+
			coeff[order-1][2]*(p[ix][iz+3]-p[ix][iz-2])+
			coeff[order-1][3]*(p[ix][iz+4]-p[ix][iz-3]);

                dpdx[ix][iz] = b_x[ix] * dpdx[ix][iz] + a_x[ix] * ux;
                ux = ux / K_x[ix] + dpdx[ix][iz];
		vx[ix][iz] = txscl*ux;

                dpdz[ix][iz] = b_z[iz] * dpdz[ix][iz] + a_z[iz] * uz;
                uz = uz / K_z[iz] + dpdz[ix][iz];
		vz[ix][iz] = tzscl*uz;

	   }
      #pragma omp parallel for default(shared) private(ix,iz,ux,uz,pc) 
	for(ix=order;ix<nx+pml_thickness-order+1;++ix)
	   for(iz=order;iz<nz+pml_thickness-order+1;++iz)
	   {
		ux = coeff[order-1][0]*(vx[ix][iz]-vx[ix-1][iz])+
			coeff[order-1][1]*(vx[ix+1][iz]-vx[ix-2][iz])+
			coeff[order-1][2]*(vx[ix+2][iz]-vx[ix-3][iz])+
			coeff[order-1][3]*(vx[ix+3][iz]-vx[ix-4][iz]);		

		uz = coeff[order-1][0]*(vz[ix][iz]-vz[ix][iz-1])+
			coeff[order-1][1]*(vz[ix][iz+1]-vz[ix][iz-2])+
			coeff[order-1][2]*(vz[ix][iz+2]-vz[ix][iz-3])+
			coeff[order-1][3]*(vz[ix][iz+3]-vz[ix][iz-4]);
		
                dpdxx[ix][iz] = b_x[ix] * dpdxx[ix][iz] + a_x[ix] * ux;
                dpdzz[ix][iz] = b_z[iz] * dpdzz[ix][iz] + a_z[iz] * uz;

                ux = ux / K_x[ix] + dpdxx[ix][iz];
                uz = uz / K_z[iz] + dpdzz[ix][iz];
		vxx[ix][iz] = txscl*ux;
		vzz[ix][iz] = tzscl*uz;

	        pc = 0.5*tscl2*dvv[ix][iz]*((1+2*epsilon[ix][iz])*vxx[ix][iz]+vzz[ix][iz])*ud[ix][iz]+2*p[ix][iz]-pm[ix][iz];

                E[ix][iz]=pc*pc;

	   }
		Emax=0;
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	   {
		
		if(E[ix][iz]>Emax)
		{
			Emax=E[ix][iz];
		}
	   }

	if(Emax>1e-20)
	{
	#pragma omp parallel for default(shared) private(ix,iz)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	   {
		
		E[ix][iz]=E[ix][iz]/Emax;
	   }
	}
	else
	{
	Emax=Emax+1e-20;
	#pragma omp parallel for default(shared) private(ix,iz)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	   {
		
		E[ix][iz]=E[ix][iz]/Emax;
	   }
	}
	
//#if 0
//optical flow
	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
	   {
		
		ux = -0.08333333*(E[ix+2][iz]-E[ix-2][iz])+0.66666667*(E[ix+1][iz]-E[ix-1][iz]);
		uz = -0.08333333*(E[ix][iz+2]-E[ix][iz-2])+0.66666667*(E[ix][iz+1]-E[ix][iz-1]);

		Ex[ix][iz] = txscl*ux;
		Ez[ix][iz] = tzscl*uz;

         	Et[ix][iz]= (E[ix][iz]-Ep[ix][iz])/dt;
		Ep[ix][iz]= E[ix][iz];
           }


	int ii=0;
	float uk;	
	float vk;	

	while(ii < iteration)
	{
		ii++;

	#pragma omp parallel for default(shared) private(ix,iz,uk,vk)
	for(ix=order-1;ix<nx+pml_thickness-order;ix++)
	   for(iz=order-1;iz<nz+pml_thickness-order;iz++)
       {               
	       uk=u[ix][iz];
	       vk=v[ix][iz];		

               um[ix][iz]=(u[ix+1][iz]+u[ix][iz]+u[ix-1][iz]+u[ix+1][iz+1]+u[ix][iz+1]+u[ix-1][iz+1]+u[ix+1][iz-1]+u[ix][iz-1]+u[ix-1][iz-1]+u[ix][iz-1]+u[ix][iz+1]+u[ix+1][iz]+u[ix-1][iz]-u[ix][iz])/12;
               vm[ix][iz]=(v[ix+1][iz]+v[ix][iz]+v[ix-1][iz]+v[ix+1][iz+1]+v[ix][iz+1]+v[ix-1][iz+1]+v[ix+1][iz-1]+v[ix][iz-1]+v[ix-1][iz-1]+v[ix][iz-1]+v[ix][iz+1]+v[ix+1][iz]+v[ix-1][iz]-v[ix][iz])/12;

               u[ix][iz]=um[ix][iz]-(Ex[ix][iz]*(Ex[ix][iz]*um[ix][iz]+Ez[ix][iz]*vm[ix][iz]+Et[ix][iz]))/(lambda*lambda+Ex[ix][iz]*Ex[ix][iz]+Ez[ix][iz]*Ez[ix][iz]);
               v[ix][iz]=vm[ix][iz]-(Ez[ix][iz]*(Ex[ix][iz]*um[ix][iz]+Ez[ix][iz]*vm[ix][iz]+Et[ix][iz]))/(lambda*lambda+Ex[ix][iz]*Ex[ix][iz]+Ez[ix][iz]*Ez[ix][iz]);

        }

		gausssmooth_2d(u,nz+pml_thickness,nx+pml_thickness,ssize);
		gausssmooth_2d(v,nz+pml_thickness,nx+pml_thickness,ssize);

//#if 0
	}
     

	//#pragma omp parallel for default(shared) private(ix,iz,uvm)
	for(ix=order;ix<nx+pml_thickness-order+1;++ix)
	   for(iz=order;iz<nz+pml_thickness-order+1;++iz)
	   {
		
                uvm=u[ix][iz]*u[ix][iz]+v[ix][iz]*v[ix][iz];
		if(uvm>maxl)
		{
			maxl=uvm;
		}

           }
	//printf("maxuv=%f\n",maxl);

//#endif
//#if 0
      #pragma omp parallel for default(shared) private(ix,iz,uit,vit,ud1,uvm) 
	//for(ix=order;ix<nx+pml_thickness-order+1;++ix)
	//   for(iz=order;iz<nz+pml_thickness-order+1;++iz)
	for(ix=pml_thick;ix<nx+pml_thick;++ix)
	   for(iz=pml_thick;iz<nz+pml_thick;++iz)
	   {
                uvm=u[ix][iz]*u[ix][iz]+v[ix][iz]*v[ix][iz];
		if (uvm>maxl*gama)
		{
             	uit=u[ix][iz]*u[ix][iz]/uvm;
             	vit=v[ix][iz]*v[ix][iz]/uvm;     
		ud1=(1 + 2*epsilon[ix][iz])*uit + vit;
          	ud[ix][iz]=1+sqrt(1-8*(epsilon[ix][iz]-delta[ix][iz])*uit*vit/(ud1*ud1));
//		printf("uit=%.20f vit=%.20f \n",uit,vit);
		}

	   }

//#endif
//#if 0
	abn=1.0/(pml_thick*pml_thick);
//A area
	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=pml_thick;ix<nx+pml_thick;ix++)
	   for(iz=order;iz<pml_thick;iz++)
	  { 
		ux=iz*iz*abn;
		uz=ux*(ud[ix][pml_thick]-ud[ix][iz])+ud[ix][iz];
		ud[ix][iz]=uz;
	  }
//C area
	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=pml_thick;ix<nx+pml_thick;ix++)
	   for(iz=nz+pml_thick;iz<pml_thickness+nz-order;iz++)
	   {
		ux=(nz+pml_thickness-iz)*(nz+pml_thickness-iz)*abn;
		uz=ux*(ud[ix][nz+pml_thick-1]-ud[ix][iz])+ud[ix][iz];
		ud[ix][iz]=uz;
	   }

//B area

	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=order;ix<pml_thick;ix++)
	   for(iz=order;iz<pml_thickness+nz;iz++)
	   {
		ux=ix*ix*abn;
		uz=ux*(ud[pml_thick][iz]-ud[ix][iz])+ud[ix][iz];
		ud[ix][iz]=uz;
	   }
//D area

	#pragma omp parallel for default(shared) private(ix,iz,ux,uz)
	for(ix=nx+pml_thick;ix<nx+pml_thickness-order;ix++)
	   for(iz=order;iz<pml_thickness+nz;iz++)
	   {
		ux=(nx+pml_thickness-ix)*(nx+pml_thickness-ix)*abn;
		uz=ux*(ud[nx+pml_thick-1][iz]-ud[ix][iz])+ud[ix][iz];
		ud[ix][iz]=uz;
	   }
//#endif

      #pragma omp parallel for default(shared) private(ix,iz,pc) 
	for(ix=order;ix<nx+pml_thickness-order+1;++ix)
	   for(iz=order;iz<nz+pml_thickness-order+1;++iz)
	   {
	        pc = 0.5*tscl2*dvv[ix][iz]*((1+2*epsilon[ix][iz])*vxx[ix][iz]+vzz[ix][iz])*ud[ix][iz]+2*p[ix][iz]-pm[ix][iz];

	        pp[ix][iz] = pc+s[ix][iz];
                pm[ix][iz] = p[ix][iz];
                p[ix][iz] = pp[ix][iz];
	   }
	free2float(Et);
	free2float(Ex);
	free2float(Ez);
	free2float(um);
	free2float(vm);
	free2float(vxx);
	free2float(vzz);
	free(uinline);
	free(udepth);
	free(vinline);
	free(vdepth);

}

/* write seismic record in su format*/
int getoutnum(int n,int dn);
static void write_su(int is,int ds,int xsn,int nx,int ndtr,float fx, 
	int nt,float dt,float dtout,float dx,float **hs,FILE*fp_out)
{
	Head head={0};
	register int itr,it;
	int ndtt,nxout,ntout;
	float dtr;
	int cdp_shot,cdp_rev;
	unsigned long offset;
	int nit=0;

	dtr = dx*ndtr;
	ndtt = (int)(dtout/dt+0.5);
	nxout = getoutnum(nx,ndtr);
	ntout = getoutnum(nt,ndtt);
	cdp_shot = xsn+(is-1)*ds;

	offset = ((is-1)*nxout)*(ntout*sizeof(float)+sizeof(Head));
	fseek(fp_out,offset,SEEK_SET);

	for(itr=0;itr<nxout;itr++)
	{

		head.tracl = (is-1)*(nxout)+itr+1;
		head.tracr = (is-1)*(nxout)+itr+1;
		head.fldr = is;
		head.trid = 1;
		head.tracf = itr+1;
		cdp_rev = fx+itr*ndtr+1+(is-1)*ds;
		head.cdp = (int)((cdp_shot+cdp_rev)/2.0+0.5);
		head.offset = dx*(cdp_rev-cdp_shot);
		head.sx = (int)(fx+(cdp_shot-1)*dx);
		head.gx = (int)(fx+(cdp_rev-1)*dx);
		head.sy = 0;
		head.gy = 0;
		head.ns = (short)ntout;
		head.dt = (short)(dtout*1000000);

		fwrite(&head,sizeof(Head),1,fp_out);
		nit = 0;
		for(it=0;it<nt;it+=ndtt)
		{	
			nit++;	
			fwrite(&hs[itr*ndtr][it],sizeof(float),1,fp_out);
		}
		if(nit!=ntout) printf("write ns error!\n"),exit(0);
	}
}

int getoutnum(int n,int dn) 
{
	int nout;
	
	nout = n/dn;
	if((n-dn*(nout-1)-1)==dn)
	nout +=1; 	
	return nout;
}

void gausssmooth_2d(float**in,int n1,int n2, float r1)
{
    int ix,iz,i1,j1;
    float *mat1,*mat2,*kernel;
    float sigma,s,tmp,sum,conv;
    int nw,hfs;
    /* define filter size */
    nw = round(r1);
    if (nw==0)	return;
	if (!(nw%2)) {
		if (nw > 0)
			nw += 1;
		else
			nw -= 1;
	}
    nw = abs(nw);
    /* parameters */
    hfs = abs(nw)/2;
	sigma = hfs/2.0;
	s = 2.0*sigma*sigma;

    kernel = alloc1float(2*hfs+1);

    /* create filter kernel */
    sum = 0.0;
	for (i1=0;i1<2*hfs+1;i1++)
	{
		tmp = 1.0*(i1-hfs);
	    kernel[i1] = exp(-(tmp*tmp)/s);
		sum += kernel[i1];
        }
    /* normalize kernel */
    for (i1=0;i1<2*hfs+1;i1++)
    {
        kernel[i1] /= sum;
    }

    #pragma omp parallel for default(shared) private(ix,mat1,i1,j1,conv)
    for(ix=0;ix<n2;ix++)	
    {
    /* copy input to mat */
    mat1 = alloc1float(n2+2*hfs);
    for(i1=0;i1<n1;i1++)
    {
    	mat1[i1+hfs] = in[ix][i1];
    }
    /* extend boundary */
    for(i1=0;i1<hfs;i1++)
    {
        mat1[i1] = mat1[hfs];
        mat1[i1+n1+hfs] = mat1[n1+hfs-1];
    }
    /* apply Gaussian filter */
    //#pragma omp parallel for default(shared) private(i1,j1,conv)
    for (i1=hfs;i1<n1+hfs;i1++)
    {
        /* loop over kernel*/
        conv = 0.0;
	  	for (j1=0;j1<2*hfs+1;j1++)
	  	{
	         conv += mat1[i1+j1-hfs]*kernel[j1];
                }
        /* output of filtered gradient */
        in[ix][i1-hfs] = conv;
    }
   free1float(mat1);
    }
	

    #pragma omp parallel for default(shared) private(iz,mat2,i1,j1,conv)
    for(iz=0;iz<n1;iz++)	
    {

    mat2 = alloc1float(n2+2*hfs);
    /* copy input to mat */
    for(i1=0;i1<n2;i1++)
    {
    	mat2[i1+hfs] = in[i1][iz];
    }
    /* extend boundary */
    for(i1=0;i1<hfs;i1++)
    {
        mat2[i1] = mat2[hfs];
        mat2[i1+n2+hfs] = mat2[n2+hfs-1];
    }
    /* apply Gaussian filter */
    for (i1=hfs;i1<n2+hfs;i1++)
    {
        /* loop over kernel*/
        conv = 0.0;
	  	for (j1=0;j1<2*hfs+1;j1++)
	  	{
	         conv += mat2[i1+j1-hfs]*kernel[j1];
                }
        /* output of filtered gradient */
        in[i1-hfs][iz] = conv;
    }
    free1float(mat2);
    }

   free1float(kernel);
}
