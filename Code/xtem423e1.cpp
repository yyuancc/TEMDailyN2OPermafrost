
/* **************************************************************
*****************************************************************
XTEM423E1.CPP - Terrestrial Ecosystem Model Version 4.2
*****************************************************************

Modifications:

19991028 - DWK added bug fixes found and implemented by Gregg
           Christopher and Jim Long
20000616 - DWK eliminates the global variables tempred[][][] and atmspred[][]
20000616 - DWK eliminates the task of writing TEM output files from
           extrapolate()
20000616 - DWK adds the global constant double ZERO
20010314 - JSC adds transient solar radiation as input possibility
20010418-  Q.Z. addition of soil thermal model
20020202 - DWK changed include from xtem423e.cpp to xtem423e1.cpp
20220404 - YYuan change to MPI
*****************************************************************
************************************************************** */

using namespace std;

#include<mpi.h>  //added by YYuan mpi
#include<stdio.h>
#include<iostream>
#include<fstream>
#include <fcntl.h>
#include<iomanip>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include<time.h>
//#include <conio.h>
#include <values.h>


const int CYCLE = 12;
const int MAXRTIME = 600;

const int MAXPRED = 92+4+10+8+1+6+2; // changed from 85 to 92 by DWK on 20000202, QZ, plus 4+10 //add +8 for N cycle by YYL //add +1 for DOC by YYL//added 6 by YYuan

const int MAXGRID = 1;
// const int MAXNUMROW = 360; // commented out by DWK on 20000202

const double MISSING = -999999.9;
const double ZERO = 0.000000;  // added by DWK on 20000615
long int kswitch; //added for soil thermal model


struct Temflags{
int stmflg;   // added for soil thermal model
int equil;
int RTIME;

int spinflag;
int numspin;
int spintime;
int ispinout;
int totsptime;
int transtime;
//int tskburden; //added by YYuan mpi

int atmstotyr;
int atmsflag;
int atmsoutfg;
int temflag;
int stateflag;

int cldflag;
int natmspred;
int ntempred;
int totpred;
char predmap[MAXPRED][9];
int kdinflg;
int kdoutflg;
}tmflgs;



ofstream flog1;
ifstream fpara;


/* *************************************************************
		 Function Declarations
************************************************************* */

int askpred(char pvarname[MAXPRED][9], int spred);
void extrapolate(void);
void startclm(void);
void starttem(void);
void subndtem(void);  //YY added for mpi
void fptname(char* fnmtemfl,int rank);  //YY added for mpi
// *************************************************************
int fatalerr;
int end1;
int icount;
int glob_count;

struct Temstac  //YY added for mpi
{
	int avlnflag;
	int nfeed;
	int initbase;
	int baseline;
	int moistlim;
	int strteq;
	int maxyears;
	int runsize;
	int maxnrun;
	int rheqflag;
	double wtol;
	double ctol;
	double ntol;
	int startyr;
	int endyr;
	int diffyr;
	
	double inittol;
	int maxit;
	long maxitmon;
	
}temstac;

#include "elmnt423.cpp"    //Elmnt Class
#include "telm423e1.cpp"     //TEMelmnt Class

Elmnt elmnt;
TEMelmnt telmnt[MAXGRID];

//YY added for mpi
struct Fnminout{
   char fnmlonlat[60];
   char fnmstxt[60];
   char fnmelev[60];
   char fnmtveg[60];
   char fnmclds[60];
   char fnmtair[60];
   char fnmprec[60];
   char fnmnirr[60];
   char fnmvap[60]; 
   char fnmpar[60];
   char fnmlulc[60];
   char fnmnpp[60];
   char fnmstate[MAXSTATE][60];
   char fnmclmpred[NUMATMS][60];
   char fnmtempred[MAXPRED][60];
   char fnmkdin[60];
   char fnmkdout[60];
   char fnmph[60];
   char fnmdens[60];
   char fnmpercn[60];   
 
}fnminout;
//end of edit

FILE* fvap;
FILE* flonlat;
FILE* fstxt;
FILE* felev;
FILE* fph;
FILE* fdens;
FILE* fpercn;
FILE* ftveg;
FILE* fclds;
FILE* ftair;
FILE* fprec;
FILE* fnirr;
FILE* fpar;
FILE* flulc;
FILE* fnpp;

FILE* fstate[MAXSTATE];
ofstream fclmpred[NUMATMS];
ofstream ftempred[MAXPRED];

FILE* fkdin;
ofstream fkdout;
KDdata kdparam;
ofstream otemls; //YY added


/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main(int argc,char* argv[])  //YY edit for mpi
//int main()
{
  	int i;
  	int dir_first;  
  	int dir_label;  
    printf("xtem192");
	int root=0;
	MPI::Init(argc,argv);
  	
	int tasks=MPI::COMM_WORLD.Get_size();
	int myrank=MPI::COMM_WORLD.Get_rank();
  		
	if(myrank==0)
	{ 
		if(argc < 2)
                {
                        flog1.open("TEM4.log");
                        fpara.open("run.go4");
                }
                else
                {
                        flog1.open(argv[2]);
                        fpara.open(argv[1]);
                }

		otemls.open("Temout.list",ofstream::out);
		if(!otemls.is_open())
		{
  			cout <<"Can not creat file temout.list"<<endl;
  			exit(-1);
		}	
		if(! fpara.is_open())
		{
			cout <<"Cannot open file run.go4 to specify data"<<endl;

  			exit(-1);
		}
		if (! flog1.is_open())
		{
			cout << endl << "Cannot open log file TEM4.log"<<endl;

  			exit(-1);
		}
 //end of editting
  // Open log file

  //flog1.open("tem4.log");

  // Open log file
  //flog1.open("tem4.log");

/* *************************************************************
  Run equilibrium simulation or transient simulation ?
************************************************************* */
//YY added for mpi
  std::cout <<endl <<"Please enter the lable of subdirectory you want to start with the simulations? "<< std::endl;
  fpara >> dir_first;

  cout <<"part-"<<dir_first<<"/"<<endl;
  flog1 <<endl <<"Please enter the lable of subdirectory you want to start with the simulations? "<< std::endl;
  flog1 <<"part-"<<dir_first<<"/"<<endl;
//end of edit

  tmflgs.equil = 1;
  std::cout << endl << "Do you want to run the model only for steady state conditions ? " << std::endl;
  std::cout << " Enter 0 for transient simulation" << std::endl;
  std::cout << " Enter 1 for steady state simulation" << std::endl;
  fpara >> tmflgs.equil;

  flog1 << std::endl << "Do you want to run the model only for steady state conditions ? " << std::endl;
  flog1 << " Enter 0 for transient simulation" << std::endl;
  flog1 << " Enter 1 for steady state simulation" << std::endl;
  flog1 << "equil = " << tmflgs.equil << std::endl << std::endl;

  tmflgs.RTIME = 1;
  tmflgs.ispinout = 0;
  if (tmflgs.equil == 0)
  {
    // Start transient TEM run from equilibrium conditions (i.e. spinflag == 0)
    //  or with a "spin up" period to remove potential artifacts associated with
    //  changing from equilibrium conditions to transient conditions from model
    //  results

    cout << "Do you want to start the transient with a spin up period? " << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    fpara  >> tmflgs.spinflag;

    flog1 << "Do you want to start the transient with a spin up period? " << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "spinflag = " << tmflgs.spinflag << endl << endl;

    tmflgs.totsptime = 0;
    if (tmflgs.spinflag == 1)
    {
      // Specify conditions for initializing TEM with a transient "spin up" period
      //   for a grid cell

      cout << "How many spins do you want in the spin up period? ";
      fpara >> tmflgs.numspin;
      flog1 << "How many spins do you want in the spin up period? " << endl;
      flog1 << "numspin = " << tmflgs.numspin << endl << endl;

      cout << "How many years per spin? ";
      fpara >> tmflgs.spintime;
      flog1 << "How many years per spin? " << endl;
      flog1 << "spintime = " << tmflgs.spintime << endl << endl;
      
      //YY edit for ispinout
			std::cout << "Do you want to output the spin up tem data? 1: yes/0: no"<< std::endl;
			fpara >> tmflgs.ispinout;
			flog1 <<"Do you want to output the spin up tem data? 1: yes/0: no  ";
			flog1 <<tmflgs.ispinout<<endl;
      
      tmflgs.totsptime = tmflgs.spintime * tmflgs.numspin;
      flog1 << "totsptime = " << tmflgs.totsptime << endl << endl;
      tmflgs.RTIME += tmflgs.totsptime;
    }

    // Specify conditions for the "non-spin up" part of the transient TEM run

    cout << endl << "How many years do you run for transient simulations ? " << endl;
    fpara >> tmflgs.transtime;

    flog1 << endl << "How many years do you run for transient simulations ? " << endl;
    flog1 << "transtime = " << tmflgs.transtime << endl << endl;

    tmflgs.RTIME += tmflgs.transtime;
    flog1 << "RTIME = " << tmflgs.RTIME << endl << endl;
  }

/* *************************************************************
What models do you want to run?  - Open associated files
************************************************************* */
 //  STM (soil thermal model)?

  cout << endl << "Do you want to run the SOIL THERMAL MODEL (STM) model for soil temperatures?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
  fpara >> tmflgs.stmflg;

  flog1 << endl << "Do you want to run the SOIL THERMAL MODEL (STM) model for soil temperatures?"<< endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << " stmflg = " << tmflgs.stmflg << endl << endl;
  if (tmflgs.stmflg == 1) {
// Get snow and soil layer dependent parameters, added for soil thermal model integration
  telmnt[0].tem.sthermal.getsnowecd(flog1);   //qtsp44a2.ecd
  telmnt[0].tem.sthermal.getsoillecd(flog1);  // qtla44a2.ecd
  telmnt[0].tem.sthermal.getsoiltecd(flog1);  // qtst44a2.ecd
  }


// POTSCLM model ?

  cout << endl << "Do you want to run the POTSCLM model for solar radiation variables?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
  fpara >> tmflgs.atmsflag;

  flog1 << endl << "Do you want to run the POTSCLM model for solar radiation variables?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << " atmsflag = " << tmflgs.atmsflag << endl << endl;

  tmflgs.totpred = tmflgs.natmspred = 0;
  if (tmflgs.atmsflag == 1) { 
  startclm(); 
  }

// TEM ?

  cout << endl << "Do you want to run the terrestrial ecosystem model (TEM)?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes" << endl;
  fpara >> tmflgs.temflag;

  flog1 << endl << "Do you want to run the terrestrial ecosystem model (TEM)?" << endl;
  flog1 << "  Enter 0 for No" << endl;
  flog1 << "  Enter 1 for Yes" << endl;
  flog1 << "temflag = " << tmflgs.temflag << endl << endl;

  tmflgs.ntempred = 0;
  if (tmflgs.temflag == 1) { starttem(); }

  // Use all records in GIS data sets (i.e. desired coverage)?

  elmnt.ask(flog1);
  fpara.close();
  flog1.close();
	otemls.close();
 }
  // Run Potsclm and TTEM modules over desired region (modules could potentially
  // run over several grid cells simultaneously with the use of extrapolate()
  
//YY added for mpi
	char strtmp[160],chtmp[10];

	//fptname(strtmp,myrank);
	//flog1.open(strtmp);

	flog1 << "Run tem on node " << myrank << std::endl;
	std::cout << "Run tem on node " << myrank << std::endl;
  
	MPI::COMM_WORLD.Bcast(&dir_first,1,MPI::INT,root);
	MPI::COMM_WORLD.Bcast(&tmflgs,sizeof(Temflags),MPI::BYTE,root);
	MPI::COMM_WORLD.Bcast(&(telmnt[0]),sizeof(TEMelmnt),MPI::BYTE,root);
	MPI::COMM_WORLD.Bcast(&elmnt,sizeof(Elmnt),MPI::BYTE,root);
	MPI::COMM_WORLD.Bcast(&fnminout,sizeof(Fnminout),MPI::BYTE,root);
	MPI::COMM_WORLD.Bcast(&temstac,sizeof(Temstac),MPI::BYTE,root);
	MPI::COMM_WORLD.Bcast(&(basicf), sizeof(Soilthermal), MPI::BYTE, root);

	dir_label = myrank + dir_first;	
	sprintf(strtmp,"./part-%d/%s-%d",dir_label,"TEM4.log",dir_label);
	flog1.open(strtmp);

	flog1 << "Run tem on node " << myrank << std::endl;
	std::cout << "Run tem on node " << myrank << std::endl;  
	if(myrank!=root)
	{
		subndtem();
	}
	
  
	if(myrank!=root)
	{
  		subndtem();
	}

	fptname(fnminout.fnmclds,dir_label);
	fclds=fopen(fnminout.fnmclds,"r");
	if(fclds==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmclds<<endl;
		exit(-1);
	}	

	if(telmnt[0].lonlatflag==0)
	{
		fptname(fnminout.fnmlonlat,dir_label);
		flonlat=fopen(fnminout.fnmlonlat,"r");
		if(flonlat==NULL)
		{
			cout <<"Can not open file "<<fnminout.fnmlonlat<<endl;
			exit(-1);
		}		
	}
	if(tmflgs.atmsoutfg==1)
	{
		for(int ii=0;ii<tmflgs.natmspred;ii++)
		{

			fptname(fnminout.fnmclmpred[i],dir_label);
			fclmpred[i].open(fnminout.fnmclmpred[ii],ofstream::out);
			if(!fclmpred[i].is_open())
			{
				cout <<"Failed to open file "<<fnminout.fnmclmpred[i]<<endl;
				exit(-1);
			}		
		}
	}	

	fptname(fnminout.fnmstxt,dir_label);
	fstxt=fopen(fnminout.fnmstxt,"r");
	if(fstxt==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmstxt<<endl;
		exit(-1);
	}	

	fptname(fnminout.fnmtveg,dir_label);
	ftveg=fopen(fnminout.fnmtveg,"r");
	if(ftveg==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmtveg<<endl;
		exit(-1);
	}	
	fptname(fnminout.fnmelev,dir_label);
	felev=fopen(fnminout.fnmelev,"r");
	if(felev==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmelev<<endl;
		exit(-1);
	}
 	fptname(fnminout.fnmph,dir_label);
	fph=fopen(fnminout.fnmph,"r");
	if(fph==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmph<<endl;
		exit(-1);
	}
 	fptname(fnminout.fnmdens,dir_label);
	fdens=fopen(fnminout.fnmdens,"r");
	if(fdens==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmdens<<endl;
		exit(-1);
	}
 	fptname(fnminout.fnmpercn,dir_label);
	fpercn=fopen(fnminout.fnmpercn,"r");
	if(fpercn==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmpercn<<endl;
		exit(-1);
	}
	fptname(fnminout.fnmtair,dir_label);
	ftair=fopen(fnminout.fnmtair,"r");
	if(ftair==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmtair<<endl;
		exit(-1);
	}		
	fptname(fnminout.fnmprec,dir_label);
	fprec=fopen(fnminout.fnmprec,"r");
	if(fprec==NULL)
	{
		cout <<"Can not open file "<<fnminout.fnmprec<<endl;
		exit(-1);
	}
 	fptname(fnminout.fnmvap,dir_label);
	fvap=fopen(fnminout.fnmvap,"r");
 	if (fvap==NULL) {
	cout <<"Can not open file "<<fnminout.fnmvap<<endl;
 		exit(-1); 
	}
	if(tmflgs.atmsflag==0)
	{
		fptname(fnminout.fnmnirr,dir_label);
		fnirr=fopen(fnminout.fnmnirr,"r");
		if(fnirr==NULL)
		{
			cout <<"Can not open file "<<fnminout.fnmnirr<<endl;
			exit(-1);
		}	
		fptname(fnminout.fnmpar,dir_label);
		fpar=fopen(fnminout.fnmpar,"r");
		if(fpar==NULL)
		{
			cout <<"Can not open file "<<fnminout.fnmpar<<endl;
			exit(-1);
		}	
	}	
	if(tmflgs.equil==0)
	{
		if(telmnt[0].tem.ag.tlulcflag==1)
		{
			fptname(fnminout.fnmlulc,dir_label);
			flulc=fopen(fnminout.fnmlulc,"r");
			if(flulc==NULL)
			{
				cout <<"Can not open file "<<fnminout.fnmlulc<<endl;
				exit(-1);
			}
			if(telmnt[0].tem.ag.RAP0flag==0)
			{
				fptname(fnminout.fnmnpp,dir_label);
				fnpp=fopen(fnminout.fnmnpp,"r");
				if(fnpp==NULL)
				{
					cout <<"Can not open file "<<fnminout.fnmnpp<<endl;
					exit(-1);
				}	
			}		
		}	
	}	
	if(tmflgs.kdinflg==1)
	{
		fptname(fnminout.fnmkdin,dir_label);
		fkdin=fopen(fnminout.fnmkdin,"r");
		if(fkdin==NULL)
		{
			cout <<"Can not open file "<<fnminout.fnmkdin<<endl;
			exit(-1);
		}	
	}

	if(tmflgs.kdoutflg==1)
	{
		fptname(fnminout.fnmkdout,dir_label);
		fkdout.open(fnminout.fnmkdout,ofstream::out);
		if(!fkdout.is_open())
		{
			cout <<"Can not create file "<<fnminout.fnmkdout<<endl;
			exit(-1);
		}	
	}	

	for (int kk = tmflgs.natmspred; kk< tmflgs.totpred; kk++)
	{
		fptname(fnminout.fnmtempred[kk],dir_label);
		ftempred[kk].open(fnminout.fnmtempred[kk],ofstream::out);
		if(!ftempred[kk].is_open())
		{
			cout <<"Failed to create file "<<fnminout.fnmtempred[kk]<<endl;
			exit(-1);
		}		
	}
	if(tmflgs.stateflag==1)
	{
		for(int kk=0;kk<13;kk++)
		{
			fptname(fnminout.fnmstate[kk],dir_label);
			fstate[kk]=fopen(fnminout.fnmstate[kk],"r");
			if(fstate[kk]==NULL)
			{
				cout <<"Can not open file "<<fnminout.fnmstate[kk]<<endl;
				exit(-1);
			}	
		}	
	}


  extrapolate();
//end of eidt

// Finished processing all elements - close open files

  if (fatalerr != 0)
  {
    if (elmnt.grdcnt != -99 && elmnt.count <= elmnt.grdcnt)
    {
      cout << "FATAL ERROR! Program Terminated" << endl;
    }
    flog1 << "FATAL ERROR! Program Terminated" << endl;
  }
  else
  {
    cout << "Extrapolation successfully completed - Congratulations!" << endl;
    flog1 << "Extrapolation successfully completed - Congratulations!" << endl;
  }

  if (tmflgs.atmsflag == 1)
  {
    fclose(fclds);
    if (telmnt[0].lonlatflag == 0) { fclose(flonlat); }
    if (tmflgs.atmsoutfg == 1)
    {
      for (i = 0; i < tmflgs.natmspred; i++) { fclmpred[i].close(); }
    }
  }
	if(tmflgs.kdinflg==1)fclose(fkdin);
	if(tmflgs.kdoutflg==1)fkdout.close();
 
  if (tmflgs.temflag == 1)
  {
    fclose(ftveg);
    fclose(fstxt);
    fclose(felev);
    fclose(fph);
    fclose(fdens);
    fclose(fpercn);
    if (tmflgs.atmsflag == 0)
    {
      fclose(fnirr);
      fclose(fpar);
    }
    fclose(ftair);
    fclose(fprec);
    fclose(fvap);

    if(telmnt[0].tem.ag.tlulcflag == 1)
    {
      fclose(flulc);
      if(telmnt[0].tem.ag.RAP0flag == 0) { fclose(fnpp); }
    }
    for (i = 0; i < tmflgs.ntempred; i++) { ftempred[i].close(); }
    if (tmflgs.stateflag == 1)
    {
      for (i = 0; i < MAXSTATE; i++) { fclose(fstate[i]); }
    }
  }

  flog1.close();
//YY added for mpi
  MPI::Finalize();
  return 0;
//end of edit
};

/* *************************************************************
*********************END OF MAIN PROGRAM************************
************************************************************* */


/* **************************************************************
************************************************************** */

int askpred(char pvarname[MAXPRED][9], int spred)
{
  int i;
  int j;
  int k;
  int t;
  int cnt;
  int numpred;
  int length;

  int posspred = MAXPRED;
  if (tmflgs.atmsoutfg == 1 && tmflgs.temflag == 0) { posspred = NUMATMS+1; }

  cout << endl << endl << "           POSSIBLE OUTPUT VARIABLES:" << endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
    for (t = 0; t < 5; t++)
    {
      cout << pvarname[(5*i)+t] << " ";
    }
    cout << endl;
  }

  flog1 << endl << endl << "           POSSIBLE OUTPUT VARIABLES:" << endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
     for (t = 0; t < 5; t++)
    {
      flog1 << pvarname[(5*i)+t] << ' ';
    }
    flog1 << endl;
  }

  cout << endl << endl << "How many variables are to be mapped (max " << posspred << ") in output files?  ";
  fpara >> numpred;
  cout << numpred << endl;

  flog1 << endl << endl << "How many variables are to be mapped (max " << posspred << ") in output files?" << numpred << endl << endl;

  cout << "Please enter output variable: " << endl;
  flog1 << "Please enter output variable: " << endl;

  for (i = 0; i < numpred; i++)
  {
    cnt = i + 1;
    k = i+spred;
    cout << cnt << " ";
    fpara >> tmflgs.predmap[k];
    length = strlen(tmflgs.predmap[k]);
    for (j = 0; j < length; j++) {tmflgs.predmap[k][j] = toupper(tmflgs.predmap[k][j]); }
    flog1 << cnt << " " << tmflgs.predmap[k] << endl;
  }

  return numpred;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void extrapolate(void)
{

  int i;
  int j;
  int k;
  int dm;
  int maxrec;
  int dyr;
  int itype;
//  double totamnt[NUMEQ+2];

  for (i = 1; i < MAXGRID; i++)
  {
    if (tmflgs.atmsflag == 1)
    {
      telmnt[i].clm.tcldsflag = telmnt[0].clm.tcldsflag;
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
    }
    if (tmflgs.temflag == 1)
    {
      telmnt[i].tem.atms.tco2flag = telmnt[0].tem.atms.tco2flag;
      telmnt[i].lonlatflag = telmnt[0].lonlatflag;
      telmnt[i].tem.atms.ttairflag = telmnt[0].tem.atms.ttairflag;
      telmnt[i].tem.atms.tprecflag = telmnt[0].tem.atms.tprecflag;
      telmnt[i].tem.ag.tlulcflag = telmnt[0].tem.ag.tlulcflag;

 // added for hydrology model
      telmnt[i].tem.hyd.vapflag = telmnt[0].tem.hyd.vapflag;

// Soil texture specific TEM parameters

      telmnt[i].tem.soil.pctpora = telmnt[0].tem.soil.pctpora;
      telmnt[i].tem.soil.pctporb = telmnt[0].tem.soil.pctporb;
      telmnt[i].tem.soil.fldcapa = telmnt[0].tem.soil.fldcapa;
      telmnt[i].tem.soil.fldcapb = telmnt[0].tem.soil.fldcapb;
      telmnt[i].tem.soil.wiltpta = telmnt[0].tem.soil.wiltpta;
      telmnt[i].tem.soil.wiltptb = telmnt[0].tem.soil.wiltptb;

// Initialize CO2 for all grid cells

      telmnt[i].tem.atms.initco2 = telmnt[0].tem.atms.initco2;
      telmnt[i].tem.atms.co2level = telmnt[0].tem.atms.co2level;
      if (tmflgs.equil == 0 && telmnt[0].tem.atms.tco2flag == 1)
      {
	for (j = 0; j < tmflgs.RTIME; j++)
        {
	  for (dm = 0; dm < CYCLE; dm++)
          {
 	    telmnt[i].tem.atms.tco2[j][dm] = telmnt[0].tem.atms.tco2[j][dm];
	  }
	}
      }
      for (j = 0; j < NUMVEG; j++)
      {
	telmnt[i].tem.soil.rootza[j] = telmnt[0].tem.soil.rootza[j];
	telmnt[i].tem.soil.rootzb[j] = telmnt[0].tem.soil.rootzb[j];
	telmnt[i].tem.soil.rootzc[j] = telmnt[0].tem.soil.rootzc[j];
	telmnt[i].tem.soil.minrootz[j] = telmnt[0].tem.soil.minrootz[j];
	telmnt[i].tem.veg.numtype[j] = telmnt[0].tem.veg.numtype[j];
	telmnt[i].tem.veg.minleaf[j] = telmnt[0].tem.veg.minleaf[j];
	telmnt[i].tem.veg.aleaf[j] = telmnt[0].tem.veg.aleaf[j];
	telmnt[i].tem.veg.bleaf[j] = telmnt[0].tem.veg.bleaf[j];
	telmnt[i].tem.veg.cleaf[j] = telmnt[0].tem.veg.cleaf[j];

//Site specific TEM parameters

	telmnt[i].tem.vegca[j] = telmnt[0].tem.vegca[j];
	telmnt[i].tem.vegcb[j] = telmnt[0].tem.vegcb[j];
	telmnt[i].tem.strna[j] = telmnt[0].tem.strna[j];
	telmnt[i].tem.strnb[j] = telmnt[0].tem.strnb[j];
	telmnt[i].tem.solca[j] = telmnt[0].tem.solca[j];
	telmnt[i].tem.solcb[j] = telmnt[0].tem.solcb[j];
	telmnt[i].tem.solna[j] = telmnt[0].tem.solna[j];
	telmnt[i].tem.solnb[j] = telmnt[0].tem.solnb[j];
	telmnt[i].tem.avlna[j] = telmnt[0].tem.avlna[j];
	telmnt[i].tem.avlnb[j] = telmnt[0].tem.avlnb[j];
	telmnt[i].tem.stona[j] = telmnt[0].tem.stona[j];
	telmnt[i].tem.stonb[j] = telmnt[0].tem.stonb[j];
	telmnt[i].tem.veg.unleaf12[j] = telmnt[0].tem.veg.unleaf12[j];
	telmnt[i].tem.veg.prvleafmx[j] = telmnt[0].tem.veg.prvleafmx[j];
	telmnt[i].tem.veg.cmaxcut[j] = telmnt[0].tem.veg.cmaxcut[j];
	telmnt[i].tem.veg.cmax1a[j] =  telmnt[0].tem.veg.cmax1a[j];
	telmnt[i].tem.veg.cmax1b[j] = telmnt[0].tem.veg.cmax1b[j];
	telmnt[i].tem.veg.cmax2a[j] = telmnt[0].tem.veg.cmax2a[j];
	telmnt[i].tem.veg.cmax2b[j] = telmnt[0].tem.veg.cmax2b[j];
	telmnt[i].tem.veg.cfall[j] = telmnt[0].tem.veg.cfall[j];
	telmnt[i].tem.veg.kra[j] = telmnt[0].tem.veg.kra[j];
	telmnt[i].tem.veg.krb[j] = telmnt[0].tem.veg.krb[j];
	telmnt[i].tem.microbe.kda[j] = telmnt[0].tem.microbe.kda[j];
	telmnt[i].tem.microbe.kdb[j] = telmnt[0].tem.microbe.kdb[j];
	telmnt[i].tem.microbe.lcclnc[j] = telmnt[0].tem.microbe.lcclnc[j];
	telmnt[i].tem.microbe.propftos[j] = telmnt[0].tem.microbe.propftos[j];
	telmnt[i].tem.veg.vsmmina[j] = telmnt[0].tem.veg.vsmmina[j];
	telmnt[i].tem.veg.vsmminb[j] = telmnt[0].tem.veg.vsmminb[j];
	telmnt[i].tem.veg.nmaxcut[j] = telmnt[0].tem.veg.nmaxcut[j];
	telmnt[i].tem.veg.nmax1a[j] = telmnt[0].tem.veg.nmax1a[j];
	telmnt[i].tem.veg.nmax1b[j] = telmnt[0].tem.veg.nmax1b[j];
	telmnt[i].tem.veg.nmax2a[j] = telmnt[0].tem.veg.nmax2a[j];
	telmnt[i].tem.veg.nmax2b[j] = telmnt[0].tem.veg.nmax2b[j];
	telmnt[i].tem.veg.nfall[j] = telmnt[0].tem.veg.nfall[j];
	telmnt[i].tem.microbe.nupa[j] = telmnt[0].tem.microbe.nupa[j];
	telmnt[i].tem.microbe.nupb[j] = telmnt[0].tem.microbe.nupb[j];
	telmnt[i].tem.soil.nloss[j] = telmnt[0].tem.soil.nloss[j];
	telmnt[i].tem.microbe.nfixpar[j] = telmnt[0].tem.microbe.nfixpar[j];
	telmnt[i].tem.veg.cneven[j] = telmnt[0].tem.veg.cneven[j];
	telmnt[i].tem.veg.cnmin[j] = telmnt[0].tem.veg.cnmin[j];
	telmnt[i].tem.veg.c2na[j] = telmnt[0].tem.veg.c2na[j];
	telmnt[i].tem.veg.c2nb[j] = telmnt[0].tem.veg.c2nb[j];
	telmnt[i].tem.veg.c2nmin[j] = telmnt[0].tem.veg.c2nmin[j];
	telmnt[i].tem.microbe.cnsoil[j] = telmnt[0].tem.microbe.cnsoil[j];

        // Vegetation specific TEM parameters

	telmnt[i].tem.veg.kc[j] = telmnt[0].tem.veg.kc[j];
	telmnt[i].tem.veg.ki[j] = telmnt[0].tem.veg.ki[j];
	telmnt[i].tem.veg.gva[j] = telmnt[0].tem.veg.gva[j];
	telmnt[i].tem.veg.tmin[j] = telmnt[0].tem.veg.tmin[j];
	telmnt[i].tem.veg.toptmin[j] = telmnt[0].tem.veg.toptmin[j];
	telmnt[i].tem.veg.toptmax[j] = telmnt[0].tem.veg.toptmax[j];
	telmnt[i].tem.veg.tmax[j] = telmnt[0].tem.veg.tmax[j];
	telmnt[i].tem.veg.raq10a0[j] = telmnt[0].tem.veg.raq10a0[j];
	telmnt[i].tem.veg.raq10a1[j] = telmnt[0].tem.veg.raq10a1[j];
	telmnt[i].tem.veg.raq10a2[j] = telmnt[0].tem.veg.raq10a2[j];
	telmnt[i].tem.veg.raq10a3[j] = telmnt[0].tem.veg.raq10a3[j];
	telmnt[i].tem.veg.kn1[j] = telmnt[0].tem.veg.kn1[j];
	telmnt[i].tem.veg.labncon[j] = telmnt[0].tem.veg.labncon[j];

        telmnt[i].tem.veg.leafmxc[j] = telmnt[0].tem.veg.leafmxc[j];
	telmnt[i].tem.veg.kleafc[j] = telmnt[0].tem.veg.kleafc[j];
	telmnt[i].tem.veg.sla[j] = telmnt[0].tem.veg.sla[j];
	telmnt[i].tem.veg.cov[j] = telmnt[0].tem.veg.cov[j];
	telmnt[i].tem.veg.fpcmax[j] = telmnt[0].tem.veg.fpcmax[j];

//Vegetation specific TEM parameters for microbes

	telmnt[i].tem.microbe.rhq10[j] = telmnt[0].tem.microbe.rhq10[j];
	telmnt[i].tem.microbe.kn2[j] = telmnt[0].tem.microbe.kn2[j];
	telmnt[i].tem.microbe.moistmin[j] = telmnt[0].tem.microbe.moistmin[j];
	telmnt[i].tem.microbe.moistopt[j] = telmnt[0].tem.microbe.moistopt[j];
	telmnt[i].tem.microbe.moistmax[j] = telmnt[0].tem.microbe.moistmax[j];
	telmnt[i].tem.microbe.k1[j] = telmnt[0].tem.microbe.k1[j];
	telmnt[i].tem.microbe.kmax[j] = telmnt[0].tem.microbe.kmax[j];
	telmnt[i].tem.microbe.fno3k[j] = telmnt[0].tem.microbe.fno3k[j];
	telmnt[i].tem.microbe.fco2k[j] = telmnt[0].tem.microbe.fco2k[j];
	telmnt[i].tem.microbe.kn[j] = telmnt[0].tem.microbe.kn[j];

        for (k=0; k < NUMMSAC; k++)
        {
	  telmnt[i].tem.veg.subtype[j][k] = telmnt[0].tem.veg.subtype[j][k];
	  telmnt[i].tem.veg.pcttype[j][k] = telmnt[0].tem.veg.pcttype[j][k];
	}

	telmnt[i].tem.ag.slashpar[j] = telmnt[0].tem.ag.slashpar[j];
	telmnt[i].tem.ag.vconvert[j] = telmnt[0].tem.ag.vconvert[j];
	telmnt[i].tem.ag.vrespar[j] = telmnt[0].tem.ag.vrespar[j];
	telmnt[i].tem.ag.sconvert[j] = telmnt[0].tem.ag.sconvert[j];
	telmnt[i].tem.ag.prod10par[j] = telmnt[0].tem.ag.prod10par[j];
	telmnt[i].tem.ag.prod100par[j] = telmnt[0].tem.ag.prod100par[j];
	telmnt[i].tem.ag.nvretconv[j] = telmnt[0].tem.ag.nvretconv[j];
	telmnt[i].tem.ag.nsretconv[j] = telmnt[0].tem.ag.nsretconv[j];
      }
    }
  }

  end1 = 1;
  fatalerr = 0;

//* *************************************************************

// Skip to the desired record in the GIS data sets

  if (elmnt.strtflag == 0)
  {
    for (i = 0; i < elmnt.numskip; i++)
    {
      if (tmflgs.atmsflag == 1)
      {
	end1 = telmnt[0].atmsgisin(fclds,tmflgs.cldflag,flonlat,telmnt[0].lonlatflag,
                                   tmflgs.numspin,tmflgs.spintime,tmflgs.RTIME);
	if (end1 == -1) { exit(0); }
      }

      if (tmflgs.temflag == 1)
      {
	end1 = telmnt[0].veggisin(flog1,fatalerr,tmflgs.atmsflag,tmflgs.temflag,ftveg);
	if (end1 == -1) { exit(0); }

        for (itype = 0; itype < telmnt[icount].maxtype; itype++)
        {
	  if (tmflgs.temflag == 1 && fatalerr == 0)
          {
	    end1 = telmnt[0].temgisin(flog1,fatalerr,tmflgs.atmsflag,
                                      itype,telmnt[icount].maxtype,
                                      tmflgs.kdinflg,tmflgs.stateflag,
	     	                      fstxt,felev, fph, fdens,fpercn,
                                      fnirr,fpar,
                                      ftair,fprec,fvap,
                                      flulc,fnpp,
                                      fkdin,fstate,
                                      tmflgs.numspin,tmflgs.spintime,tmflgs.RTIME);
            if (end1 == -1) { exit(0); }
	  }
	}
      }
    }
  }

//************************************************************** */

  cout << endl;
  flog1 << endl << endl;

  telmnt[0].col = MISSING;
  telmnt[0].row = MISSING;
  telmnt[0].tem.totyr = -99;

  elmnt.show(flog1, telmnt[0].col, telmnt[0].row,
             telmnt[0].tem.totyr, telmnt[0].tem.inittol);

  while (end1 != -1 && fatalerr == 0)
  {

    icount = 0;
    while (icount < MAXGRID && end1 != -1)
    {

// Get spatially-explicit cloudiness data   *********************

      if (tmflgs.atmsflag == 1)
      {
 	end1 = telmnt[icount].atmsgisin(fclds,tmflgs.cldflag,flonlat,
                                        telmnt[icount].lonlatflag,
                                        tmflgs.numspin,tmflgs.spintime,tmflgs.RTIME);
      }

// Get spatially-explicit vegetation data   *********************

      if (end1 != -1 && tmflgs.temflag == 1)
      {
	end1 = telmnt[icount].veggisin(flog1,fatalerr,tmflgs.atmsflag,tmflgs.temflag,ftveg);
      }

/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

      if (end1 != -1 && tmflgs.temflag == 1)
      {
        for (itype = 0; itype < telmnt[icount].maxtype; itype++)
        {
          // Get spatially-explicit input data for TEM

	  if (end1 != -1)
          {
	    end1 = telmnt[icount].temgisin(flog1,fatalerr,
                                           tmflgs.atmsflag,itype,
                                           telmnt[icount].maxtype,
                                           tmflgs.kdinflg,tmflgs.stateflag,
					   fstxt,felev, fph, fdens, fpercn,
                                           fnirr,fpar,
                                           ftair,fprec,fvap,
                                           flulc,fnpp,
                                           fkdin,fstate,
                                           tmflgs.numspin,tmflgs.spintime,tmflgs.RTIME);
          }
	}
      }

      if (end1 != -1) { ++icount; }
    }

    maxrec = icount;

    for (i = 0; i < maxrec; i++)
    {

      // Run TEM for "i" grid cells

      telmnt[i].runtem(flog1,tmflgs.predmap,
                       tmflgs.cldflag,tmflgs.atmsflag,tmflgs.atmsoutfg, tmflgs.natmspred, fclmpred,
                       tmflgs.temflag,tmflgs.kdinflg,tmflgs.ntempred,ftempred,
                       tmflgs.stateflag,tmflgs.equil,tmflgs.totsptime, tmflgs.RTIME,tmflgs.ispinout);

    }


/* *************************************************************
	  Write georeferenced data to output files
************************************************************* */

    for (i = 0; i < maxrec; i++)
    {
      for (dyr = 0; dyr < telmnt[i].outyr; dyr++)
      {
	     for (itype = 0; itype < telmnt[i].maxtype; itype++)
        {
	       if (tmflgs.temflag == 1)
          {
            // Write out spatially explicit kd parameters

	         if (dyr == 0 && tmflgs.kdoutflg == 1)
            {
	           kdparam.outdel(fkdout, telmnt[i].col, telmnt[i].row,
                             telmnt[i].tem.veg.temveg,
                             telmnt[i].tem.veg.subtype[telmnt[i].mez][itype],
                             telmnt[i].tem.microbe.kdsave[itype],
                             telmnt[i].contnent);
	         }
	       }
	     }
      }

      if (telmnt[i].tem.intflag > 0)
      {
        if (elmnt.count < elmnt.grdcnt)
        {
	       cout << "Integration terminated before attaining tolerance level" << endl;
	     }
	     flog1 << "Integration terminated before attaining tolerance level" << endl;
      }
      elmnt.show(flog1, telmnt[i].col, telmnt[i].row,
                 telmnt[i].tem.totyr, telmnt[i].tem.tol);
      if (elmnt.stopflag == 1) { end1 = -1; }
    }
  }

};

/* *************************************************************
************************************************************** */
//YY added for mpi
void subndtem(void)
{
	telmnt[0].tem.initrun(temstac);
	telmnt[0].tem.ask(temstac);
};
//end of edit

/* *************************************************************
************************************************************** */

void startclm(void)
{

  int i;
  char ifilename[25];
  char clmpredfile[25];

  tmflgs.cldflag = 1;
  cout << "Do you have spatially explicit solar radiation data or cloudiness data?:" << endl;
  cout << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
  cout << "Enter 1 for cloudiness data (percent cloudiness): ";
  fpara >> tmflgs.cldflag;

  flog1 << "Do you have spatially explicit solar radiation data or cloudiness data?:" << endl;
  flog1 << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
  flog1 << "Enter 1 for cloudiness data (percent cloudiness): " << endl;
  flog1 << "cldflag = " << tmflgs.cldflag << endl;

  telmnt[0].clm.tcldsflag = 0;
  if (tmflgs.cldflag == 1)
  {
    cout << "Do you have transient cloudiness data?:" << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    fpara >> telmnt[0].clm.tcldsflag;

    flog1 << endl << "Do you have transient cloudiness data?:" << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "telmnt[0].clm.tcldsflag = " << telmnt[0].clm.tcldsflag << endl << endl;

    cout << "Please enter the name of the file containing the mean monthly cloudiness data: " << endl;
    cout << "               (e.g., CLDS.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmclds;
    flog1 << "Please enter the name of the file containing the mean monthly cloudiness data: " << endl;
    flog1 << "               (e.g., CLDS.GIS) " << endl;
    flog1 << fnminout.fnmclds << std::endl << std::endl;
  }
  else
  {
    cout << "Do you have transient mean monthly solar radiation data?:" << endl;
    cout << "Enter 0 for no:" << endl;
    cout << "Enter 1 for yes: ";
    cin >> telmnt[0].clm.tcldsflag;

    flog1 << endl << "Do you have transient mean monthly solar radiation data?:" << endl;
    flog1 << "Enter 0 for no:" << endl;
    flog1 << "Enter 1 for yes: " << endl;
    flog1 << "telmnt[0].clm.tcldsflag = " << telmnt[0].clm.tcldsflag << endl << endl;

    cout << "Please enter the name of the file containing the mean monthly solar radiation data: " << endl;
    cout << "               (e.g., NIRR.GIS) " << endl;
    cin >> ifilename;
    flog1 << "Please enter the name of the file containing the mean monthly solar radiation data: " << endl;
    flog1 << "               (e.g., NIRR.GIS) " << endl;
    flog1 << ifilename << endl << endl;
  }

  //fclds = fopen(ifilename, "r");

  //if (!fclds)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

  cout << "How do you locate your grid cells?" << endl;
  cout << "Enter 0 for column/row:" << endl;
  cout << "Enter 1 for longitude/latitude: ";
  fpara >> telmnt[0].lonlatflag;

  flog1 << "How do you locate your grid cells?" << endl;
  flog1 << "Enter 0 for column/row:" << endl;
  flog1 << "Enter 1 for longitude/latitude: " << telmnt[0].lonlatflag << endl << endl;

  if (telmnt[0].lonlatflag == 0)
  {
    cout << "Please enter the name of the file containing the latitude data: " << endl;
    cout << "               (e.g., LAT.GIS) " << endl;
    cin >> ifilename;

    flog1 << "Please enter the name of the file containing the latitude data: " << endl;
    flog1 << "               (e.g., LAT.GIS) " << endl;
    flog1 << ifilename << endl << endl;

    //flonlat = fopen(ifilename, "r");

    //if (!flonlat)
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}
  }

  tmflgs.atmsoutfg = 0;
  cout << endl << "Do you wish output data from the irradiance model? " << endl;
  cout << "  Enter 0 for no" << endl;
  cout << "  Enter 1 for yes: ";
  fpara >> tmflgs.atmsoutfg;

  flog1 << endl << "Do you wish output data from the irradiance model? " << endl;
  flog1 << "  Enter 0 for no" << endl;
  flog1 << "  Enter 1 for yes: " << endl;
  flog1 << "atmsoutfg = " << tmflgs.atmsoutfg << endl << endl;

  if (tmflgs.atmsoutfg == 1)
  {
    tmflgs.totpred = tmflgs.natmspred = askpred(telmnt[0].clm.predstr, tmflgs.totpred);

    for (i = 0; i < tmflgs.natmspred; i++)
    {
      cout << endl << "Enter the name of the OUTPUT file to contain ";
      cout << tmflgs.predmap[i] << ":  ";
      //cin >> clmpredfile;
      fpara >> fnminout.fnmclmpred[i];
      flog1 << endl << "Enter the name of the OUTPUT file to contain ";
      flog1 << tmflgs.predmap[i] << ":  " << fnminout.fnmclmpred[i] << endl;
      //fclmpred[i].open(clmpredfile, ios::app);
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void starttem(void)
{

  int i;
  int numcmnt;
  char ifilename[40];
  char tempredfile[25];

  telmnt[0].tem.initrun(flog1,tmflgs.equil,temstac);
  telmnt[0].tem.ask(flog1,temstac);


// Get soil texture dependent parameters

  telmnt[0].tem.soil.getecd(flog1);


// Get vegetation type dependent parameters

  telmnt[0].tem.soil.getrootz(flog1);
  telmnt[0].tem.veg.getecd(flog1);
  telmnt[0].tem.veg.getleafecd(flog1);
  telmnt[0].tem.microbe.getvegecd(flog1);
// YYL added for penmon equation hydrology model, Q. Zhuang 
  telmnt[0].tem.hyd.getecd(flog1);
  telmnt[0].tem.hydm.getecdsoil(flog1);

//Get calibration site specific parameters

  cout << "Please enter the number of community types with calibration data:";
  //cin >> numcmnt;
  fpara >> numcmnt;
  flog1 << endl << endl << "Please enter the number of community types with calibration data:";
  flog1 << numcmnt << endl;

  telmnt[0].tem.getsitecd(numcmnt, flog1);


// Get vegetation mosaic information

  telmnt[0].tem.veg.getvtype(flog1);


  cout << "Please enter the name of the file containing the soil texture data:";
  cout << endl;
  cout << "               (e.g., TEXTURE.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmstxt;
  
  flog1 << "Please enter the name of the file containing the soil texture data:";
  flog1 << endl;
  flog1 << "               (e.g., TEXTURE.GIS) " << endl;
  flog1 << fnminout.fnmstxt << endl << endl;
  //fstxt = fopen(ifilename, "r");

  //if (!fstxt)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

  cout << "Please enter the name of the file containing the vegetation data: " << endl;
  cout << "               (e.g., TEMVEG.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmtveg; 
  
  flog1 << "Please enter the name of the file containing the vegetation data: " << endl;
  flog1 << "               (e.g., TEMVEG.GIS) " << endl;
  flog1 << fnminout.fnmtveg << endl << endl;
  //ftveg = fopen(ifilename, "r");

  //if (!ftveg)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

  cout << "Please enter the name of the file containing the elevation data: " << endl;
  cout << "               (e.g., ELEV.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmelev;
  flog1 << "Please enter the name of the file containing the elevation data: " << endl;
  flog1 << "               (e.g., ELEV.GIS) " << endl;
  flog1 << fnminout.fnmelev << endl << endl;
  //felev = fopen(ifilename, "r");

  //if (!felev)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

// added for reading pH value
  cout << "Please enter the name of the file containing the pH data: " << endl;
  cout << "               (e.g., pH.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmph;  
  
  flog1 << "Please enter the name of the file containing the pH data: " << endl;
  flog1 << "               (e.g., pH.GIS) " << endl;
  flog1 << fnminout.fnmph << endl << endl;
  //fph = fopen(ifilename, "r");

  //if (!fph)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}
  
  // added for reading Density value
  cout << "Please enter the name of the file containing the Density data: " << endl;
  cout << "               (e.g.,  Density.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmdens;
  
  flog1 << "Please enter the name of the file containing the  Density data: " << endl;
  flog1 << "               (e.g.,  Density.GIS) " << endl;
  flog1 << fnminout.fnmdens << endl << endl;
  //fdens = fopen(ifilename, "r");

  //if (!fdens)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}
  // added for reading permafrost CN value
  cout << "Please enter the name of the file containing the permafrost CN data: " << endl;
  cout << "               (e.g.,  percn.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmpercn;
  
  flog1 << "Please enter the name of the file containing the permafrost CN data: " << endl;
  flog1 << "               (e.g.,  percn.GIS) " << endl;
  flog1 << fnminout.fnmpercn << endl << endl;

  telmnt[0].tem.atms.ttairflag = 0;
  if (tmflgs.equil == 0)
  {
    cout << "Do you have transient air temperature data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.ttairflag;

    flog1 << endl << "Do you have transient air temperature data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.atms.ttairflag = " << telmnt[0].tem.atms.ttairflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the mean monthly air temperature data: " << endl;
  cout << "               (e.g., TAIR.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmtair;
  flog1 << "Please enter the name of the file containing the mean monthly air temperature data: " << endl;
  flog1 << "               (e.g., TAIR.GIS) " << endl;
  flog1 << fnminout.fnmtair << endl << endl;
  //ftair = fopen(ifilename, "r");

  //if (!ftair)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

  telmnt[0].tem.atms.tprecflag = 0;
  if (tmflgs.equil == 0)
  {
    cout << "Do you have transient precipitation data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.tprecflag;

    flog1 << endl << "Do you have transient precipitation data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "temnt[0].tem.atms.tprecflag = " << telmnt[0].tem.atms.tprecflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the monthly precipitation data: " << endl;
  cout << "               (e.g., PREC.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmprec;
  
  flog1 << "Please enter the name of the file containing the monthly precipitation data: " << endl;
  flog1 << "               (e.g., PREC.GIS) " << endl;
  flog1 << fnminout.fnmprec << endl << endl;
  //fprec = fopen(ifilename, "r");

  //if (!fprec)
  //{
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

 // added for hydrology model, Q. Zhuang

  telmnt[0].tem.hyd.vapflag = 0;
  if (tmflgs.equil == 0) {
    cout << "Do you have transient vapor pressure data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
  	fpara >> telmnt[0].tem.hyd.vapflag;

    flog1 << endl << "Do you have transient vapor pressure data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.hyd.vapflag = " << telmnt[0].tem.hyd.vapflag << endl << endl;
  }

  cout << "Please enter the name of the file containing the mean monthly vapor pressure data: " << endl;
  cout << "               (e.g., tvap.GIS) " << endl;
  //cin >> ifilename;
  fpara >> fnminout.fnmvap;  
  flog1 << "Please enter the name of the file containing the mean monthly vapor pressure data: " << endl;
  flog1 << "               (e.g., TVAP.GIS) " << endl;
  flog1 << fnminout.fnmvap << endl << endl;
  //fvap = fopen(ifilename, "r");

  //if (!fvap) {
    //cerr << "\nCannot open " << ifilename << " for data input" << endl;
    //exit(-1);
  //}

  // end of adding for hydrology model
  if (tmflgs.atmsflag == 0)
  {
    cout << "Please enter the name of the file containing the net irradiance data: " << endl;
    cout << "               (e.g., NIRR.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmnirr;
    flog1 << "Please enter the name of the file containing the net irradiance data: " << endl;
    flog1 << "               (e.g., NIRR.GIS) " << endl;
    flog1 <<  fnminout.fnmnirr << endl << endl;
    //fnirr = fopen(ifilename, "r");

    //if (!fnirr)
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the mean monthly photosynthetically active radiation data: " << endl;
    cout << "               (e.g., PAR.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmpar;
    flog1 << "Please enter the name of the file containing the mean monthly photosynthetically active radiation data: " << endl;
    flog1 << "               (e.g., PAR.GIS) " << endl;
    flog1 << fnminout.fnmpar << endl << endl;
    //fpar = fopen(ifilename, "r");

    //if (!fpar)
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}
  }

  cout << endl << endl << "Enter the initial concentration of carbon dioxide in ppmv: ";
  fpara >> telmnt[0].tem.atms.initco2;
  flog1 << endl << endl << "Enter the initial concentration of carbon dioxide in ppmv: " << telmnt[0].tem.atms.initco2 << endl << endl;

  cout << endl << endl << "Enter the final equilibrium concentration of carbon dioxide in ppmv: ";
  fpara >> telmnt[0].tem.atms.co2level;
  flog1 << endl << endl << "Enter the final equilibrium concentration of carbon dioxide in ppmv: " << telmnt[0].tem.atms.co2level << endl << endl;

  telmnt[0].tem.atms.tco2flag = 0;
  if (tmflgs.equil == 0)
  {
    cout << "Do you have transient CO2 data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.atms.tco2flag;

    flog1 << "Do you have transient CO2 data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: " << endl;
    flog1 << "telmnt[0].tem.atms.tco2flag = " << telmnt[0].tem.atms.tco2flag << endl << endl;

    if (telmnt[0].tem.atms.tco2flag == 1)
    {
      telmnt[0].tem.atms.loadmpitCO2(flog1, tmflgs.totsptime, tmflgs.RTIME);
    }
  }

  cout << endl << endl << "Enter the factor for changing C:N per ppmv of enhanced CO2:" << endl;
  cout << "                     (Enter 0.0 for no change): " << endl;
  fpara >> telmnt[0].tem.veg.dc2n;

  flog1 << endl << "Enter the factor for changing C:N per ppmv of enhanced CO2:" << endl;
  flog1 << "                     (Enter 0.0 for no change): " << endl;
  flog1 << "telmnt[0].tem.veg.dc2n = " << telmnt[0].tem.veg.dc2n << endl << endl;

  telmnt[0].tem.ag.tlulcflag = 0;
  if (tmflgs.equil == 0)
  {
    cout << "Do you have transient land use data?:" << endl;
    cout << "Enter 0 for No:" << endl;
    cout << "Enter 1 for Yes: ";
    fpara >> telmnt[0].tem.ag.tlulcflag;

    flog1 << "Do you have transient land use data?:" << endl;
    flog1 << "Enter 0 for No:" << endl;
    flog1 << "Enter 1 for Yes: ";
    flog1 << "telmnt[0].tem.ag.tlulcflag = " << telmnt[0].tem.ag.tlulcflag << endl << endl;;

    if (telmnt[0].tem.ag.tlulcflag == 1)
    {
      cout << "Please enter the name of the file containing the land use data: " << endl;
      cout << "               (e.g., LULC.GIS) " << endl;
      //cin >> ifilename;
      fpara >> fnminout.fnmlulc;
      flog1 << "Please enter the name of the file containing the land use data: " << endl;
      flog1 << "               (e.g., LULC.GIS) " << endl;
      flog1 << fnminout.fnmlulc << endl << endl;
      //flulc = fopen(ifilename, "r");

      //if (!flulc)
      //{
	//cerr << "\nCannot open " << ifilename << " for data input" << endl;
	//exit(-1);
      //}

      cout << "Will Relative Agricultural Production (RAP) always be 0.0?" << endl;
      cout << "Enter 0 for No:" << endl;
      cout << "Enter 1 for Yes: ";
      fpara >> telmnt[0].tem.ag.RAP0flag;

      flog1 << "Will Relative Agricultural Production (RAP) always be 0.0?" << endl;
      flog1 << "Enter 0 for No:" << endl;
      flog1 << "Enter 1 for Yes: " << endl;
      flog1 << "telmnt[0].tem.ag.RAP0flag = " << telmnt[0].tem.ag.RAP0flag << endl;

      if(telmnt[0].tem.ag.RAP0flag == 0)
      {
 	cout << "Please enter the name of the file containing the monthly potential NPP data: " << endl;
 	cout << "               (e.g., NPP.GIS) " << endl;
 	//cin >> ifilename;
  fpara >> fnminout.fnmnpp;
 	flog1 << "Please enter the name of the file containing the monthly potential NPP data: " << endl;
 	flog1 << "               (e.g., NPP.GIS) " << endl;
	flog1 << fnminout.fnmnpp << endl << endl;
	//fnpp = fopen(ifilename, "r");

	//if (!fnpp)
        //{
	  //cerr << "\nCannot open " << ifilename << " for data input" << endl;
	  //exit(-1);
	//}
      }

      telmnt[0].tem.ag.getecd(flog1);
    }
  }

  tmflgs.kdinflg = 0;
  cout << "Do you want to use spatially explicit values of Kd?" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
  fpara >> tmflgs.kdinflg;

  flog1 << "Do you want to use spatially explicit values of Kd?" << endl;
  flog1 << "Enter 0 for No:" << endl;
  flog1 << "Enter 1 for Yes: " << endl;
  flog1 << "kdinflg = " << tmflgs.kdinflg << endl << endl;

  if (tmflgs.kdinflg == 1)
  {
    cout << "Please enter the name of the file containing the Kd values: " << endl;
    cout << "               (e.g., KDIN.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmkdin;
    flog1 << "Please enter the name of the file containing the Kd values: " << endl;
    flog1 << "               (e.g., KDIN.GIS) " << endl;
    flog1 << fnminout.fnmkdin << endl << endl;

    //fkdin = fopen(ifilename, "r");

    //if (!fkdin)
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}
  }

  tmflgs.kdoutflg = 0;
  cout << "Do you want to generate spatially explicit values of Kd?" << endl;
  cout << "Enter 0 for No:" << endl;
  cout << "Enter 1 for Yes: ";
  fpara >> tmflgs.kdoutflg;

  flog1 << endl << "Do you want to generate spatially explicit values of Kd?" << endl;
  flog1 << "Enter 0 for No:" << endl;
  flog1 << "Enter 1 for Yes: " << endl;
  flog1 << "kdoutflg = " << tmflgs.kdoutflg << endl << endl;

  if (tmflgs.kdoutflg == 1)
  {
    cout << "Please enter the name of the file to contain the Kd values: " << endl;
    cout << "               (e.g., KDOUT.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmkdout;
    flog1 << "Please enter the name of the file to contain the Kd values: " << endl;
    flog1 << "               (e.g., KDOUT.GIS) " << endl;
    flog1 << fnminout.fnmkdout << endl << endl;

    //fkdout.open(ifilename, ios::app);
  }

  tmflgs.ntempred = askpred(telmnt[0].tem.predstr,tmflgs.totpred);
  tmflgs.totpred += tmflgs.ntempred;

  for (i = tmflgs.natmspred; i < tmflgs.totpred; i++)
  {
    cout << endl << "Enter the name of the OUTPUT file to contain " << tmflgs.predmap[i] << ":  ";
    //cin >> tempredfile;
    fpara >> fnminout.fnmtempred[i];
    flog1 << "Enter the name of the OUTPUT file to contain " << tmflgs.predmap[i] << ":  " << fnminout.fnmtempred[i] << endl;
    //ftempred[i].open(tempredfile, ios::app);
    otemls << fnminout.fnmtempred[i] << std::endl;
  }

  tmflgs.stateflag = 0;
  cout << endl << "Do you want to use spatially explicit data for initial conditions? " << endl;
  cout << "  Enter 1 for YES:" << endl;
  cout << "  Enter 0 for NO:" << endl;
  cout << "  NOTE:  If YES, you will need spatially explicit data for VEGC," << endl;
  cout << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE, AVAILH2O, RGRNDH2O," << endl;
  cout << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET, EET,(i.e. 13 data files)" << endl;
  fpara >> tmflgs.stateflag;

  flog1 << endl << "Do you want to use spatially explicit data for intial conditions? " << endl;
  flog1 << "  Enter 1 for YES:" << endl;
  flog1 << "  Enter 0 for NO:" << endl;
  flog1 << "  NOTE:  If YES, you will need spatially explicit data for VEGC," << endl;
  flog1 << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE, AVAILH2O, RGRNDH2O," << endl;
  flog1 << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET, EET,(i.e. 13 data files)" << endl;
  flog1 << "stateflag = " << tmflgs.stateflag << endl << endl;

  if (tmflgs.stateflag == 1)
  {
    cout << "Please enter the name of the file containing the vegetation carbon data: " << endl;
    cout << "               (e.g., VEGC.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[0];
    flog1 << "Please enter the name of the file containing the vegetation carbon data: " << endl;
    flog1 << "               (e.g., VEGC.GIS) " << endl;
    flog1 << fnminout.fnmstate[0] << endl << endl;
    //fstate[0] = fopen(ifilename, "r");

    //if (!fstate[0])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the vegetation structural nitrogen data: " << endl;
    cout << "               (e.g., STRN.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[1];
    flog1 << "Please enter the name of the file containing the vegetation structural nitrogen data: " << endl;
    flog1 << "               (e.g., STRN.GIS) " << endl;
    flog1 << fnminout.fnmstate[1] << endl << endl;
    //fstate[1] = fopen(ifilename, "r");

    //if (!fstate[1])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the soil organic carbon data: " << endl;
    cout << "               (e.g., SOLC.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[2];
    flog1 << "Please enter the name of the file containing the soil organic carbon data: " << endl;
    flog1 << "               (e.g., SOLC.GIS) " << endl;
    flog1 << fnminout.fnmstate[2] << endl << endl;
    //fstate[2] = fopen(ifilename, "r");

    //if (!fstate[2])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the soil organic nitrogen data: " << endl;
    cout << "               (e.g., SOLN.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[3];
    flog1 << "Please enter the name of the file containing the soil organic nitrogen data: " << endl;
    flog1 << "               (e.g., SOLN.GIS) " << endl;
    flog1 << fnminout.fnmstate[3] << endl << endl;
    //fstate[3] = fopen(ifilename, "r");

    //if (!fstate[3])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the soil available nitrogen data: " << endl;
    cout << "               (e.g., AVLN.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[4];
    flog1 << "Please enter the name of the file containing the soil availablenitrogen data: " << endl;
    flog1 << "               (e.g., AVLN.GIS) " << endl;
    flog1 << fnminout.fnmstate[4] << endl << endl;
    //fstate[4] = fopen(ifilename, "r");

    //if (!fstate[4])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the vegetation labile nitrogen data: " << endl;
    cout << "               (e.g., STON.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[5];
    flog1 << "Please enter the name of the file containing the vegetation labile nitrogen data: " << endl;
    flog1 << "               (e.g., STON.GIS) " << endl;
    flog1 << fnminout.fnmstate[5] << endl << endl;
    //fstate[5] = fopen(ifilename, "r");

    //if (!fstate[5])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the available soil moisture data: " << endl;
    cout << "               (e.g., AVLW.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[6];
    flog1 << "Please enter the name of the file containing the available soil moisture data: " << endl;
    flog1 << "               (e.g., AVLW.GIS) " << endl;
    flog1 << fnminout.fnmstate[6]  << endl << endl;
    //fstate[6] = fopen(ifilename, "r");

    //if (!fstate[6])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the rain-derived groundwater data: " << endl;
    cout << "               (e.g., RGRW.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[7];
    flog1 << "Please enter the name of the file containing the rain-derived groundwater data: " << endl;
    flog1 << "               (e.g., RGRW.GIS) " << endl;
    flog1 << fnminout.fnmstate[7]  << endl << endl;
    //fstate[7] = fopen(ifilename, "r");

    //if (!fstate[7])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the snowpack data: " << endl;
    cout << "               (e.g., SNWP.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[8];
    flog1 << "Please enter the name of the file containing the snowpack data: " << endl;
    flog1 << "               (e.g., SNWP.GIS) " << endl;
    flog1 << fnminout.fnmstate[8] << endl << endl;
    //fstate[8] = fopen(ifilename, "r");

    //if (!fstate[8])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the snowpack-derived groundwater data: " << endl;
    cout << "               (e.g., SGRW.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[9];
    flog1 << "Please enter the name of the file containing the snowpack-derived groundwater data: " << endl;
    flog1 << "               (e.g., SGRW.GIS) " << endl;
    flog1 << fnminout.fnmstate[9] << endl << endl;
    //fstate[9] = fopen(ifilename, "r");

    //if (!fstate[9])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the unnormalized leaf phenology data: " << endl;
    cout << "               (e.g., UNLF.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[10];
    flog1 << "Please enter the name of the file containing the unnormalized leaf phenology data: " << endl;
    flog1 << "               (e.g., UNLF.GIS) " << endl;
    flog1 << fnminout.fnmstate[10] << endl << endl;
    //fstate[10] = fopen(ifilename, "r");

    //if (!fstate[10])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the potential evapotranspiration data: " << endl;
    cout << "               (e.g., PET.GIS) " << endl;
    //cin >> ifilename;
    fpara >> fnminout.fnmstate[11];
    flog1 << "Please enter the name of the file containing the potential evapotranspiration data: " << endl;
    flog1 << "               (e.g., PET.GIS) " << endl;
    flog1 << fnminout.fnmstate[11] << endl << endl;
    //fstate[11] = fopen(ifilename, "r");

    //if (!fstate[11])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}

    cout << "Please enter the name of the file containing the estimated evapotranspiration data: " << endl;
    cout << "               (e.g., EET.GIS) " << endl;
  	fpara >> fnminout.fnmstate[12];
    flog1 << "Please enter the name of the file containing the estimated evapotranspiration data: " << endl;
    flog1 << "               (e.g., EET.GIS) " << endl;
    flog1 << fnminout.fnmstate[12] << endl << endl;
    //fstate[12] = fopen(ifilename, "r");

    //if (!fstate[12])
    //{
      //cerr << "\nCannot open " << ifilename << " for data input" << endl;
      //exit(-1);
    //}
  }

};

//YY added for mpi
void fptname(char* fnmtemfl,int rank)
{
	char chtmp[220];
	strcpy(chtmp, fnmtemfl);	
	sprintf(fnmtemfl,"./part-%d/%s-%d", rank, chtmp, rank);
	return ;
};
