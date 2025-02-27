#if !defined(QSOILTEMP_H)
  #include "qsoiltemp.hpp"
#endif

using namespace std;
// Soil thermal model soiltemp.cpp

ifstream snowdepth, initfile;
Soilthermal basicf;  // initializing parameters
char dummy[60];

struct {
    double dt, dtday, per, first, time, hlat, tf, sigma, emiss, htop, hbot,
	    tair, sflux, tbot, gflux, smass, dtfaz, final, theta, theta1,
	    spare[6], x[211], dx[211], xfa[211], xfb[211], xhet[211], water[
	    211], ddry[211], condt[30], condf[30], spht[30], sphf[30], conx[
	    211], capx[211], e[211], s[211], t[211], told[211], ht[211],
	    htold[211], eltim[150000], caltim[150000],arr[3000];   //zbl add arr[3000] to prevent infinite loop;
    integer ism1, is, igm1, ig, imax1, imax, mmax, mater[211], nan, nst;
}  basicdata;
#define basicdata1 basicdata

// soil water
struct {
    double soilwater1, soilwater2, soilwater3;
} blkwater_;
#define blkw_ blkwater_

// snow classification
struct {
    double cdsnow;
} blk1_;
#define blk1_1 blk1_

// air temperature
struct {
    double airt1, airt2, airt3;
} blk2_;
#define blk2_1 blk2_

// snow height
struct {
    double hsnow1, hsnow2, hsnow3;
} blk3_;
#define blk3_1 blk3_

int Soilthermal::soiltemp_(double *water2, double *calcu_cdsnow, double *airt19, double *airt29, double *airt39, double
	*hsnow19, double *hsnow29, double *hsnow39, double *weight9, double *smass9, integer *is9, integer *ism19, double *x9, 
  double *dx9, double *xfb9, double *xfa9, double *water9, double *t9, double *tsoil, double *frontd, double *thawbegin, double *thawend, double diffsoilt[10], int& cmnt,	double& tmax, double& tmin, double& tave, double& longmean) //YY 072822 ,	double& tmax, double& tmin, double& tave, double& longmean,  double& vsm1, double& vsm2, double& vsm3, double& pthick
{
    static integer i__1= 0; //YY 072822 add=0 
    static integer first_phase=0,second_phase=0; 
    double r__1 =0;  
    static double adif[211] = {0.0}, dep10 = 0.0, dep11 = 0.0, dep12 = 0.0, dep13 = 0.0, dep14 = 0.0, dep15 = 0.0, dep16 = 0.0, dep17 = 0.0, heat[6000] = {0.0}, dep18 = 0.0;  
    static integer kall = 0;
    static double acum = 0.0, tdif = 0.0, sden[3000] = {0.0}, taer[3000] = {0.0}, tber[3000] = {0.0};
    static integer kint[2] = {0}, maxm = 0, kiso = 0;
    static double tiso[10] = {0.0};
    static integer liso = 0, lmax = 0;
    static double vdep = 0.0, turb = 0.0;
    static integer kfor = 0, ktop = 0 ;
	  static double snow[3000] = {0.0}, comp = 0.0, frza = 0.0, thwa = 0.0, tslo = 0.0, tshi = 0.0, thws = 0.0, frzs = 0.0, thwg = 0.0, frzg = 0.0;
    static integer nsum = 0;
    extern  int surf_(integer *, double *, double *, double *, double*), pram_(integer *, double *, double *, double *, double *, double *, double*), tblo_(integer *);
    static double tset = 0, sflx1 = 0, sflx2 = 0, sflx3 = 0, g = 0;
    static integer j=0, l=0, m=0;
    extern  int neige_(double *, double *, double *, double *, double *,double *, double *, double *);
    static integer nread = 0;
    static double tanee[211] = {0.0};
    static integer ldepf = 0;
    extern  int bcbot_(integer *, double *);
    static integer iheat = 0;
    static double thalt[211] = {0.0}, heatt[3000] = {0.0};
    static integer ktemp = 0, kenvl = 0, kflux = 0, imaxp = 0;
    extern  int bctop_(integer *, double *, integer *, double *,double *, double *);
	  static double tgrlo = 0.0, tgrhi = 0.0, tgran = 0.0;
    static integer ksnow = 0;
	  static double tsurf = 0.0, tgrnd = 0.0, denfac = 0.0;
    static integer jj = 0;
	  static double calend = 0.0, td = 0.0, tt = 0.0, fcmelt = 0.0, tu = 0.0, adifmx = 0.0;
    extern  int asmblo_(integer *);
	  static double deptem[20] = {0.0}, tairan = 0.0, tairhi = 0.0, depflx[20] = {0.0};
    static integer nanmax = 0, knodes = 0;
    extern  int snofal_(double *, double *, double *, double *, double * , double *);

    static double snoden = 0.0, weight[10] = {0.0};
    extern  int asmone_(integer *);
    static double tairlo = 0.0;
    extern  int inthet_(integer *, double *), sethet_(integer *, integer *, double *, double *);

    static integer kallyr = 0;
    static double tdrive = 0.0;
    static integer ibzero = 0;
    static double topold = 0.0, thalto = 0.0;
    static integer nwrite = 0, itzero = 0;
    static double topnew = 0.0;
    extern  int asmtwo_(integer *, integer *);
    static double dif = 0.0, cnd = 0.0, rad = 0.0;
    static integer jjj, kkk;
   	static double thi[211] = {0.0}, sph = 0.0, tlo[211] = {0.0}, dxx = 0.0, dep1 = 0.0, dep2 = 0.0, dep3 = 0.0, dep4 = 0.0, dep5 = 0.0, dep6 = 0.0, dep7 = 0.0, dep8 = 0.0, dep9 = 0.0, eta0 = 0.0;
   
    static integer isp1= 0,itz1= 0;
    double phase_count=0.0;
	  double maxtemp = 0.0;
	  double mintemp = 0.0;	
	  double avetemp = 0.0;
	  static double newtemp[25] = {0.0};
     
   	int debug1=0;
    --t9;
    --water9;
    --xfa9;
    --xfb9;
    --dx9;
    --x9;
    --weight9;

//	printf(" 2# %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n ", t9, water9, xfa9,dx9,x9,weight9);
//printf(" %d", cmnt);

   blk1_1.cdsnow= *calcu_cdsnow; // passing from the snow classification Q. Z. 28/mar/2001

    blk2_1.airt1 = *airt19;
    blk2_1.airt2 = *airt29;
    blk2_1.airt3 = *airt39;
    blk3_1.hsnow1 = *hsnow19;
    blk3_1.hsnow2 = *hsnow29;
    blk3_1.hsnow3 = *hsnow39;


  // linkage of STM with HYM

    blkw_.soilwater1 = *water2;
    blkw_.soilwater2 = *water2;
//    blkw_.soilwater3 = *water3;

    for (i__ = 1; i__ <= 210; ++i__) {
	basicdata1.x[i__ - 1] = x9[i__];
	basicdata1.dx[i__ - 1] = dx9[i__];
	basicdata1.xfa[i__ - 1] = xfa9[i__];
	basicdata1.xfb[i__ - 1] = xfb9[i__];
	basicdata1.water[i__ - 1] = water9[i__];
	basicdata1.t[i__ - 1] = t9[i__];
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	weight[i__ - 1] = weight9[i__];
    }
    maxm = 3000;
        	//zbl add arr[3000] to prevent infinite loop
    	for (i__ = 1; i__<= 3000; ++i__) {
		  basicdata1.arr[i__ - 1] = -99.0;
    	}  //end of add
    jjj = 98;
//    if (kswitch == 0) {

//    initfile.open("initl.dat",ios::in);
//    snowdepth.open("snoda.dat",ios::in);
//    snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
//    snowdepth >>nanmax >> basicdata1.nst >> kallyr >> knodes >> kiso >> ktemp >> ksnow >> kenvl >> kflux;
    nanmax=basicf.MAX[cmnt]; 
    basicdata1.nst= basicf.NST[cmnt]; 
    kallyr= basicf.KALLYR[cmnt];
    knodes = basicf.KNODES[cmnt]; 
    kiso = basicf.KISO[cmnt];
    ktemp= basicf.KTEMP[cmnt]; 
    ksnow= basicf.KSNOW[cmnt];
    kenvl=basicf.KENVL[cmnt]; 
    kflux=basicf.KFLUX[cmnt];
//    printf("*** %d",nanmax);
//    getch();

    mynanmax = nanmax;
    mykallyr = kallyr;
    myknodes =knodes;
    mykiso =kiso;
    myktemp=ktemp;
    myksnow = ksnow;
    mykenvl=kenvl;
    mykflux =kflux;
    my_nst=basicdata1.nst;
//     }
/*   else {
    nanmax = mynanmax;
    kallyr = mykallyr;
    knodes =myknodes;
    kiso =mykiso;
    ktemp=myktemp;
    ksnow = myksnow;
    kenvl=mykenvl;
    kflux =mykflux;
    basicdata1.nst=my_nst;

   }
 */
     kall = kallyr;
 //  if (kswitch ==0) {
//    snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
//    snowdepth >> liso;
    liso= basicf.LISO[cmnt];

    myliso = liso;
//   }
 //  else {   liso = myliso;
 //  }

    i__1 = liso;
//    if (kswitch ==0 ) {
         for (l = 1; l <= i__1; ++l) {
       //  snowdepth >> tiso[l - 1];
           tiso[l - 1] = basicf.TISO[cmnt];

         mytiso[l-1]=tiso[l-1]; }
  //   }
 //    else {
 //         for (l = 1; l <= i__1; ++l) {
//           tiso[l-1]=mytiso[l-1]; }
//          }
//    if (kswitch ==0) {
      //  snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
     //   snowdepth >> lmax >> vdep;
      lmax= basicf.LMAX[cmnt];
      vdep = basicf.VDEPTH[cmnt];

       // cout << endl <<lmax << " " << vdep;
       mylmax = lmax; myvdep=vdep;
  //    }
 //   else { lmax = mylmax; vdep=myvdep; }

    i__1 = lmax;
//    if (kswitch == 0) {
        for (l = 1; l <= i__1; ++l) {
      //   snowdepth >>deptem[l - 1];
       deptem[0]= basicf.DEPTEM1[cmnt];
       deptem[1]= basicf.DEPTEM2[cmnt];
       deptem[2]= basicf.DEPTEM3[cmnt];
       deptem[3]= basicf.DEPTEM4[cmnt];
       deptem[4]= basicf.DEPTEM5[cmnt];

       mydeptem[l-1] = deptem[l-1];
        }
//     }
//     else {
//      for (l = 1; l <= i__1; ++l)  deptem[l-1] = mydeptem[l-1];
//      }

    i__1 = lmax;
    for (l = 1; l <= i__1; ++l) {
   	deptem[l - 1] = vdep * deptem[l - 1];
    }

//   if (kswitch ==0 ) {
//   snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
//   snowdepth >> ldepf >> vdep;
      ldepf = basicf.NDEPF[cmnt];
      vdep= basicf.VDEP[cmnt];

   // cout << endl<<ldepf<<" " << vdep <<endl;
   myldepf= ldepf; myvdep1= vdep;
//     }
//    else {ldepf= myldepf; vdep= myvdep1; }

   i__1 = ldepf;

//   if (kswitch ==0 ) {
      for (l = 1; l <= i__1; ++l) {
     //	snowdepth >>depflx[l - 1];
       depflx[0]=basicf.DEPFLX1[cmnt];
       depflx[1]=basicf.DEPFLX2[cmnt];
       depflx[2]=basicf.DEPFLX3[cmnt];
       depflx[3]=basicf.DEPFLX4[cmnt];
       depflx[4]=basicf.DEPFLX5[cmnt];

       mydepflx[l - 1]=depflx[l - 1];
      }
//    }
//    else {  for (l = 1; l <= i__1; ++l)	depflx[l - 1]=mydepflx[l - 1];
//         }

    i__1 = ldepf;
    for (l = 1; l <= i__1; ++l) {
	depflx[l - 1] = vdep * depflx[l - 1];
    }

    basicdata1.time = 0.0;
    basicdata1.sigma = 0.0;
    rad = 0.0;
    turb = 0.0;

//   if (kswitch ==0 ) {
//      snowdepth >> dummy >> dummy >> dummy >> dummy;
//      snowdepth >>basicdata1.hlat >>basicdata1.tf >>basicdata1.gflux >>blk1_1.cdsnow;
      basicdata1.hlat= basicf.HLAT[cmnt];
      basicdata1.tf= basicf.TF[cmnt];
      basicdata1.gflux=basicf.gflux[cmnt];
//      blk1_1.cdsnow = basicf.cdsnow[cmnt];

      my_hlat=basicdata1.hlat;
      my_tf=basicdata1.tf;
      my_gflux=basicdata1.gflux;
//      my1_cdsnow=blk1_1.cdsnow;
      // cout <<endl<<blk1_1.cdsnow<<endl;
//    }
//    else {
//       basicdata1.hlat=my_hlat;
//       basicdata1.tf=my_tf;
//       basicdata1.gflux=my_gflux;
//       blk1_1.cdsnow=my1_cdsnow;
//    }

// Function to Set grid spacing, node depths, material type, density, 
// water content and explicit thermal properties
    basicf.grid_(&maxm,cmnt); //
     
// Function to intialize Soil Temp profile
    basicf.tinitl_(cmnt, tmax, tmin, tave, longmean, newtemp, profileflag);  // From Bailu   
     

    imaxp = basicdata1.imax + 1;
    dep1 = basicdata1.x[9];
    dep2 = basicdata1.x[11];
    dep3 = basicdata1.x[14];
    dep4 = basicdata1.x[15];
    dep5 = basicdata1.x[16];
    dep6 = basicdata1.x[17];
    dep7 = basicdata1.x[18];
    dep8 = basicdata1.x[19];
    dep9 = basicdata1.x[20];
    dep10 = basicdata1.x[24];
    dep11 = basicdata1.x[59];
    dep12 = basicdata1.x[64];
    dep13 = basicdata1.x[69];
    dep14 = basicdata1.x[74];
    dep15 = basicdata1.x[79];
    dep16 = basicdata1.x[84];
    dep17 = basicdata1.x[89];
    dep18 = basicdata1.x[94];
    kfor = kswitch + 10;
    
// Function to assign top boundary condition
    bctop_(&maxm, taer, &ktop, &sflx1, &sflx2, &sflx3);
// Function to compute snow fall
    basicf.snofal_(snow, sden, &comp, &eta0, &denfac, &fcmelt,cmnt);    
// Function to assign internal heat sources
	  basicf.sethet_(&maxm, kint, heat, heatt,cmnt);

    bcbot_(&maxm, tber);

//    if (kswitch == 0) {
// Function to intialize Soil Temp profile
	  //basicf.tinitl_(cmnt); //, tmax, tmin, tave, longmean, newtemp, profileflag  //From Bailu

//    }

/* --- INITIALIZATION */
    basicdata1.is = *is9;
    basicdata1.ism1 = *ism19;
    basicdata1.smass = *smass9;
    if (kswitch != 0) {
	     goto L70;
    }
    basicdata1.is = basicdata1.ig;
    basicdata1.ism1 = basicdata1.igm1;
    basicdata1.smass = 0.f;
    i__1 = basicdata1.igm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
      basicdata1.xfa[i__ - 1] = -1e10f;
	    basicdata1.xfb[i__ - 1] = 0.f;
	    basicdata1.water[i__ - 1] = 1.f;
	    basicdata1.dx[i__ - 1] = 0.f;
	    weight[i__ - 1] = 0.f;
    }
    weight[basicdata1.ig - 1] = 0.f;

/* --- SET XFA VALUE FOR INITIAL PHASE BOUNDARY */
    basicdata1.dx[basicdata1.imax - 1] = 0.f;
    
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <=basicdata1.imax; ++i__) {
      dxx = basicdata1.dx[i__ - 1];
      tu = basicdata1.t[i__ - 1] - basicdata1.tf;

	      if (i__ < basicdata1.imax) {
	          td = basicdata1.t[i__] - basicdata1.tf;
        }
	      basicdata1.xfa[i__ - 1] = -1e10f;
        basicdata1.xfb[i__ - 1] = 0.f;

        if (tu * td >= 0.f) {
	         goto L51;
	      }

        if (abs(tu - td) < 0.0001)  {tu = td + 0.01; } // to protect
        g = dxx * tu / (tu - td);
	      if (basicdata1.hlat != 0.f) {
	         basicdata1.xfa[i__ - 1] = g;
	      }
L51:
	;
    }
L70:
    i__1 = basicdata1.imax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	      basicdata1.ht[i__ - 1] = 0.f;
	      thalt[i__ - 1] = -999.9f;
	      basicdata1.htold[i__ - 1] = 0.f;
    }

/* --- CALCULATE AIR FREcmntING & THAWING INDEX */
    frza = 0.f;
    thwa = 0.f;
    tairlo = 1e3f;
    tairhi = -1e3f;
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      basicdata1.tair = taer[m - 1];
	      if (basicdata1.tair >= tairhi) {
	         tairhi = basicdata1.tair;
	      }
	      if (basicdata1.tair <= tairlo) {
	         tairlo = basicdata1.tair;
	      }
	      if (basicdata1.tair >= 0.f) {
	         goto L202;
     	  } else {
	         goto L201;
	}
L201:
	frza -= basicdata1.tair;
	goto L200;
 
L202:
	thwa += basicdata1.tair;
L200:
	;
    }
    thwa *= basicdata1.dtday;
    frza *= basicdata1.dtday;


//* ****************** BEGIN ANNUAL CYCLE # NAN *************************
    basicdata1.nan = 0;
    --nanmax;
L100:
    ++basicdata1.nan;
    if (basicdata1.nan > nanmax) {
	     kallyr = 1;
    }
    tslo = 1e3f;
    tshi = -1e3f;
    tgrlo = 1e3f;
    tgrhi = -1e3f;
    thws = 0.f;
    frzs = 0.f;
    thwg = 0.f;
    frzg = 0.f;
    tairan = 0.f;
    tgran = 0.f;
    nsum = 0;
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      thi[i__ - 1] = -1e3f;
	      tlo[i__ - 1] = 1e3f;
	      tanee[i__ - 1] = 0.f;
    }

/* ****************** BEGIN TIME STEPPING ********************************
 */
    kkk = (integer) basicdata1.first;
    basicdata1.time = 0.f;
    m = 0;
/* --- TIME = ELAPSED TIME IN DAYS AT END OF STEP */
L10:
    basicdata1.time += basicdata1.dtday;
    ++m;
    calend = basicdata1.first + basicdata1.time;
    if (calend > basicdata1.per) {
	     calend -= basicdata1.per;
    }

/* ---EVALUATE NEW AIR TEMPERATURE (TAIR) */
    basicdata1.tair = taer[m - 1];
    tdrive = basicdata1.tair;

/* --- EVALUATE NEW SURFACE CONDITIONS */
    surf_(&ktop, &sflx1, &sflx2, &sflx3, &tdrive);
    basicdata1.dtfaz = basicdata1.dt;

/* --- SNOW COVER GRID */
    acum = snow[m - 1];
    if (basicdata1.tair > 0.f) {
	     acum = 0.f;
    }
    snoden = sden[m - 1];

    neige_(&acum, &snoden, weight, &comp, &eta0, &denfac, &fcmelt, &tdrive);
    topold = basicdata1.t[basicdata1.is - 1];
    basicdata1.ism1 = basicdata1.is - 1;
    isp1 = basicdata1.is + 1;
    basicdata1.t[basicdata1.ism1 - 1] = tdrive;
    basicdata1.e[basicdata1.ism1 - 1] = tdrive;
    basicdata1.s[basicdata1.ism1 - 1] = 0.f;


/* --- UPDATE INTERNAL HEAT SOURCES */
    for (iheat = 1; iheat <= 2; ++iheat) {
	      heatt[iheat - 1] = heat[m + iheat * 3000 - 3001];
    }
    inthet_(kint, heatt);

/* --- EVALUATE NEW BOTTOM CONDITIONS */
    basicdata1.tbot = tber[m - 1];
    basicdata1.t[basicdata1.imax - 1] = basicdata1.tbot;

/* --- EVALUATE THERMAL PROPERTIES */
    tt = basicdata1.t[basicdata1.is - 1];
    i__1 = basicdata1.imax1;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {
	      dxx = basicdata1.dx[i__- 1];
	      td = basicdata1.t[i__];

	      pram_(&basicdata1.mater[i__ - 1], &tt, &td, &basicdata1.ddry[i__ - 1], &basicdata1.water[i__ - 1], &sph, &cnd);
	      tt = td;

        if (abs(dxx) < 0.01)  { dxx =0.01; };// to protect

	      basicdata1.conx[i__ - 1] = cnd / dxx;
	      basicdata1.capx[i__ - 1] = basicdata1.ddry[i__ - 1] * sph * dxx;
    }
    basicdata1.conx[basicdata1.ism1 - 1] = basicdata1.htop;
    basicdata1.capx[basicdata1.ism1 - 1] = 0.f;
    if (basicdata1.hlat > 0.f) {
	     goto L300;
    }
/* --- NO LATENT HEAT FOR ANY ELEMENT */

    asmblo_(&basicdata1.is);
    tblo_(&basicdata1.is);

    goto L60;

/* --- LOCATE UPPER & LOWER PHASE PLANES */
L300:

    i__1 = basicdata1.imax;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {
	      itzero = i__;
        if (basicdata1.xfa[i__ - 1] >= 0.f) {
	         goto L310;
	      }
    }
L310:
    jj = basicdata1.ism1 + basicdata1.imax1;
    i__1 = basicdata1.imax1;
    for (j = basicdata1.ism1; j <= i__1; ++j) {
	      i__ = jj - j;
	      ibzero = i__;
	      if (basicdata1.xfa[i__ - 1] >= 0.f) {
	         goto L320;
	      }
    }
L320:
    i__1 = ibzero - itzero;
    if (i__1 < 0) {
	     goto L33;
    } else if (i__1 == 0) {
	     goto L34;
    } else {
	     goto L35;
    }

/* --- NO PHASE PLANE */
L33:

    asmblo_(&basicdata1.is);

    topnew = tdrive * basicdata1.s[basicdata1.is - 1] + basicdata1.e[basicdata1.is -	1];
    if ((topnew - basicdata1.tf) * (topold - basicdata1.tf) < 0.f) {
	     goto L40;
    }
    tblo_(&basicdata1.is);
    goto L60;

/* --- SINGLE PHASE PLANE */
L34:
    if (basicdata1.xfb[itzero - 1] > 0.f) {
	     goto L35;
    }

    asmone_(&itzero);
    topnew = basicdata1.t[basicdata1.is - 1];
    goto L42;

/* --- MORE THAN ONE PHASE PLANE */
  phase_count =0.0;
L35:
//  if (phase_count < 1e6f) {  // to deal with multiple phase, 17/Feb/2000
    asmtwo_(&itzero, &ibzero);
//    phase_count = phase_count+1;
//    }
// 	else goto L60;
    topnew = basicdata1.t[basicdata1.is - 1];
    goto L42;

/* --- START OF SINGLE PHASE PLANE AT SURFACE */
L40:
    itzero = basicdata1.is;
    ibzero = basicdata1.is;
    basicdata1.xfa[itzero - 1] = 0.f;
    if (basicdata1.per > 0.f) {
	     goto L41;
    }
/* --- TRANSIENT PROBLEMS ONLY */
    topold = topnew;
    basicdata1.t[itzero - 1] = topold;
    goto L34;
L41:

    if (abs(topold - topnew) < 0.0001) { topold =  topnew + 0.01;}; // to protect

    basicdata1.dtfaz = basicdata1.dt * (topold - basicdata1.tf) / (topold - topnew);
    itz1 = itzero + 1;
    basicdata1.t[itzero - 1] = basicdata1.tf;

    asmblo_(&itz1);
    tblo_(&itz1);

    r__1 = topnew - basicdata1.tf;
    if (r__1 >=0)   topold = basicdata1.tf + c_b127;
    else   topold = basicdata1.tf - c_b127;

    basicdata1.t[itzero - 1] = topold;
    basicdata1.dtfaz = basicdata1.dt - basicdata1.dtfaz;
    goto L34;

/* --- START OF SECOND PHASE PLANE AT SURFACE */
/* --- TEMPERATURES SET TO ZERO BETWEEN PHASE PLANES */

L42:
	++debug1;
	if(debug1 >= 50){
		debug1=0;
		topnew=0;
		topold=0;
		goto L60;
	}
    if ((topnew - basicdata1.tf) * (topold - basicdata1.tf) >= 0.f) {
	     goto L60;
    }
    r__1 = topnew - basicdata1.tf;

    if (r__1 >=0) {tset = basicdata1.tf - c_b127;}
    else tset = basicdata1.tf + c_b127;

    i__1 = itzero;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {

	      basicdata1.e[i__ - 1] = 0.f;
	      basicdata1.s[i__ - 1] = 0.f;
	      basicdata1.t[i__ - 1] = tset;       // ***** should be careful .........
    }

    if (abs(topnew - topold)< 0.0001) { topnew =  topold + 0.01;}; // to protect

    basicdata1.dtfaz = basicdata1.dt * (topnew - basicdata1.tf) / (topnew - topold);
    ibzero = itzero;
/* --- TO ENSURE 2ND CALL TO ASMTWO IGNORES LOWER PHASE PLANE */
/*    IN CASE PREVIOUSLY HAD MULTIPLE PHASE PLANES */
    if (ibzero > itzero) {
       ibzero = 1;
    }
    itzero = basicdata1.is;
    r__1 = topnew - basicdata1.tf;
    if (r__1 >=0) topold = basicdata1.tf + c_b127;
    else topold = basicdata1.tf - c_b127;

    basicdata1.t[basicdata1.is - 1] = topold;
    if (basicdata1.xfa[basicdata1.is - 1] > 0.f) {
	     basicdata1.xfb[basicdata1.is - 1] = basicdata1.xfa[basicdata1.is - 1];
    }
    basicdata1.xfa[basicdata1.is - 1] = 0.f;

    goto L35;


/* **** CALCULATE RESULTS FOR TIME STEP M */
L60:
/* --- CALCULATE SURFACE FREcmntING AND THAWING INDICES */
    tsurf = basicdata1.t[basicdata1.is - 1];
    if (tsurf >= 0.f) {
	     goto L205;
    } else {
	     goto L204;
    }
L204:
    frzs -= tsurf;
    goto L206;
L205:
    thws += tsurf;
L206:

/* --- CALCULATE GROUND FREcmntING AND THAWING INDICES */
    tgrnd = basicdata1.t[basicdata1.ig - 1];
    if (tgrnd >= 0.f) {
       goto L208;
    } else {
	     goto L207;
    }
L207:
    frzg -= tgrnd;
    goto L209;
L208:
    thwg += tgrnd;
L209:

/* --- CALCULATE MAXIMUM AND MINIMUM TEMPERATURES */
    if (tsurf >= tshi) {
	tshi = tsurf;
    }
    if (tsurf <= tslo) {
	tslo = tsurf;
    }
    if (tgrnd >= tgrhi) {
	tgrhi = tgrnd;
    }
    if (tgrnd <= tgrlo) {
	tgrlo = tgrnd;
    }

    if (abs(basicdata1.mmax) < 0.001) { basicdata1.mmax = 5; }; // to protect
    
    tairan += taer[m - 1] / basicdata1.mmax;
    tgran += tgrnd / basicdata1.mmax;
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      if (basicdata1.t[i__ - 1] >= thi[i__ - 1]) {
	         thi[i__ - 1] = basicdata1.t[i__ - 1];
	      }
	      if (basicdata1.t[i__ - 1] <= tlo[i__ - 1]) {
	         tlo[i__ - 1] = basicdata1.t[i__ - 1];
        }
 //  printf(" ## %10.4f ", basicdata1.t[i__ - 1]);
 //  getch();
 		if(basicdata1.t[i__ - 1] > 50.0) { basicdata1.t[i__ - 1] = 50.0;}
		if(basicdata1.t[i__ - 1] < -50.0) { basicdata1.t[i__ - 1] = -50.0;}
    tanee[i__ - 1] += basicdata1.t[i__ - 1]/ (basicdata1.mmax);
    }

/* --- OUTPUT */
    if (kallyr == 0) {
	     goto L1000;
    }
    ++nsum;
    if (nsum < basicdata1.nst) {
	     goto L1000;
    }
    nsum = 0;
    ++kkk;

L1000:
    if (m < basicdata1.mmax) {
	     goto L10;
    }


/* ****************** ANNUAL CYCLE  NAN COMPLETE ***********************
*/

    thws *= basicdata1.dtday;
    frzs *= basicdata1.dtday;
    thwg *= basicdata1.dtday;
    frzg *= basicdata1.dtday;
    tdif = tgran - tairan;

    if (kenvl == 0) {
	     goto L9901;
    }
    if (kallyr == 0 && basicdata1.nan < nanmax) {
	     goto L9901;
    }
/* --- EQUILIBRIUM CHECK */
L9901:
    if (basicdata1.nan <= nanmax) {
	     goto L9902;
    }

// Output the active layer, including front frozen depth (frontd) and deep forzen depth (deepd)
first_phase=0;
int i___max = 0;
for (i__= 9; i__ <37; i__++) {
    if (tanee[i__] <basicdata1.tf) {
 			 if (i__ == 9) { // frozen top layer (example winter temperature profile: -10,-8,-6,...)
          *frontd =0.0; i___max = 9; /*printf("sthermal 792: frontd = %.3f \n", *frontd);*/}
       else {
            *frontd = (-1.0) * (basicdata1.x[i__-1]+ (basicdata1.x[i__]-basicdata1.x[i__-1])*((tanee[i__-1]-basicdata1.tf)/(tanee[i__-1]-tanee[i__]))); //printf("sthermal 794: basicdata1.xi-1 = %.3f, basicdata1.xi=%.3f,taneei-1=%.3f,taneei=%.3f\n", basicdata1.x[i__-1],basicdata1.x[i__],tanee[i__-1],tanee[i__]);
       			i___max = i__;
//*frontd = (-1.0) * (basicdata1.x[i__-1]+ basicdata1.x[i__])/2.0;
            if (fabs(*frontd) < 0.01) *frontd =0.0;
             }
       first_phase=i__;
 if (( *frontd >0 ) || (fabs(*frontd) < 0.05)) *frontd =0.0;
       i__=37; //exit();
   }
  }
  //for (i__= 9; i__ <37; i__++) {
  //printf("sthermal 805: d = %.3f, t = %.3f\n", basicdata1.x[i__], tanee[i__]); }
 //printf("sthermal 804: frontd = %.3f \n", *frontd);
//if (i___max == 0) {*frontd = nan;} 
if (i___max == 0) {*frontd = -99.0;}
  //printf("qsoil 807: frontd = %.3f\n", *frontd);
second_phase=0;
*thawbegin=0.0;
if (!(first_phase ==0 )) {
   for (i__= first_phase+1; i__ <37; i__++) {
       if (tanee[i__] >basicdata1.tf) {
          *thawbegin=(-1.0) *basicdata1.x[i__];
          second_phase=i__;
          i__=37;
       }
   }
   *thawend=0.0;
if (!(second_phase ==0)) {
   for (i__= second_phase+1; i__ <37; i__++) {
       if (tanee[i__] <basicdata1.tf) {
          *thawend=(-1.0) *basicdata1.x[i__];
          i__=37;
   }
  }
 }
}

//end of the active layer module

//    *tsoil = (tanee[9] + tanee[10] + tanee[11] + tanee[12] + tanee[13])/ 5.f;
//  *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
//  tanee[11]*basicdata1.x[11] + tanee[12]*basicdata1.x[12]+ tanee[13]*basicdata1.x[13]
//  +tanee[14]*basicdata1.x[14]  +tanee[15]*basicdata1.x[15])/
//  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]+basicdata1.x[13]
//  +basicdata1.x[14]+ basicdata1.x[15]); //define the top 20cm is from the bottom of the moss layer, assuming the moss thickness 10cm

// for veg type 2, alpine tundra

  switch (cmnt)
     {
     case 2: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
              tanee[11]*basicdata1.x[11])  /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;
     case 3: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
             tanee[11]*basicdata1.x[11])  /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;
     case 4:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] + tanee[11]*basicdata1.x[11]+ tanee[12]*basicdata1.x[12]+ tanee[13]*basicdata1.x[13])
             /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]+basicdata1.x[13]);break;

     case 5: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
              tanee[11]*basicdata1.x[11])  /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;
     case 6: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
              tanee[11]*basicdata1.x[11])  /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;

     case 7: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10])  /  (basicdata1.x[9]+basicdata1.x[10]); break;

     case 8:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] + tanee[11]*basicdata1.x[11]+ tanee[12]*basicdata1.x[12])
             /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]);break;
     case 9:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] + tanee[11]*basicdata1.x[11]+ tanee[12]*basicdata1.x[12])
              /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]);break;

     case 10:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] + tanee[11]*basicdata1.x[11]+ tanee[12]*basicdata1.x[12])
             /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]);break;

     case 11: *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +
              tanee[11]*basicdata1.x[11])  /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;

     case 12:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] + tanee[11]*basicdata1.x[11]+ tanee[12]*basicdata1.x[12])
             /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]+basicdata1.x[12]);break;

     default:    *tsoil = (tanee[9]*basicdata1.x[9] + tanee[10]*basicdata1.x[10] +   tanee[11]*basicdata1.x[11])
              /  (basicdata1.x[9]+basicdata1.x[10]+basicdata1.x[11]); break;
     }


 // to get frozen index temperature; using tanee[11] or [10] to determine index 28/Nov/2001
  switch (cmnt)
    {
     case 2: diffsoilt[0]=tanee[10]; break;
     case 3: diffsoilt[0]=tanee[10]; break;
     case 4: diffsoilt[0]=tanee[11]; break;
     case 5: diffsoilt[0]=tanee[10]; break;
     case 6: diffsoilt[0]=tanee[10]; break;
     case 7: diffsoilt[0]=tanee[10]; break;
     case 8: diffsoilt[0]=tanee[11]; break;
     case 9: diffsoilt[0]=tanee[10]; break;
     case 10: diffsoilt[0]=tanee[10]; break;
     case 11: diffsoilt[0]=tanee[10]; break;
     case 12: diffsoilt[0]=tanee[10]; break;
     default: diffsoilt[0]=tanee[10]; break;
    }


	diffsoilt[0] = cal_out_soilt(tanee, basicdata1.x,0.05); // 5cm
	//diffsoilt[0] = cal_out_soilt(tanee, basicdata1.x,0.0); 
	diffsoilt[1] = cal_out_soilt(tanee, basicdata1.x,0.1); 	// 10cm
	diffsoilt[2] = cal_out_soilt(tanee, basicdata1.x,0.2); 	// 20cm
	diffsoilt[3] = cal_out_soilt(tanee, basicdata1.x,0.5); 	// 50cm
	diffsoilt[4] = cal_out_soilt(tanee, basicdata1.x,1); 	// 100cm
	diffsoilt[5] = cal_out_soilt(tanee, basicdata1.x,2); 	// 200cm
	diffsoilt[6] = cal_out_soilt(tanee, basicdata1.x,3); 	// 300cm  
	diffsoilt[7] = cal_out_soilt(tanee, basicdata1.x,5); 	// 500cm 
	diffsoilt[8] = cal_out_soilt(tanee, basicdata1.x,8); 	// 800cm 
	diffsoilt[9] = cal_out_soilt(tanee, basicdata1.x,10); 	// 1000cm 


	for (i__= 0; i__ <37; i__++) { 
		//printf("stm 667: inod = %i, depth = %.3f, temp = %.3f \n",  inod, basicdata1.x[inod], tanee[inod]);
	}
	int last_layer;
	for (i__ = 10; i__ < 210; i__ ++)
	{	
		//printf("stm 667: inod = %i, depth = %.3f, temp = %.3f \n",  inod, basicdata1.x[inod], tanee[inod]);
		if (basicdata1.x[i__] == 0) {
		   last_layer = i__;
		   break;
		}
	}
		//printf("676: i = 0, newtemp = %.3f\n", newtemp[0]);
	i___max = 0;
	for (i__ = 0; i__ < 24; i__ ++){	
		//if (inod == 0) {printf("i = 0, newtemp = %.3f \n", newtemp[inod]);}
		  if (basicf.DEPTH[cmnt][i__] <= basicdata1.x[last_layer - 1]){		
		     newtemp[i__] = cal_out_soilt(tanee, basicdata1.x, basicf.DEPTH[cmnt][i__]); 
		     profileflag = 1;
		     i___max = i__;
		//printf("inod = %i, depth = %.3f, temp = %.3f, x = %.3f \n",  inod, basicf.DEPTH[cmnt][inod], newtemp[inod], basicdata1.x[last_layer]);
		}
	}
		newtemp[i___max + 1] = newtemp[i___max];
		//printf("inod_max + 1 = %i, depth = %.3f, temp = %.3f, x = %.3f \n",  inod_max + 1, basicf.DEPTH[cmnt][inod_max + 1], newtemp[inod_max + 1], basicdata1.x[last_layer]);
		//printf("691: i = 0, newtemp = %.3f\n", newtemp[0]);
	//this part gives soil temperature at 20cm depth for Rh
	*tsoil = diffsoilt[2];

//    diffsoilt[0]=tanee[9]; // surface
/*    diffsoilt[1]=tanee[10];// depth 0.12m
    diffsoilt[2]=tanee[11]; // depth 0.22m
    diffsoilt[3]=tanee[12]; // depth 0.32m
    diffsoilt[4]=tanee[13]; // depth 0.42
    diffsoilt[5]=tanee[14]; // depth 0.52 */

  //  diffsoilt[6]=tanee[15]; // depth 0.62
  //  diffsoilt[7]=tanee[16]; // depth 0.72m

    *is9 = basicdata1.is;
    *ism19 = basicdata1.ism1;
    *smass9 = basicdata1.smass;
    for (i__ = 1; i__ <= 210; ++i__) {
	      x9[i__] = basicdata1.x[i__ - 1];
	      dx9[i__] = basicdata1.dx[i__ - 1];
	      xfa9[i__] = basicdata1.xfa[i__ - 1];
	      xfb9[i__] = basicdata1.xfb[i__ - 1];
	      water9[i__] = basicdata1.water[i__ - 1];
	      t9[i__] = basicdata1.t[i__ - 1];
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	      weight9[i__] = weight[i__ - 1];
    }
    goto L999;
    
L9902:
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      thalto = thalt[i__ - 1];
	      thalt[i__ - 1] = basicdata1.t[i__ - 1];
	      dif = (r__1 = thalt[i__ - 1], dabs(r__1)) - dabs(thalto);
	      adif[i__ - 1] = dabs(dif);
    }
    adifmx = 0.f;
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      if (adif[i__ - 1] > adifmx) {
	         adifmx = adif[i__ - 1];
	      }
    }
    if (adifmx > .005f) {
	     goto L9904;
    }
    nanmax = 0;
    goto L100;
    
L9904:
    if (basicdata1.nan <= nanmax) {
	     goto L100;
    }

L999:
  //  initfile.close();
//    snowdepth.close();
//    kswitch=1;

    return 0;
} /* soiltemp_ */


// for surface heat flux
int surf_(integer *ktop, double *sflx1, double *sflx2, double *sflx3, double *tdrive)
{
    static double turb= 0.0, time0= 0.0, onoff2= 0.0, onoff3= 0.0, sigta3= 0.0, tairab= 0.0, rlwfix= 0.0, rad= 0.0;

    switch (*ktop) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
    }

L10:
    return 0;
L20:
    onoff2 = 0.f;
    onoff3 = 0.f;
    if (basicdata1.time - basicdata1.first < basicdata1.per / 2.f) {
       onoff2 = 1.f;
    }
    if (basicdata1.time - basicdata1.first >= basicdata1.per / 2.f) {
      	onoff3 = 1.f;
    }
    basicdata1.sflux = *sflx1 + onoff2 * *sflx2 - onoff3 * *sflx3;
    *tdrive = basicdata1.tf + onoff2 - onoff3;
    return 0;

L30:
/* --- LINEARIZED HEAT BALANCE TREATMENT */
 if (abs(basicdata1.per) < 0.0001) {  basicdata1.per = 0.01; }; // to protect

    rad = cos((basicdata1.time - time0) * 6.28f / basicdata1.per) * 100.f + 150.f;
    basicdata1.emiss = .7f;

    tairab = basicdata1.tair + 273.15f;
    sigta3 = basicdata1.sigma * tairab * tairab * tairab;
    rlwfix = sigta3 * (basicdata1.emiss - 1.f) * tairab;
    basicdata1.sflux = rad + rlwfix;

/* --- EVALUATE HEAT TRANSFER COEFFICIENT   TURB */
    turb = 6.f;
    basicdata1.htop = turb + sigta3 * 4.f;
    return 0;
} /* surf_ */

int inthet_(integer *kint, double *source)
{
    integer i__1;

    static double htip1= 0.0;
    static integer i__= 0, nheat= 0, ip1= 0;
    static double hti= 0.0;

    --source;
    --kint;

    i__1 = basicdata1.imax;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {
	      basicdata1.htold[i__ - 1] = basicdata1.ht[i__ - 1];
	      basicdata1.ht[i__ - 1] = 0.f;
    }

/* ---ASSEMBLE NEW NODAL HEATING VECTOR */
    nheat = 0;
    i__1 = basicdata1.imax1;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      ip1 = i__ + 1;
	      if (basicdata1.xhet[i__ - 1] == 999.9f) {
	         goto L20;
	      }
/* --- ELEMENT CONTAINS A POINT HEAT SOURCE */
	      ++nheat;
	      if (nheat > 2) {
	         return 0;
	      }
	      if (kint[nheat] != 1) {
	         goto L20;
	      }
	      if (basicdata1.t[i__ - 1] < -10.f && basicdata1.t[i__] < -10.f) {
	         source[nheat] = 0.f;
	      }

        if (abs(basicdata1.dx[i__ - 1]) < 0.0001) { basicdata1.dx[i__ - 1]=0.01;}; // to protect

	hti = source[nheat] * (basicdata1.x[ip1 - 1] - basicdata1.xhet[i__ - 1]) /
		 basicdata1.dx[i__ - 1];
	htip1 = source[nheat] * (basicdata1.xhet[i__ - 1] - basicdata1.x[i__ - 1])
		 / basicdata1.dx[i__ - 1];
	basicdata1.ht[i__ - 1] += hti;
	basicdata1.ht[ip1 - 1] += htip1;
L20:
	;
    }
    return 0;
} /* inthet_ */


int Soilthermal::grid_(integer *maxm, int& cmnt) //, double& vsm1, double& vsm2, double& vsm3, double& pthick
{
    integer i__1= 0;
    double r__1 = 0.0;
    doubledouble d__1= 0.0, d__2= 0.0;

    extern  int pram_(integer *, double *, double *, double *, double *, double *, double *);
    static double tempsoil = 0.0, rwas = 0.0, fwas = 0.0, vsph = 0.0, rnew = 0.0, f = 0.0;
    static integer i__ = 0, j = 0, k = 0, n = 0;
    static double r__ = 0.0, dtbig = 0.0, dense = 0.0, thick = 0.0, dtgiv = 0.0;
    static integer layer = 0;
    static double vdens = 0.0, vcond = 0.0, total = 0.0, dxmin = 0.0, dtmax = 0.0;
    static integer nodes = 0;
    static double ratio = 0.0, power = 0.0, df = 0.0;
    static integer nn =0;
    static double xx = 0.0, vspace = 0.0, dxagiv = 0.0, thmdxa = 0.0;
    static integer lmimax = 0;
    static double usedxa = 0.0, epsmin = 0.0, cnd = 0.0, dxa = 0.0, dxb = 0.0;
    static integer mat = 0;
    static double eps = 0.0, err = 0.0, sph = 0.0, wat = 0.0, top = 0.0, dxx = 0.0, comment[25];

	  for (i__1 = 0; i__1 < 25; i__1 ++){
	      comment[i__1] = 0.0;	
	  }
    lmimax = 210;
    --lmimax;
    
    if (kswitch == 0) {
//    snowdepth >> dummy >> dummy ;
//    snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy;
//    snowdepth >> basicdata1.first>>basicdata1.final>>basicdata1.per>>basicdata1.dtday>>basicdata1.theta;
      basicdata1.first= basicf.FIRST[cmnt];
      basicdata1.final= basicf.FINAL[cmnt];
      basicdata1.per = basicf.PER[cmnt];
      basicdata1.dtday= basicf.DTDAY[cmnt];
      basicdata1.theta = basicf.THETA[cmnt];

/*     my_first=basicdata1.first;
     my_final=basicdata1.final;
     my_per=basicdata1.per;
     my_dtday=basicdata1.dtday;
     my_theta=basicdata1.theta; */
   //  snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy;
   //  snowdepth >>top>>basicdata1.ig>>epsmin>>vspace>>vdens;
     top= basicf.TOP[cmnt];
     basicdata1.ig= basicf.IG[cmnt];
     epsmin=basicf.EPSMIN[cmnt];
     vspace= basicf.VSPACE[cmnt];
     vdens= basicf.VDEN[cmnt];

     my_ig= basicdata1.ig;

    // cout << endl << basicdata1.ig <<endl;
    mytop=top; myepsmin = epsmin; myvspace= vspace; myvdens=vspace;
    }
    else { top=mytop; epsmin = myepsmin; vspace= myvspace; vdens=myvspace;
     basicdata1.first=my_first;
     basicdata1.final=my_final;
     basicdata1.per=my_per;
     basicdata1.dtday=my_dtday;
     basicdata1.theta=my_theta;
     basicdata1.ig=my_ig;
      }

    if (basicdata1.theta <= .5f) {
	     basicdata1.theta = .5f;
    }
    basicdata1.theta *= 2.f;
    basicdata1.theta1 = 2.f - basicdata1.theta;
    basicdata1.igm1 = basicdata1.ig - 1;
    if (basicdata1.per <= 0.f && basicdata1.final <= basicdata1.first) {
       return 0;
    }
    if (basicdata1.per > 0.f && basicdata1.final <= 0.f) {
	     basicdata1.final = basicdata1.first + basicdata1.per;
    }
    total = basicdata1.final - basicdata1.first;

   if (abs(basicdata1.mmax) < 0.0001) {  basicdata1.mmax = 5; }; // to protect
   if (abs(basicdata1.dtday) < 0.0001) { basicdata1.dtday = 0.5; }; // to protect

    basicdata1.mmax = integer (total / basicdata1.dtday);
    basicdata1.dtday = total / basicdata1.mmax;
    dtgiv = basicdata1.dtday;
    xx = top * vspace;
    basicdata1.x[basicdata1.ig - 1] = xx;
    i__ = basicdata1.ig - 1;
    layer = 0;
    dtbig = 1e10f;
    
L10:
    ++layer;
    if (kswitch ==0) {
//     snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >>dummy;
//     snowdepth >> thick>>dxa>>dxb>>mat>>dense>>wat;
    if (layer==1) {
     thick= basicf.THICK1[cmnt]; dxa = basicf.DXA1[cmnt]; dxb= basicf.DXB1[cmnt];mat= basicf.MAT1[cmnt];dense=basicf.DENSE1[cmnt]; wat=basicf.WATER1[cmnt];
     //wat = vsm1;
    }
    if (layer==2) {
     thick= basicf.THICK2[cmnt]; dxa = basicf.DXA2[cmnt]; dxb= basicf.DXB2[cmnt];mat= basicf.MAT2[cmnt];dense=basicf.DENSE2[cmnt]; wat=basicf.WATER2[cmnt];
     //wat = vsm1; 
    }
/*    if (layer==3) {
     	 if(pthick > basicf.THICK2[cmnt]){
		     thick= pthick - basicf.THICK2[cmnt];
		     dxa = basicf.DXA2[cmnt]; dxb= basicf.DXB2[cmnt];mat= basicf.MAT2[cmnt];dense=basicf.DENSE2[cmnt]; wat=basicf.WATER2[cmnt];
		   }
		   else {
		        thick= basicf.THICK3[cmnt]; dxa = basicf.DXA3[cmnt]; dxb= basicf.DXB3[cmnt];mat= basicf.MAT3[cmnt];dense=basicf.DENSE3[cmnt]; wat=basicf.WATER3[cmnt];
		        wat = vsm3;
		   }
    } */
    if (layer==3) {
      thick= basicf.THICK3[cmnt]; dxa = basicf.DXA3[cmnt]; dxb= basicf.DXB3[cmnt];mat= basicf.MAT3[cmnt];dense=basicf.DENSE3[cmnt]; wat=basicf.WATER3[cmnt];
    }
   if (layer==4) {
     thick= basicf.THICK4[cmnt]; dxa = basicf.DXA4[cmnt]; dxb= basicf.DXB4[cmnt];mat= basicf.MAT4[cmnt];dense=basicf.DENSE4[cmnt]; wat=basicf.WATER4[cmnt];
    }
   if (layer==5) {
     thick= basicf.THICK5[cmnt]; dxa = basicf.DXA5[cmnt]; dxb= basicf.DXB5[cmnt];mat= basicf.MAT5[cmnt];dense=basicf.DENSE5[cmnt]; wat=basicf.WATER5[cmnt];
    }
   if (layer==6) {
     thick= basicf.THICK6[cmnt]; dxa = basicf.DXA6[cmnt]; dxb= basicf.DXB6[cmnt];mat= basicf.MAT6[cmnt];dense=basicf.DENSE6[cmnt]; wat=basicf.WATER6[cmnt];
    }


    // cout << endl <<wat <<endl;
    mythick[layer-1]=thick;mydxa[layer-1]= dxa; mydxb[layer-1]= dxb;
    mymat[layer-1]=mat;mydense[layer-1]=dense;
    mywat[layer-1]=wat;
   }
    else
     {thick=mythick[layer-1];dxa= mydxa[layer-1]; dxb= mydxb[layer-1];
     mat=mymat[layer-1];dense=mydense[layer-1];

   // adding for using soil water from TEM, 29/March/2001
   if ((layer ==3))
     {
         wat = blkw_.soilwater2;
//          printf(" %5.2f", wat);
     }
    else wat=mywat[layer-1];
 
     wat=mywat[layer-1];

    }

    if (thick <= 0.f) {goto L50;}

    thick = vspace * thick;
//        printf(" $$ = %5.2f", dxa);

    dxa = vspace * dxa;
    dxb = vspace * dxb;
    dense = vdens * dense;

    if (thick < dxa) { return 0;}

    if (mat > 10) { goto L12;}

 if (kswitch ==0)  {
//  snowdepth >> dummy >> dummy >> dummy >> dummy >> dummy >>dummy;
//  snowdepth >> vcond >> vsph >>basicdata1.condt[mat - 1]>>basicdata1.spht[mat - 1]>>basicdata1.condf[mat - 1]>>basicdata1.sphf[mat - 1];
    if (layer==1) {
    vcond=basicf.vcond1[cmnt]; vsph=basicf.vsph1[cmnt]; basicdata1.condt[mat - 1]=basicf.cond1[cmnt];
    basicdata1.spht[mat - 1]=basicf.spht1[cmnt]; basicdata1.condf[mat - 1]=basicf.condf1[cmnt];
    basicdata1.sphf[mat - 1]=basicf.sphf1[cmnt];
    }
    if (layer==2) {
    vcond=basicf.vcond2[cmnt]; vsph=basicf.vsph2[cmnt]; basicdata1.condt[mat - 1]=basicf.cond2[cmnt];
    basicdata1.spht[mat - 1]=basicf.spht2[cmnt]; basicdata1.condf[mat - 1]=basicf.condf2[cmnt];
    basicdata1.sphf[mat - 1]=basicf.sphf2[cmnt];
    }
    if (layer==3) {
    vcond=basicf.vcond3[cmnt]; vsph=basicf.vsph3[cmnt]; basicdata1.condt[mat - 1]=basicf.cond3[cmnt];
    basicdata1.spht[mat - 1]=basicf.spht3[cmnt]; basicdata1.condf[mat - 1]=basicf.condf3[cmnt];
    basicdata1.sphf[mat - 1]=basicf.sphf3[cmnt];
    }
    if (layer==4) {
    vcond=basicf.vcond4[cmnt]; vsph=basicf.vsph4[cmnt]; basicdata1.condt[mat - 1]=basicf.cond4[cmnt];
    basicdata1.spht[mat - 1]=basicf.spht4[cmnt]; basicdata1.condf[mat - 1]=basicf.condf4[cmnt];
    basicdata1.sphf[mat - 1]=basicf.sphf4[cmnt];
    }
    if (layer==5) {
    vcond=basicf.vcond5[cmnt]; vsph=basicf.vsph5[cmnt]; basicdata1.condt[mat - 1]=basicf.cond5[cmnt];
    basicdata1.spht[mat - 1]=basicf.spht5[cmnt]; basicdata1.condf[mat - 1]=basicf.condf5[cmnt];
    basicdata1.sphf[mat - 1]=basicf.sphf5[cmnt];
    }

  // cout <<endl << vsph<<endl;
  myvcond[layer-1]=vcond; myvsph[layer-1]=vsph;
  mycondt[mat - 1][layer-1]=basicdata1.condt[mat - 1];
  myspht[mat - 1][layer-1] =basicdata1.spht[mat - 1];
  mycondf[mat - 1][layer-1]=basicdata1.condf[mat - 1];
  mysphf[mat - 1][layer-1]=basicdata1.sphf[mat - 1];
   }
  else {vcond=myvcond[layer-1]; vsph=myvsph[layer-1];
  basicdata1.condt[mat - 1]=mycondt[mat - 1][layer-1];
  basicdata1.spht[mat - 1]=myspht[mat - 1][layer-1];
  basicdata1.condf[mat - 1]=mycondf[mat - 1][layer-1];
  basicdata1.sphf[mat - 1]=mysphf[mat - 1][layer-1];
    }

   if (vcond == 0.f) {
      vcond = 1.f;
   }
   if (vsph == 0.f) {
	     vsph = 1.f;
   }
   
    basicdata1.condt[mat - 1] = vcond * basicdata1.condt[mat - 1];
    basicdata1.spht[mat - 1] = vsph * basicdata1.spht[mat - 1];
    basicdata1.condf[mat - 1] = vcond * basicdata1.condf[mat - 1];
    basicdata1.sphf[mat - 1] = vsph * basicdata1.sphf[mat - 1];

    cnd = basicdata1.condf[mat - 1];
    sph = basicdata1.sphf[mat - 1];
    goto L14;

L12:
    if (mat > 20) {
	     goto L13;
    }
    pram_(&mat, &basicdata1.tf, &basicdata1.tf, &dense, &wat, &sph, &cnd);

    goto L14;

L13:
    if (mat > 30) {
       return 0;
    }
    tempsoil = basicdata1.tf + 1.f;
    for (n = 1; n <= 10; ++n) {
	      pram_(&mat, &tempsoil, &tempsoil, &dense, &wat, &sph, &cnd);
	      tempsoil += -1.f;
    }

/* --- CALCULATE MINIMUM ALLOWED GRID SPACING */
L14:
    double sph_dense;
    if ((sph * dense) ==0.0)  {  sph_dense = 0.01;}
    else  sph_dense = sph * dense;
//    eps = epsmin * 86400.f * cnd / (sph * dense);
    eps = epsmin * 86400.f * cnd / sph_dense;         // to protect
    dxagiv = dxa;
L20:
    if (( eps * basicdata1.dtday) < 0.0) {  eps =0.01; }; // to protect

    dxmin = sqrt(eps * basicdata1.dtday);
    if (dxa >= dxmin) {
	     goto L21;
    }
    dxa = dxmin;
    if (dxa < thick) {
	     goto L21;
    }
    dxa = dxagiv;
    basicdata1.dtday *= .5f;
    basicdata1.mmax <<= 1;
    goto L20;
L21:
    if (abs(eps) < 0.0001) { eps = 0.01;}; // to protect

    dtmax = dxa * dxa / eps;
    if (dtmax < dtbig) {
	     dtbig = dtmax;
    }

/* --- SET VALUES AT CALCULATION NODES FOR LAYER N */
    if (dxa < thick * .5f && dxb > dxa) {
	     goto L400;
    }
/* --- LINEAR GRID SPACING WITHIN LAYER */

    if (abs(dxa) < 0.0001) { dxa = 0.01;}; //to protect
    if (abs(nodes) < 0.00010) { nodes = 1;};
    if (abs(vspace) < 0.0001) { vspace = 0.01;};

    nodes = integer (thick / dxa + .01f);
    dxa = thick / nodes;
    usedxa = dxa / vspace;

    i__1 = nodes;
    for (j = 1; j <= i__1; ++j) {
	      ++i__;
	      if (i__ > lmimax) {
	         goto L50;
	      }
	      basicdata1.x[i__ - 1] = xx;
	      basicdata1.dx[i__ - 1] = dxa;
	      xx += dxa;
	      basicdata1.mater[i__ - 1] = mat;
	      basicdata1.ddry[i__ - 1] = dense;
	      basicdata1.water[i__ - 1] = wat;
    }
    goto L10;

L400:
    if (abs(vspace) < 0.0001) {  vspace = 0.01;};
    usedxa = dxa / vspace;
    if (abs(thick - dxb) < 0.0001) { thick = dxb + 0.01;}; // to protect
    r__ = (thick - dxa) / (thick - dxb);

    if (dxa ==0.0) { dxa = 0.01;}; //to protect
      ratio = dxb / dxa;

    if (r__ ==1.0) { r__ = 0.01;};  //to protect
    if (r__ < 0.0) {  r__=0.01;};
    if (abs(dxa) < 0.0001) {  dxa = 0.01;}; //to protect
    nn = integer (log(fabs(dxb / dxa)) * 10.f / log(r__) + 5.f);

    nn /= 10;
    nodes = nn + 1;
/* --- REVISE VALUE OF R BY SECANT METHOD */
    thmdxa = thick - dxa;
    power = 1.f / nn;
    d__1 = (doubledouble) ratio;
    d__2 = (doubledouble) power;
    rwas = pow(d__1,d__2);

    fwas = dxa * pow(rwas, (double)nodes) - rwas * thick + thmdxa;
    err = 1e-6f;
    for (k = 1; k <= 20; ++k) {
	      f = dxa * pow(r__, (double)nodes) - r__ * thick + thmdxa;
	      if (fabs(f) <= err) {
	         goto L402;
	      }
	      df = f - fwas;
	      if (df == 0.f) {
	         goto L402;
	      }
	      rnew = r__ - (r__ - rwas) * f / df;
	      rwas = r__;
	      fwas = f;
	      r__ = rnew;
    }
    
L402:
    r__ = rnew;
    dxx = dxa;
    i__1 = nodes;
    for (j = 1; j <= i__1; ++j) {
	      ++i__;
	      if (i__ > lmimax) {
	         goto L50;
	      }
	      basicdata1.x[i__ - 1] = xx;
	      basicdata1.dx[i__ - 1] = dxx;
	      xx += dxx;
	      dxx *= r__;
	      basicdata1.mater[i__ - 1] = mat;
	      basicdata1.ddry[i__ - 1] = dense;
	      basicdata1.water[i__ - 1] = wat;
    }
    goto L10;

L50:
    basicdata1.imax1 = i__;
    basicdata1.imax = i__ + 1;
    basicdata1.dt = basicdata1.dtday * 86400.f;
    i__1 = basicdata1.imax1;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	  if (i__ <= lmimax) {
	         goto L800;
	 }
   
L800:
 return 0;
	;
    }
    basicdata1.x[basicdata1.imax - 1] = xx;
    return 0;
} /* grid_ */


 int bctop_(integer *maxm, double *taer, integer *ktop, double *sflx1, double *sflx2, double *sflx3){
    integer i__1 = 0;
    extern  int data_(integer *, integer *, double *);
    static double tint = 0.0;
    static integer m = 0, index = 0, np = 0;
    static double thnfac = 0.0, fznfac = 0.0, fac = 0.0;

    --taer;
    *ktop = 1;
    basicdata1.htop = 1e20f;
    basicdata1.sflux = 0.f;
    index = 3;
    np = -99;
    data_(&index, &np, &taer[1]);

    thnfac = 1.f;
    fznfac = 1.f;
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      tint = taer[m];
	      fac = thnfac;
	      if (tint <= 0.f) {
	         fac = fznfac;
	      }
	  taer[m] = fac * tint;
    }
    return 0;
} /* bctop_ */


 int Soilthermal::sethet_(integer *maxm, integer *kint, double *heat, double *source, int& cmnt)
{
    integer heat_dim1 = 0, heat_offset = 0, i__1 = 0;
    extern  int data_(integer *, integer *, double *);
    static double vdep = 0.0;
    static integer i__ = 0, m = 0, iheat = 0, index = 0, kinti = 0;
    static double header[80] = {0.0};
    static integer np = 0;
    static double dephet = 0.0;
	  for (i__1 = 0; i__1 < 80; i__1 ++){
	      header[i__1] = 0.0;
	  }
      --source;
    heat_dim1 = *maxm;
    heat_offset = heat_dim1 + 1;
    heat -= heat_offset;
    --kint;

    for (iheat = 1; iheat <= 2; ++iheat) {
	      kint[iheat] = 1;
    }
    iheat = 0;
    ++iheat;
    if (kswitch ==0 ) {
//    snowdepth >> dummy >> dummy >> dummy;
//c-- read information on successive heat sources
 //   snowdepth >> dummy >> dummy >> dummy >>dummy>>dummy>>dummy;
 //   snowdepth >> dummy >> dummy >> dummy;
 //   snowdepth >> kinti >> vdep >> dephet;
       kinti= basicf.kint[cmnt];
       vdep=basicf.VDEP1[cmnt];
       dephet=basicf.DEPHET[cmnt]; // probably problem here ??

       mykinti =kinti; myvdep2=vdep; mydephet=dephet;
    // cout << " " << kinti << endl;
    }
    else { kinti =mykinti; vdep=myvdep2; dephet=mydephet;  }

    if (kinti == 99999) {
	     return 0;
    }
    dephet = vdep * dephet + basicdata1.x[basicdata1.ig - 1];
    kint[iheat] = kinti;

    i__1 = basicdata1.imax1;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      basicdata1.xhet[i__ - 1] = 999.9f;
	      if (basicdata1.x[i__ - 1] <= dephet && dephet <= basicdata1.x[i__]) {
	         goto L710;
	      }
	     goto L700;
           
L710:
	basicdata1.xhet[i__ - 1] = dephet;
 
L700:
	;
    }
    data_(&index, &np, &source[1]);
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      heat[m + iheat * heat_dim1] = source[m];
    }
    return 0;
} /* sethet_ */

 int bcbot_(integer *maxm, double *tber)
{
    integer i__1 = 0;

    static integer kbot = 0, m = 0;

    --tber;

    kbot = 2;
    basicdata1.hbot = 0.f;
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      tber[m] = 0.f;
    }
    return 0;
} /* bcbot_ */


 int data_(integer *index, integer *np, double *fitted)
{
    integer i__1 = 0;
    static double phas = 0.0;
    static integer nmax = 0;
    static double dayo = 0.0;
    static integer m = 0, n = 0;
    static double omega = 0.0, delta = 0.0, tmean = 0.0, elraw[1500], value = 0.0, vtime = 0.0;
    static integer nn = 0, no = 0;
    static double tampcs = 0.0, tampsn = 0.0;
    extern  int intrpl_(double *, double *, integer *, integer *,double *, double *, integer *);
    static integer nn2 = 0;
    static double cal = 0.0, pie = 0.0, day = 0.0, raw[1500], omt = 0.0, comment[30];
	  for (i__1 = 0; i__1 < 1500; i__1 ++){
	      elraw[i__1] = 0.0;
	      raw[i__1] = 0.0;	
	  }
	  for (i__1 = 0; i__1 < 30; i__1 ++){
	      comment[i__1] = 0.0;	
	  }
    --fitted;

    vtime = 1.f;
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      basicdata1.eltim[m - 1] = m * basicdata1.dtday;
	      cal = basicdata1.eltim[m - 1] + basicdata1.first;
	      if (cal > basicdata1.per) {
	         cal -= basicdata1.per;
	      }
	      basicdata1.caltim[m - 1] = cal;
    }
    switch (*index) {
	  case 1:  goto L1;
	  case 2:  goto L2;
	  case 3:  goto L3;
  	case 4:  goto L4;
    }
    
L1:
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
      	fitted[m] = value;
    }
    return 0;
L2:
    pie = 3.141593f;
    if ( abs(basicdata1.per) < 0.0001) { basicdata1.per = 0.01;}; // to protect
    if ( abs(basicdata1.mmax) < 0.0001) { basicdata1.mmax = 1;}; // to protect

    omega = pie * 2.f / basicdata1.per;
    i__1 = basicdata1.mmax;
    for (m = 1; m <= i__1; ++m) {
	      omt = omega * (basicdata1.eltim[m - 1] + phas);
	      fitted[m] = tmean + tampsn * sin(omt) + tampcs * cos(omt);
    }
    goto L600;
L3:
    elraw[0] = 0.f;
    raw[0] = blk2_1.airt1;
    elraw[1] = 15.f;
    raw[1] = blk2_1.airt2;
    elraw[2] = 30.f;
    raw[2] = blk2_1.airt3;
    nmax = 3;
    goto L500;
L4:
    n = 0;
    day = -10.f;
    delta = 0.f;
L500:
/* --- INTERPOLATE */
    intrpl_(elraw, raw, &nmax, np, basicdata1.eltim, &fitted[1], &basicdata1.mmax)
	    ;
L600:
    return 0;
} /* data_ */


int pram_(integer *mat, double *ta, double *tb, double *dense,double *wat, double *sph, double *cnd)
{
    doubledouble d__1 = 0.0, d__2 = 0.0;
    static double clca = 0.0, cice = 0.0, cacl = 0.0, diff = 0.0, csca = 0.0, cicl = 0.0;
    static integer mat10 = 0;
    static double cscl = 0.0, tempsoil = 0.0, exbm1 = 0.0, a = 0.0, b = 0.0, c__ = 0.0, f = 0.0, q = 0.0, r__ = 0.0, v = 0.0, cpore = 0.0, xwsat = 0.0, ca = 0.0,
	    fa = 0.0, ga = 0.0, ha = 0.0, gc = 0.0, ci = 0.0, cl = 0.0, fi = 0.0, fl = 0.0, hi = 0.0, hl = 0.0, cs = 0.0, du = 0.0, ri = 0.0, rl = 0.0, xa = 0.0, xi = 0.0,
	    fs = 0.0, hs = 0.0, xl = 0.0, rs = 0.0, ra = 0.0, cp = 0.0, xp = 0.0, vt = 0.0, xs = 0.0, fp = 0.0, hp = 0.0, xw = 0.0, am1 = 0.0, ap1 = 0.0, um1 = 0.0,
	    um2 = 0.0, um4 = 0.0, vt0 = 0.0, del = 0.0, gam = 0.0, tab = 0.0, pie = 0.0, unf = 0.0, sqt = 0.0, cica = 0.0;

    if (*mat <= 10) {
	     goto L1000;
    }
    tempsoil = (*ta + *tb) * .5f;
    if (*mat >= 31) {
	     goto L999;
    }
    mat10 = *mat - 10;
    if (*mat <= 20 && basicdata1.time > basicdata1.dtday) {
	     goto L1000;
    }
    switch (mat10) {
	case 1:  goto L11;
	case 2:  goto L12;
	case 3:  goto L13;
	case 4:  goto L14;
	case 5:  goto L15;
	case 6:  goto L16;
	case 7:  goto L17;
	case 8:  goto L18;
	case 9:  goto L19;
	case 10:  goto L20;
	case 11:  goto L21;
	case 12:  goto L22;
	case 13:  goto L23;
	case 14:  goto L24;
	case 15:  goto L25;
	case 16:  goto L26;
	case 17:  goto L27;
	case 18:  goto L28;
	case 19:  goto L29;
	case 20:  goto L30;
    }

/* ***********************************************************************
 */
L1000:
    if (*ta <= basicdata1.tf) {
       goto L1001;
    }
    *cnd = basicdata1.condt[*mat - 1];
    *sph = basicdata1.spht[*mat - 1];
    return 0;
L1001:
    *cnd = basicdata1.condf[*mat - 1];
    *sph = basicdata1.sphf[*mat - 1];
    return 0;

/* **********************************************************************
*/
L11:
    basicdata1.condt[*mat - 1] = .03f;
    basicdata1.spht[*mat - 1] = 1255.8f;
    basicdata1.condf[*mat - 1] = basicdata1.condt[*mat - 1];
    basicdata1.sphf[*mat - 1] = basicdata1.spht[*mat - 1];
    goto L1000;

/* **********************************************************************
*/
L12:
    basicdata1.condt[*mat - 1] = 0.f;
    if (basicdata1.condt[*mat - 1] == 0.f) {
       return 0;
    }
    goto L1000;

/* ***********************************************************************
 */
L13:
    basicdata1.sphf[*mat - 1] = *wat * 2093.f + 745.1f;
    d__1 = (doubledouble) (*dense * 8.12e-4f);
    d__2 = (doubledouble) (*dense * 9.12e-4f);
    basicdata1.condf[*mat - 1] = pow(c_b358, d__1) * .0109f + *wat * .46f *
    pow(c_b358, d__2);
/* --- THAWED */
    basicdata1.spht[*mat - 1] = *wat * 4186.f + 745.1f;
    d__1 = (doubledouble) (*dense * 6.25e-4f);

    //if ((*wat ==0.0) || (*wat < 0.0)) *wat =0.01; // 28/04/2001
	  if ((*wat < 0.0001)) *wat =0.01;
    basicdata1.condt[*mat - 1] = (log(*wat) * .044f + .26f) * pow(c_b358, d__1);
    goto L1000;


/* ***********************************************************************
 */
L14:
    basicdata1.sphf[*mat - 1] = *wat * 2093.f + 745.1f;
    d__1 = (doubledouble) (*dense * .001375f);
    d__2 = (doubledouble) (*dense * 5e-4f);
    basicdata1.condf[*mat - 1] = pow(c_b358, d__1) * .00143f + *wat *
	    1.22f * pow(c_b358, d__2);
/* --- THAWED */
    basicdata1.spht[*mat - 1] = *wat * 4186.f + 745.1f;
    d__1 = (doubledouble) (*dense * 6.25e-4f);

    //if ((*wat ==0.0) || (*wat < 0.0)) *wat =0.01; // 28/04/2001
	  if ((*wat < 0.0001)) *wat =0.01;
    basicdata1.condt[*mat - 1] = (log(*wat) * .0565f + .2304f) * pow(c_b358,
	     d__1);
    goto L1000;

L15:
    goto L999;
L16:
    goto L999;
L17:
    goto L999;
L18:
    goto L999;
L19:
    goto L999;
L20:
    goto L999;


/* ***********************************************************************
 */
L21:
    *sph = 2093.f;
    *cnd = blk1_1.cdsnow;
    if (*ta < 0.f) {
	     return 0;
    }
    *sph = 4086.f;
    *cnd = 10.f;
    return 0;
/* ---FROZEN */
    v = 6120.f;
    if (abs(tempsoil + 273.15f) < 0.01)  { tempsoil = 273.15f + 0.01;}; // to protect

    vt = v / (tempsoil + 273.15f);
    vt0 = v / 273.15f;
    cpore = (vt - 1.f) * .00110295f * exp(vt0 - vt) + .0242f;
    cice = 2.1f;
    pie = 3.1415926536f;
    q = *dense / 917.f;
    d__1 = (doubledouble) (q * 6.f / pie);
    r__ = pow(d__1, c_b365);
    if (r__ <= 1.f) {
	     goto L218;
    }
    *cnd = *dense * 2.9e-6f * *dense;
    return 0;
    
L218:
    if (abs(cpore) < 0.01) {  cpore = 0.01;}; // to protect

    a = sqrt(4.f / (pie * r__ * r__ * (cice / cpore - 1.f)) + 1.f);
    ap1 = a + 1.f;
    am1 = a - 1.f;
    del = r__ * .1f;
    gam = pie * del * del / 4.f;
    diff = cice - cpore;

    //if (diff ==0.0) diff =0.01; // 28/04/2001

    if (abs(am1) < 0.01) {  am1 = 0.01; } ; // to protect
    if (abs(diff) < 0.01) {  diff = 0.01; }
    if (am1 < 0.0)  {  am1 = 0.01; } ; // to protect

    f = (1.f - r__) / (diff * gam + cpore) + log(fabs(ap1 / am1)) * 2.f / (a * pie * r__ * diff);

    if (abs(f) < 0.01) {  f = 0.01;}; // to protect

    *cnd = 1.f / f;
    return 0;


/* ***********************************************************************
 */
L22:
    *cnd = 10.f;
    *sph = 4086.f;
    if (*ta >= 0.f) {
	     return 0;
    }
    *sph = 2093.f;
    *cnd = *dense * 2.9e-6f * *dense + .029f;
    return 0;
    
L23:
/* --- TYPE #23 */
    unf = *wat;
    du = 0.f;
    if (tempsoil >= basicdata1.tf) {
	     goto L235;
    }
    um1 = *wat * .15f;
    um4 = *wat * .01f;

    if (abs(*wat - um1) < 0.01) {   *wat = um1+ 0.01;}; // to protect
    exbm1 = (um1 - um4) / (*wat - um1);

    if (abs(1.f - exbm1) < 0.01)  { exbm1 = 1.f - 0.01;}; // to protect
    a = (*wat - um1) / (1.f - exbm1);

    if ((exbm1 <0.0) || (exbm1 ==0.0))  { exbm1 = 0.01;}; // to protect
    b = -log(exbm1);
    c__ = *wat - a;

    if ( tempsoil >= 0.0) {  tempsoil = -0.01;}; // to protect

    sqt = sqrt(-tempsoil);
    unf = a * exp(-b * sqt) + c__;
    du = b * .5f * (unf - c__) / sqt;
    if (unf - c__ > 1e-5f) {
	     goto L235;
    }
    unf = c__;
    du = 0.f;

/* --- APPARANT HEAT CAPACITY PER UNIT MASS */
L235:
    *sph = (*wat + unf) * 2093.f + 711.62f + du * 3.34e5f;
    rs = 2700.f;
    rl = 1e3f;
    ri = 917.f;
    xs = *dense / rs;
    xl = *dense * unf / rl;
    xi = *dense * (*wat - unf) / ri;
    xw = xl + xi;
    xwsat = 1.f - xs;
    xa = xwsat - xw;
    if (xw <= xwsat) {
	     goto L236;
    }
    xa = 0.f;
    xs = 1 - xw;
L236:
    cs = 7.f;
    cl = .56f;
    ci = 2.1f;
    if ( abs(tempsoil + 273.2f) < 0.0001 ) { tempsoil = -0.01;}; // to protect

    ca = exp(tempsoil * .082f) * 5.96f / (tempsoil + 273.2f) + .0242f;
    ga = .1f;
    gc = .1f;

    if ( abs(ca) < 0.01) { ca =0.01; }; // to protect

    csca = cs / ca;
	  if(csca == 1.0f) {csca = 1.1f;}    //to protect, zbl edit
    	if (abs(ga * (csca - 1.f) + 1.f) < 0.00001) {ga = (0.01 - 1.0f)/(csca - 1.f) ;}  // to protect, zbl edit
    	if (abs(gc * (csca - 1.f) + 1.f) < 0.00001) {gc = (0.01 - 1.0f)/(csca - 1.f) ;}  // to protect, zbl edit
    	if (abs(ga * (cica - 1.f) + 1.f) < 0.00001) {ga = (0.01 - 1.0f)/(cica - 1.f) ;}  // to protect, zbl edit
    	if (abs(gc * (cica - 1.f) + 1.f) < 0.00001) {gc = (0.01 - 1.0f)/(cica - 1.f) ;}  // to protect, zbl edit
    fs = .666667f / (ga * (csca - 1.f) + 1.f) + .333333f / (gc * (csca - 1.f)
	    + 1.f);
    cica = ci / ca;
    fi = .666667f / (ga * (cica - 1.f) + 1.f) + .333333f / (gc * (cica - 1.f)
	    + 1.f);
    clca = cl / ca;
    fl = .666667f / (ga * (clca - 1.f) + 1.f) + .333333f / (gc * (clca - 1.f)
	    + 1.f);
    hs = xs * fs;
    hi = xi * fi;
    hl = xl * fl;
    if (abs(xa + hs + hi + hl) < 0.01) {  xa = -(hs +hi + hl) + 0.01; } ; // to protect

    *cnd = (xa * ca + hs * cs + hi * ci + hl * cl) / (xa + hs + hi + hl);
    return 0;

/* --- TYPE #24 */
L24:
    unf = *wat;
    du = 0.f;
    if (tempsoil >= basicdata1.tf) {
	     goto L245;
    }
    um1 = *wat * .35f;
    um4 = *wat * .17f;
    if (abs(*wat - um1) < 0.000) {  *wat = um1 + 0.01; }; //to protect
    exbm1 = (um1 - um4) / (*wat - um1);

    if (abs(1.f - exbm1) < 0.0001) {  exbm1 = 0.01;}; // to protect

    a = (*wat - um1) / (1.f - exbm1);

    if (exbm1 <= 0.0) {  exbm1 = 0.01;}; // to protect

    b = -log(exbm1);
    c__ = *wat - a;

    if ( tempsoil >= 0.0) {  tempsoil = -0.01;}; // to protect

    sqt = sqrt(-tempsoil);
    unf = a * exp(-b * sqt) + c__;

    du = b * .5f * (unf - c__) / sqt;
    if (unf - c__ > 1e-5f) {
	     goto L245;
    }
    unf = c__;
    du = 0.f;

/* --- APPARANT HEAT CAPACITY PER UNIT MASS */
L245:
    *sph = (*wat + unf) * 2093.f + 711.62f + du * 3.34e5f;
    rs = 2700.f;
    rl = 1e3f;
    ri = 917.f;
    xs = *dense / rs;
    xl = *dense * unf / rl;
    xi = *dense * (*wat - unf) / ri;
    xw = xl + xi;
    xwsat = 1.f - xs;
    xa = xwsat - xw;
    if (xw <= xwsat) {
	     goto L246;
    }
    xa = 0.f;
    xs = 1 - xw;
L246:
    cs = 2.9f;
    cl = .56f;
    ci = 2.1f;
    ca = 0.f;
    ga = .125f;
    gc = .75f;
    cscl = cs / cl;

    double gaga = 0.0;
    double gcgc = 0.0;
    
	  if (abs(ga * (cscl - 1.f) + 1.f) < 0.0001){ gaga = (ga * (cscl - 1.f) + 1.f) + 0.01;}
	  else gaga = (ga * (cscl - 1.f) + 1.f);
	  if (abs(gc * (cscl - 1.f)+ 1.f) < 0.0001){ gcgc = (gc * (cscl - 1.f)+ 1.f)+ 0.01;}
	  else gcgc = (gc * (cscl - 1.f) + 1.f);

    fs = .666667f / gaga + .333333f / gcgc;

    cicl = ci / cl;

  //  fi = .666667f / (ga * (cicl - 1.f) + 1.f) + .333333f / (gc * (cicl - 1.f)    + 1.f);
    double gacic1 = 0.0;
    double gccic1 = 0.0;

    if (abs(ga * (cicl - 1.f) + 1.f) < 0.0001) { gacic1 = (ga * (cicl - 1.f) + 1.f) + 0.01;}
    else gacic1 = (ga * (cicl - 1.f) + 1.f);

    if (abs(gc * (cicl - 1.f)    + 1.f) < 0.0001) { gccic1 = (gc * (cicl - 1.f)    + 1.f) + 0.01; }
    else gccic1 = (gc * (cicl - 1.f)    + 1.f);

    fi = .666667f / gacic1 + .333333f / gccic1;

    fa = 0.f;
    hs = xs * fs;
    hi = xi * fi;
    ha = xa * fa;
    if (abs(xl + hs + hi + ha) < 0.0001) {  xl = -(hs +hi + ha) + 0.01; } ; // to protect

    *cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / (xl + hs + hi + ha);
    return 0;
    
L25:
    unf = *wat;
    du = 0.f;
    tab = tempsoil + 273.15f;

    if (tab ==0.0) {  tab =0.01;}; // to protect

    v = 5400.f;
    ca = (v / tab - 1.f) * .00110295f * exp(v * (.0036609921288669233f - 1.f /
	     tab)) + .0242f;
    if (tempsoil >= 0.f) {
	     goto L255;
    }
    v = 6120.f;

    if (abs(tab) < 0.0001) {  tab =0.01;}; // to protect
    ca = (v / tab - 1.f) * .00110295f * exp(v * (.0036609921288669233f - 1.f /
	     tab)) + .0242f;
    um1 = .57998f;
    um2 = .55797f;
    if (abs(*wat - um1) < 0.0001) {  *wat = um1 + 0.01; }; //to protect

    exbm1 = (um1 - um2) / (*wat - um1);

    if ((1.f - exbm1) ==0.0) {  exbm1 = 0.01;}; // to protect
    a = (*wat - um1) / (1.f - exbm1);

    if (abs(1.f - exbm1) < 0.0001) { exbm1 = 0.01;}; // to protect
    b = -log(exbm1);
    c__ = *wat - a;
    unf = a * exp(b * tempsoil) + c__;
    du = b * (unf - c__);

/* --- HEAT CAPACITY PER UNIT MASS */
L255:
    *sph = (*wat + unf) * 2093.f + 1925.6f + basicdata1.hlat * du;
    rs = 1550.f;
    rl = 1e3f;
    ri = 917.f;
    ra = 1.25f;
    xs = *dense / rs;
    xl = *dense * unf / rl;
    xi = *dense * (*wat - unf) / ri;
    xw = xl + xi;
    xp = xs + xw;
    xa = 1.f - xp;
    xwsat = 1.f - xs;
    if (xw <= xwsat) {
	     goto L256;
    }
    return 0;

L256:
    cs = .25f;
    cl = .56f;
    ci = 2.1f;
    fs = 1.25f;
    fi = .53f;
    hs = xs * fs;
    hi = xi * fi;
    if (xa >= .5f) {
	     goto L257;
    }
    fl = 1.f;
    if (abs( cl+ca) < 0.0001) { cl = -ca + 0.01; }; // to protect

    fa = (cl * 5.f + ca) / ((cl + ca) * 3.f);
    ha = xa * fa;

    if (abs(xl + hs + hi + ha) < 0.0001) { xl = -(hs +hi + ha) + 0.01; } ; // to protect
    *cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / (xl + hs + hi + ha);
    return 0;
L257:

    if ( abs(xl + hs + hi) < 0.0001) {  xl = -( hs + hi) + 0.01; }; // to protect
    cp = (xl * cl + hs * cs + hi * ci) / (xl + hs + hi);

    if (abs(ca + cp) < 0.0001) { ca = -cp + 0.01; }; // to protect
    fp = (ca * 5.f + cp) / ((ca + cp) * 3.f);

    hp = xp * fp;
    if (abs(xa + hp) < 0.0001) {  xa = -hp + 0.01; }; // to protect
    *cnd = (xa * ca + hp * cp) / (xa + hp);
    return 0;

L26:
    unf = *wat;
    du = 0.f;
    if (tempsoil >= basicdata1.tf) {
	     goto L265;
    }
    um1 = *wat * .15f;
    um4 = *wat * .01f;

    if (abs(*wat - um1) < 0.0001) { *wat = um1 + 0.01; }; //to protect
    exbm1 = (um1 - um4) / (*wat - um1);

    if (abs(1.f - exbm1) < 0.0001) { exbm1 = 0.01;}; // to protect
    a = (*wat - um1) / (1.f - exbm1);

    if (exbm1 <= 0.0) {  exbm1 = 0.01;}; // to protect
    b = -log(exbm1);
    c__ = *wat - a;
    sqt = sqrt(-tempsoil);
    unf = a * exp(-b * sqt) + c__;
	  if(sqt < 0.000001) {sqt = 0.001;}     
    du = b * .5f * (unf - c__) / sqt;
    if (unf - c__ > 1e-5f) {
	     goto L265;
    }
    unf = c__;
    du = 0.f;

/* --- APPARANT HEAT CAPACITY PER UNIT MASS */
L265:
    *sph = (*wat + unf) * 2093.f + 711.62f + du * 3.34e5f;
    rs = 2700.f;
    rl = 1e3f;
    ri = 917.f;
    xs = *dense / rs;
    xl = *dense * unf / rl;
    xi = *dense * (*wat - unf) / ri;
    xw = xl + xi;
    xwsat = 1.f - xs;
    xa = xwsat - xw;
    if (xw <= xwsat) {
	     goto L266;
    }

L266:
    cs = 7.f;
    cl = .56f;
    ci = 2.1f;
    ca = exp(tempsoil * .082f) * 5.96f / (tempsoil + 273.2f) + .0242f;
    ga = .1f;
    gc = .1f;
    cscl = cs / cl;

    double gagaga = 0.0;
    double gcgcgc = 0.0;
    
    	gagaga = (ga * (cscl - 1.f) + 1.f);
    	if (abs(gagaga) < 0.0001) { gagaga = 0.01;}  //to protect, zbl edit
    	gcgcgc = (gc * (cscl - 1.f) + 1.f);
    	if (abs(gcgcgc) < 0.0001) { gcgcgc = 0.01;}  //to protect, zbl edit
     
/*    if ((ga * (cscl - 1.f) + 1.f) ==0.0) { gaga = (ga * (cscl - 1.f) + 1.f) + 0.01;}
    else gagaga = (ga * (cscl - 1.f) + 1.f);
    if ((gc * (cscl - 1.f)+ 1.f) ==0.0) { gcgc = (gc * (cscl - 1.f)+ 1.f)+ 0.01;}
    else gcgcgc = (gc * (cscl - 1.f) + 1.f); */

//    fs = .666667f / (ga * (cscl - 1.f) + 1.f) + .333333f / (gc * (cscl - 1.f)    + 1.f);
    fs = .666667f / gagaga + .333333f / gcgcgc;

    cicl = ci / cl;

    double gagacic1 = 0.0;
    double gcgccic1 = 0.0;

    if (abs(ga * (cicl - 1.f) + 1.f) < 0.0001) { gagacic1 = (ga * (cicl - 1.f) + 1.f) + 0.01;}
    else gagacic1 = (ga * (cicl - 1.f) + 1.f);

    if (abs(gc * (cicl - 1.f)    + 1.f) < 0.0001) { gcgccic1 = (gc * (cicl - 1.f)    + 1.f) + 0.01; }
    else gcgccic1 = (gc * (cicl - 1.f)    + 1.f);

//    fi = .666667f / (ga * (cicl - 1.f) + 1.f) + .333333f / (gc * (cicl - 1.f)	    + 1.f);
     fi = .666667f / gagacic1 + .333333f / gcgccic1;

    cacl = ca / cl;

    double gacac1 = 0.0;
    double gccac1 = 0.0;

    if (abs(ga * (cacl - 1.f) + 1.f) < 0.0001) { gacac1 = (ga * (cacl - 1.f) + 1.f) + 0.01;}
    else gacac1 = (ga * (cacl - 1.f) + 1.f);

    if (abs(gc * (cacl - 1.f)    + 1.f) < 0.0001) { gccac1 = (gc * (cacl - 1.f)  + 1.f) + 0.01; }
    else gccac1 = (gc * (cacl - 1.f)    + 1.f);

//    fa = .666667f / (ga * (cacl - 1.f) + 1.f) + .333333f / (gc * (cacl - 1.f)	    + 1.f);
    fa = .666667f / gacac1 + .333333f / gccac1;

    hs = xs * fs;
    hi = xi * fi;
    ha = xa * fa;
    if (abs(xl + hs + hi + ha) < 0.0001) { xl = -(hs +hi + ha) + 0.01; } ; // to protect

    *cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / (xl + hs + hi + ha);
    return 0;


L27:
    unf = *wat;
    du = 0.f;
    if (tempsoil >= basicdata1.tf) {
	     goto L275;
    }
    um1 = *wat * .15f;
    um4 = *wat * .01f;

    if (abs(*wat - um1) < 0.0001) {  *wat = um1 + 0.01; }; //to protect
    exbm1 = (um1 - um4) / (*wat - um1);

    if (abs(1.f - exbm1) < 0.0001) {  exbm1 = 0.01;}; // to protect
    a = (*wat - um1) / (1.f - exbm1);

    if (exbm1 <= 0.0) {  exbm1 = 0.01;}; // to protect
    b = -log(exbm1);
    c__ = *wat - a;
    sqt = sqrt(-tempsoil);
    unf = a * exp(-b * sqt) + c__;
    if(sqt < 0.000001)  {sqt = 0.001;} // to protect, zbl edit
    du = b * .5f * (unf - c__) / sqt;
    if (unf - c__ > 1e-5f) {
	     goto L275;
    }
    unf = c__;
    du = 0.f;

 L275:
    *sph = (*wat + unf) * 2093.f + 711.62f + du * 3.34e5f;
    rs = 2700.f;
    rl = 1e3f;
    ri = 917.f;
    xs = *dense / rs;
    xl = *dense * unf / rl;
    xi = *dense * (*wat - unf) / ri;
    xw = xl + xi;
    xwsat = 1.f - xs;
    xa = xwsat - xw;
    if (xw <= xwsat) {
	     goto L276;
    }
L276:
    if (tempsoil < basicdata1.tf) {
	     d__1 = (doubledouble) (*dense * 8.12e-4f);
	     d__2 = (doubledouble) (*dense * 9.12e-4f);
	     *cnd = pow(c_b358, d__1) * .0109f + *wat * .46f * pow(c_b358, d__2);
    }
    if (tempsoil >= basicdata1.tf) {
       d__1 = (doubledouble) (*dense * 6.25e-4f);

       if (*wat <=0.0) *wat =0.01; // 28/04/2001

	     *cnd = (log(*wat) * .044f + .26f) * pow(c_b358, d__1);
    }
    return 0;

L28:
    goto L999;
L29:
    goto L999;
L30:
    goto L999;

L999:
    return 0;
}









// int Soilthermal::tinitl_(int& cmnt) //, double& tmax, double& tmin, double& tave, double& longmean, double newtemp[25], int& profilefla
int Soilthermal::tinitl_(int& cmnt, double& tmax, double& tmin, double& tave, double& longmean, double newtemp[25], int& profileflag)
{
    integer i__1;
    static double tempsoil = 0.0;
    static integer nmax = 0, jmax = 0;
    static double xbot = 0.0;
    static integer i__ = 0, j = 0, n = 0, index = 0;
    static double dtemp = 0.0, dxstr = 0.0, x0 = 0.0;
    static integer np = 0, nn = 0;
    static double tt[211], xx[211], vdepth = 0.0;
    extern  int intrpl_(double *, double *, integer *, integer *,
	    double *, double *, integer *);
    static integer nn1 = 0;
    static double tstart[211], xstart[211];
    static integer num = 0;
    double a = 0.0;
    double b = 0.0;
   
	for (i__1 = 0; i__1 < 211; i__1 ++)
	{
	tt[i__1] = 0.0;
	xx[i__1] = 0.0;
	tstart[i__1] = 0.0;
	xstart[i__1] = 0.0;	
	}    
    
  if (kswitch ==0 ) {
  // initfile >> dummy;
  // initfile >> dummy >> dummy>> dummy;
  // initfile>> index >>vdepth>>np;
      index = integer (basicf.INDEX[cmnt]);  
      vdepth = basicf.VDEPP[cmnt];
      np = basicf.NP[cmnt];
    // cout << " " << index <<endl;
    myindex=index;  myvdepth=vdepth; mynp=np;
    //printf("cmnt = %i, index = %i, depth = %.3f, np = %.3f\n", cmnt, index, vdepth, np);
    }
    else {
     index=myindex;  vdepth=myvdepth; np=mynp;
      }
      //printf("tsoil 2333: tmax = %.3f, tmin = %.3f, tave = %.3f\n", tmax, tmin, tave);
	double pi = 3.141592653f;
	double dm1 = 6.0 * acos(2.0 * tave / (tmax - tmin)) / pi; // month when average temp first exceed 0 degree
	double dm2 = 4.0 * pi - dm1; // month when average temp last exceed 0 degree
	// surface thaw index, defined as the accumulative surface temperature (C*day) above 0
	double sti = - 6.0 * (tmax - tmin) / 2.0 / pi * (sin(pi * dm2 / 6) - sin(pi * dm1 / 6)) + tave * (dm2 - dm1);
	sti = sti * 31;
	double stiroot;
		if (sti >= 0)
		{stiroot = sqrt(sti);}
		else {stiroot = 25.0;}
	// maximum active layer depth, porprotional to the square root of surface thaw index. 2.338814 comes from observation in Zackenberg site
	double aldmax = 2.338814f * stiroot / 100.0;
	//printf("2346: index = %i, profileflag = %i, aldmax = %.3f, dm1 = %.3f, dm2 = %.3f, sti = %.3f \n", index, profileflag, aldmax, dm1, dm2, sti);
   switch (index) {
	case 1:  goto L1000;
	case 2:  goto L2000;
	case 3:  goto L2000;
    }

L1000:
    return 0;
L2000:
//  if (kswitch ==0 ) {
//   initfile >> dummy >> dummy >>dummy;
//   initfile >> dummy>>dummy;
//   initfile >> dummy;
//  }
  n = 0;
/*<  100  N=N+1 >*/
L100:
    ++n;
     if (kswitch ==0 ) {
//     initfile >>xstart[n - 1]>>tstart[n - 1];
//     for (j=0;j<25; j++) {
//       xstart[n - 1] = basicf.DEPTH[cmnt][j];
//       tstart[n - 1] = basicf.TEMP[cmnt][j];
//      }
     // cout << tstart[n - 1]<<endl;
      xstart[n - 1] = basicf.DEPTH[cmnt][n-1];
      tstart[n - 1] = basicf.TEMP[cmnt][n-1];

     if(profileflag == 0){
		//tstart[n - 1] = longmean + (xstart[n - 1] - 10.0) * 20.0 / 1000.0;
		 if (longmean < 0){
			  if ((xstart[n - 1] < aldmax) && (tstart[n - 1] < 0.0))
        {
				     tstart[n - 1] = tstart[n - 2] / 1.5;
			  }
			  if ((xstart[n - 1] > aldmax) && (xstart[n - 2] < aldmax))
        {
				   tstart[n - 1] = -tstart[n - 2] / 1.5;
				//tstart[n - 1] = 0.5;
			  }
			  if ((xstart[n - 1] > aldmax) && (xstart[n - 2] > aldmax))
        {
				   a = longmean / (10.0 - aldmax);
				   b = - aldmax * longmean / (10.0 - aldmax);
				   tstart[n - 1] = xstart[n - 1] * a + b;
			  }
			  if (xstart[n - 1] > 10.0)
        {
				   tstart[n - 1] = (xstart[n - 1]  - 10.0) * 20.0 / 1000.0 + longmean;
			  }
		  }
		else  // if longmean >= 0
		{
		tstart[n - 1] = longmean + (xstart[n - 1] - 10.0) * 20.0 / 1000.0;
		}
		newtemp[n - 1] = tstart[n - 1]; 
		//printf("2403: n = %i, aldmax = %.3f, longmean = %.3f, depth = %.3f, temp = %.3f \n", n - 1, aldmax, longmean, xstart[n - 1], tstart[n - 1]);
		} // end if profileflag == 0
		else { // if profileflag != 0
			tstart[n - 1] = newtemp[n - 1]; //printf("2059: n = %i, aldmax = %.3f, depth = %.3f, temp = %.3f \n", n - 1, aldmax, xstart[n - 1], tstart[n - 1]);
			//if (n == 1)	
			//{printf("2066: n = %i, x = %.3f, t = %.3f \n", n - 1, xstart[n - 1], tstart[n - 1]);}
		}
		if (xstart[n - 1] < -0.0)
		{
			tstart[n - 1] = 0.0;
		} 
   
    myxstart[n - 1]=xstart[n - 1];
    mytstart[n - 1]=tstart[n - 1];
     }
  else {
       xstart[n - 1]=myxstart[n - 1];
       tstart[n - 1]=mytstart[n - 1];
     }
    //printf("2422: n = %i, aldmax = %.3f, depth = %.3f, temp = %.3f \n", n - 1, aldmax, xstart[n - 1], tstart[n - 1]);
    xstart[n - 1] = vdepth * xstart[n - 1];
    if (xstart[n - 1] > vdepth * -1.f) {
	     goto L100;
    }

    nmax = n - 1;

     if (nmax > 2) {
	      goto L3000;
    }

    if (abs(xstart[1] - xstart[0]) < 0.0001) {  xstart[1] = xstart[0] + 0.01;}; // to protect

    dtemp = (tstart[1] - tstart[0]) / (xstart[1] - xstart[0]);
    tempsoil = tstart[0];
    x0 = xstart[0];
    i__1 = basicdata1.imax;
    for (i__ = basicdata1.ig; i__ <= i__1; ++i__) {
	      basicdata1.t[i__ - 1] = tempsoil + dtemp * (basicdata1.x[i__ - 1] - x0);
    }
    return 0;

L3000:
    nn = n - 1;
    xbot = basicdata1.x[basicdata1.imax - 1];
    if (xstart[nn - 1] >= xbot) {
	     goto L106;
    }
    dxstr = xstart[nn - 1] - xstart[nn - 2];

    if (abs(dxstr) < 0.0001) {  dxstr =0.01;}; // to protect
    num = integer ((xbot - xstart[nn - 1]) / dxstr);
    if (num <= 1) {
	     num = 1;
    }
    dxstr = (xbot - xstart[nn - 1]) / num;
    dtemp = (tstart[nn] - tstart[nn - 1]) / num;
    nn1 = nn + 1;
    nmax = num + nn1;
    i__1 = nmax;
    for (n = nn1; n <= i__1; ++n) {
	      xstart[n - 1] = xstart[n - 2] + dxstr;
	      tstart[n - 1] = tstart[n - 2] + dtemp;
    }
L106:
    jmax = basicdata1.imax - basicdata1.igm1;
    i__1 = jmax;
    for (j = 1; j <= i__1; ++j) {
	      i__ = basicdata1.igm1 + j;
	      xx[j - 1] = basicdata1.x[i__ - 1];
    }
    intrpl_(xstart, tstart, &nmax, &np, xx, tt, &jmax);
    i__1 = jmax;
    for (j = 1; j <= i__1; ++j) {
	      i__ = basicdata1.igm1 + j;
	      basicdata1.t[i__ - 1] = tt[j - 1];
   }
    return 0;
} /* tinitl_ */



 int intrpl_(double *x, double *y, integer *nmax, integer *np,
	double *xx, double *yy, integer *nnmax)
{

    integer i__1;
    static integer m = 0, n = 0;
    static double slope = 0.0, x1 = 0, y1 = 0, x2 = 0;

    --yy;
    --xx;
    --y;
    --x;

    if (*nmax > 1500 || *nnmax > 3000) { return 0;
    }
    //xx[*nnmax + 1] = 1e10f;
    xx[*nnmax + 1] = xx[*nnmax];
    m = 1;
    i__1 = *nmax;
    for (n = 2; n <= i__1; ++n) {
	      y1 = y[n - 1];
	      x1 = x[n - 1];
	      x2 = x[n];
	      if (abs(x2 - x1) < 0.0001)  {x2 = x1 + 0.001; } // added

	  slope = (y[n] - y1) / (x2 - x1);
     
L301:
	yy[m] = y1 + slope * (xx[m] - x1);
	if (xx[m] > x2) {
	    goto L300;
	}
	++m;
	if (m <= *nnmax) {
	    goto L301;
	}
L300:
	;
    }
    return 0;
} /* intrpl_ */

 int Soilthermal::snofal_(double *snow, double *sden, double *comp, double *eta0,
	double *denfac, double *fcmelt, int& cmnt)
{

    integer i__1 = 0, i__2 = 0;

    static double abel = 0.0, dend = 0.0, elap = 0.0, dnew[365], rate = 0.0;
    static integer nmax = 0, mdat1 = 0, mdat2 = 0, m = 0, n = 0;
    static double arate[365], elraw[365];
    static integer index = 0;
    static double daton = 0.;
    static integer mhmax = 0, m1 = 0, m2 = 0;
    static double denmax = 0.0, datmax = 0.0, depmax = 0.0, snoden = 0.0, cummax = 0.0, epssno = 0.0, sphsno = 0.0,
	    convrt = 0.0, cnd = 0.0, day = 0.0;
    for (i__1 = 0; i__1 < 365; i__1 ++) {
	      arate[i__1] = 0.0;
	      elraw[i__1] = 0.0;
	      dnew[i__1] = 0.0;	
	  }
    --sden;
    --snow;
     if (kswitch ==0) {
//      snowdepth >> dummy >> dummy>> dummy>> dummy;
//      snowdepth >> index >> dummy;
      index=basicf.SNOFAL[cmnt];

      myindex2 = index;
     }
     else { index = myindex2; }

    if (index <= 0) {
	     return 0;
    }

    if (kswitch ==0) {
 //      snowdepth >> dummy >> dummy>>dummy >> dummy>> dummy>> dummy;
 //      snowdepth >>epssno>>convrt>>*eta0>>*denfac>>*fcmelt>>denmax;
         epssno= basicf.EPSSNO[cmnt];  // problem ??
         convrt=basicf.CONVRT[cmnt];
         *eta0= basicf.ETAO[cmnt];
         *denfac= basicf.DENFAC[cmnt];
         *fcmelt= basicf.FCMELT[cmnt];
         denmax= basicf.DENMAX[cmnt];

       // cout << *eta0 <<endl;
       myepssno=epssno; myconvrt=convrt; myeta0=*eta0; mydenfac =*denfac;
       myfcmelt = *fcmelt; mydenmax=denmax;
       // cout << *eta0 <<endl;
     }
     else {
       epssno=myepssno; convrt=myconvrt; *eta0=myeta0; *denfac =mydenfac;
       *fcmelt = myfcmelt; denmax=mydenmax;
     }
    abel = 2.9e-6f;
    cnd = abel * denmax * denmax;
    sphsno = 2e3f;
    *comp = sqrt(cnd * denmax * basicdata1.dt * epssno / sphsno);
    if (epssno < .1f) {
	     *comp = 1e10f;
    }
  switch (index) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
    }
L1:
    return 0;
L2:
    return 0;

L3:
    mhmax = basicdata1.mmax / 2 - 1;
    i__2 = mhmax;
    if (abs(basicdata1.mmax) < 0.0001) { 	basicdata1.mmax = 1; } // to protect

    for (m = 1; m <= i__2; ++m) {
	      snow[m] = (blk3_1.hsnow2 - blk3_1.hsnow1) * 500.f / (double)basicdata1.mmax;
	      sden[m] = 250.f;
    }
    ++mhmax;
    i__2 = basicdata1.mmax;
    for (m = mhmax; m <= i__2; ++m) {
	      snow[m] = (blk3_1.hsnow3 - blk3_1.hsnow2) * 500.f / (double) basicdata1.mmax;
	      sden[m] = 250.f;
    }
    return 0;
} /* snofal_ */


 int neige_(double *acum, double *snoden, double *weight, double *
	comp, double *eta0, double *denfac, double *fcmelt, double *tdrive)
{
    integer i__1;
    static double tsno = 0.0, wtis = 0.0;
    static integer i__ = 0, j = 0;
    static double ablat = 0.0, fonte = 0.0, dxold = 0.0, ratio = 0.0, wtisd = 0.0;
    static integer jj = 0;
    static double remain = 0.0, xx = 0.0, dsmass = 0.0, rit = 0.0;

    --weight;

    if (basicdata1.is > 2 && basicdata1.is <= basicdata1.ig) {
	     goto L10;
    }
L10:
    dsmass = *acum * basicdata1.dtday;
    if (basicdata1.smass <= 0.f && dsmass <= 0.f) {
	     return 0;
    }
    basicdata1.smass += dsmass;
    if (dsmass <= 0.f) {
	     goto L50;
    }
    tsno = basicdata1.t[basicdata1.is - 1];
    if (basicdata1.is == basicdata1.ig) {
	     basicdata1.is = basicdata1.igm1;
    }
    wtis = weight[basicdata1.is] + dsmass;
    if (wtis > *comp) {
	     goto L40;
    }
    if (abs(*snoden) < 0.0001) {  *snoden =0.01;}; // to protect

    if (abs(basicdata1.dx[basicdata1.is - 1] + dsmass / *snoden)) { basicdata1.dx[basicdata1.is - 1]= -(dsmass / *snoden)+ 0.01;}; // to protect

    basicdata1.ddry[basicdata1.is - 1] = wtis / (basicdata1.dx[basicdata1.is - 1] + dsmass / *snoden);
    if (weight[basicdata1.is] <= 0.f) {
	     basicdata1.xfa[basicdata1.is - 1] = -1e10f;
    }
    weight[basicdata1.is] = wtis;
    basicdata1.t[basicdata1.is - 1] = tsno;
    goto L60;

L40:
    wtisd = *comp - weight[basicdata1.is];

    if (abs(*snoden) < 0.0001) {   *snoden =0.01;}; // to protect
    if (abs(basicdata1.dx[basicdata1.is - 1] +  wtisd / *snoden) < 0.0001) { basicdata1.dx[basicdata1.is - 1] =  -( wtisd / *snoden) + 0.01;}; // to protect

    basicdata1.ddry[basicdata1.is - 1] = *comp / (basicdata1.dx[basicdata1.is - 1] +  wtisd / *snoden);
    weight[basicdata1.is] = *comp;
    remain = dsmass - wtisd;
    basicdata1.t[basicdata1.is - 1] = tsno;
L41:
    --basicdata1.is;
    weight[basicdata1.is] = *comp;
    basicdata1.ddry[basicdata1.is - 1] = *snoden;
    basicdata1.t[basicdata1.is - 1] = tsno;
    basicdata1.xfa[basicdata1.is - 1] = -1e10f;
    if (remain <= *comp) {
       goto L42;
    }
    remain -= *comp;
    goto L41;
L42:
    weight[basicdata1.is] = remain;
    goto L60;

L50:
    ablat = -dsmass;
L51:
    if (ablat < weight[basicdata1.is]) {
	     goto L52;
    }
    if (basicdata1.is == basicdata1.ig) {
	     goto L54;
    }
    ablat -= weight[basicdata1.is];
    weight[basicdata1.is] = 0.f;
    basicdata1.dx[basicdata1.is - 1] = 0.f;
    basicdata1.xfa[basicdata1.is - 1] = -1e10f;
    ++basicdata1.is;
    goto L51;
L52:
    weight[basicdata1.is] -= ablat;
    goto L60;
L54:
    basicdata1.smass = 0.f;
    weight[basicdata1.igm1] = 0.f;
    basicdata1.xfa[basicdata1.igm1 - 1] = -1e10f;
    return 0;
L60:
    if (*tdrive <= 0.f) {
	     goto L70;
    }
    fonte = *fcmelt * basicdata1.dtday * basicdata1.tair;
    basicdata1.smass -= fonte;
    if (basicdata1.smass > 0.f) {
	     *tdrive = .01f;
    }
L61:
    if (fonte < weight[basicdata1.is]) {
	     goto L62;
    }
    if (basicdata1.is == basicdata1.ig) {
	     goto L64;
    }
    fonte -= weight[basicdata1.is];
    weight[basicdata1.is] = 0.f;
    basicdata1.dx[basicdata1.is - 1] = 0.f;
    basicdata1.xfa[basicdata1.is - 1] = -1e10f;
    ++basicdata1.is;
    goto L61;
L62:
    weight[basicdata1.is] -= fonte;
    goto L70;
L64:
    basicdata1.smass = 0.f;
    weight[basicdata1.igm1] = 0.f;
    basicdata1.xfa[basicdata1.igm1 - 1] = -1e10f;
    return 0;
L70:
    i__1 = basicdata1.igm1;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {
	      basicdata1.mater[i__ - 1] = 21;
	      rit = basicdata1.ddry[i__ - 1];
	      dxold = basicdata1.dx[i__ - 1];

   if (abs(rit) < 0.0001) {  rit = 0.01;}; // to protect
	 basicdata1.dx[i__ - 1] = weight[i__] / rit;
	 if (basicdata1.xfa[i__ - 1] < 0.f) {
	    goto L71;
	 }
   if (abs(dxold) < 0.0001) {  dxold = 0.01;}; // to protect
	 ratio = basicdata1.dx[i__ - 1] / dxold;
	 if (basicdata1.xfb[i__ - 1] == 0.f) {
	    goto L888;
	 }
	 basicdata1.xfb[i__ - 1] += (ratio - 1.f) * (basicdata1.dx[i__ - 1] -
		basicdata1.xfa[i__ - 1]);
	 if (basicdata1.xfb[i__ - 1] < 0.f) {
	    basicdata1.xfb[i__ - 1] = basicdata1.dx[i__ - 1] * .01f;
	 }
   
L888:
	basicdata1.xfa[i__ - 1] = ratio * basicdata1.xfa[i__ - 1];
L71:
	;
    }
    xx = basicdata1.x[basicdata1.ig - 1];
    jj = basicdata1.is + basicdata1.igm1;
    i__1 = basicdata1.igm1;
    for (j = basicdata1.is; j <= i__1; ++j) {
	      i__ = jj - j;
	      xx -= basicdata1.dx[i__ - 1];
	      basicdata1.x[i__ - 1] = xx;
    }

    i__1 = basicdata1.ig;
    for (i__ = basicdata1.is; i__ <= i__1; ++i__) {
	      if (basicdata1.t[i__ - 1] > 0.f) {
	         basicdata1.t[i__ - 1] = .01f;
	      }
    }
    return 0;
} /* neige_ */

 int asmone_(integer *i__)
{
    double r__1;

    static double gold = 0.0, dens = 0.0;
    extern  int tabv_(integer *), pram_(integer *, double *,
	    double *, double *, double *, double *, double *), tblo_(integer *);
    static double g = 0.0;
    extern  int phase_(integer *, double *, double *, double *,
	    double *, double *, double *, double *, double *, double *, double *, double *,
	    double *);
    static double dtwas = 0.0, t1 = 0.0, t2 = 0.0, cd = 0.0, cu = 0.0, td = 0.0, tt = 0.0;
    extern  int asmabv_(integer *), asmblo_(integer *);
    static integer im1 = 0, ip1 = 0, ip2 = 0;
    static double cnd = 0.0, rsd = 0.0, sph = 0.0, wat = 0.0, dxx = 0.0, rlw = 0.0, rsu = 0.0;

    integer counter_asmone;

    dtwas = basicdata1.dtfaz;
    gold = basicdata1.xfa[*i__ - 1];
    im1 = *i__ - 1;
    ip1 = *i__ + 1;
    ip2 = *i__ + 2;
    counter_asmone =0;

L10:
    if (counter_asmone < 1e5)
    {

    asmabv_(&im1);
    asmblo_(&ip2);

    dxx = basicdata1.dx[*i__ - 1];
    wat = basicdata1.water[*i__ - 1];
    dens = basicdata1.ddry[*i__ - 1];
    tt = basicdata1.t[*i__ - 1];
    pram_(&basicdata1.mater[*i__ - 1], &tt, &basicdata1.tf, &basicdata1.ddry[*i__ -
	    1], &basicdata1.water[*i__ - 1], &sph, &cnd);

    cu = cnd;
    rsu = dens * sph;
    td = basicdata1.t[ip1 - 1];
    pram_(&basicdata1.mater[*i__ - 1], &basicdata1.tf, &td, &basicdata1.ddry[*i__ -
	    1], &basicdata1.water[*i__ - 1], &sph, &cnd);

    cd = cnd;
    rsd = dens * sph;
    rlw = dens * basicdata1.hlat * wat;

    phase_(i__, &rlw, &cu, &rsu, &cd, &rsd, &tt, &td, &dxx, &gold, &g, &t1, &
	    t2);

    counter_asmone ++;

   }

    if (g >= 0.f) {
	     goto L12;
    }
    basicdata1.capx[*i__ - 1] = rsd * dxx;

    if (abs(dxx) < 0.01) { dxx =0.01;}; // to protect
    basicdata1.conx[*i__ - 1] = cd / dxx;
    basicdata1.xfa[*i__ - 1] = -1e10f;
    if (*i__ == basicdata1.is) {
	     return 0;
    }
    --(*i__);
    im1 = *i__ - 1;
    ip1 = *i__ + 1;
    ip2 = *i__ + 2;
    r__1 = tt - basicdata1.tf;
//    basicdata1.t[ip1 - 1] = basicdata1.tf - r_sign(&c_b127, &r__1);
  if (r__1 >=0 )  basicdata1.t[ip1 - 1] = basicdata1.tf - c_b127;
  else basicdata1.t[ip1 - 1] = basicdata1.tf + c_b127;

    asmabv_(i__);
    tabv_(i__);
    asmblo_(&ip2);
    tblo_(&ip2);
    basicdata1.dtfaz = dtwas - basicdata1.dtfaz;
    dtwas = basicdata1.dtfaz;
    gold = basicdata1.dx[*i__ - 1];
    goto L10;

L12:
    //zbl added to prevent infinite loop
        //printf("basicdata.arr[0] = %f \n", basicdata.arr[0]);  //zbl
	if (basicdata1.arr[0] == -99.0) {basicdata1.arr[0] = g; /*printf("***01\n"); */}
     else {int i;
  	      for (i = 1; i <= 3000; i++){
	 	         if (basicdata1.arr[i]== g)  { g = g - 0.001; basicdata1.arr[i + 1]== g; /*printf("***02\n");*/}
		         else {
		             if (basicdata1.arr[i] == -99.0) {basicdata1.arr[i] = g; /*printf("***03\n");*/break; }	; 				
		       }
     }
   }
   //end of add
   
    if (g > dxx) {
	     goto L13;
    }
    basicdata1.t[*i__ - 1] = t1;
    tabv_(&im1);
    basicdata1.t[ip1 - 1] = t2;
    tblo_(&ip2);
    basicdata1.xfa[*i__ - 1] = g;
    return 0;

L13:

   if ( abs(dxx) < 0.01) {dxx =0.01; }; // to protect

    basicdata1.capx[*i__ - 1] = rsu * dxx;
    basicdata1.conx[*i__ - 1] = cu / dxx;
    basicdata1.xfa[*i__ - 1] = -1e10f;
    if (*i__ == basicdata1.imax1) {
	     return 0;
    }
    ++(*i__);
    im1 = *i__ - 1;
    ip1 = *i__ + 1;
    ip2 = *i__ + 2;
    r__1 = td - basicdata1.tf;
//    basicdata1.t[*i__ - 1] = basicdata1.tf - r_sign(&c_b127, &r__1);
    if (r__1 >=0) basicdata1.t[*i__ - 1] = basicdata1.tf - c_b127;
    else  basicdata1.t[*i__ - 1] = basicdata1.tf +c_b127;

    asmabv_(&im1);
    tabv_(&im1);
    asmblo_(&ip1);
    tblo_(&ip1);
    basicdata1.dtfaz = dtwas - basicdata1.dtfaz;
    dtwas = basicdata1.dtfaz;
    gold = 0.f;
    goto L10;
} /* asmone_ */

 int asmtwo_(integer *itzero, integer *ibzero)
{
    double r__1 = 0.0;

    static double gold = 0.0, dens = 0.0;
    extern  int tabv_(integer *), pram_(integer *, double *,
	    double *, double *, double *, double *, double *), tblo_(integer *);
    static integer i__ = 0;
    static double gaold = 0.0, gbold = 0.0, v = 0.0;
    extern  int phase_(integer *, double *, double *, double *,
	    double *, double *, double *, double *, double *, double *, double *, double *,
	    double *);
    static double dtwas = 0.0, t1 = 0.0, t2 = 0.0, cd = 0.0, ga = 0.0, gb = 0.0, td = 0.0, cu = 0.0, gv = 0.0, tt = 0.0;
    extern  int asmabv_(integer *), asmblo_(integer *);
    static integer im1 = 0, ip1 = 0, ip2 = 0;
    static double cnd = 0.0, gdn = 0.0, vdn = 0.0, rsd = 0.0, sph = 0.0, wat = 0.0, gup = 0.0, dxx = 0.0, rlw = 0.0, rsu = 0.0, vup = 0.0;

    gaold = basicdata1.xfa[*itzero - 1];
    gbold = basicdata1.xfa[*ibzero - 1] + basicdata1.xfb[*ibzero - 1];
    i__ = *itzero;
    dtwas = basicdata1.dtfaz;
    v = basicdata1.xfb[i__ - 1];
    im1 = i__ - 1;
    ip1 = i__ + 1;
    ip2 = i__ + 2;
    
L10:

    asmabv_(&im1);

    dxx = basicdata1.dx[i__ - 1];
    wat = basicdata1.water[i__ - 1];
    dens = basicdata1.ddry[i__ - 1];
    tt = basicdata1.t[i__ - 1];
    td = basicdata1.t[ip1 - 1];
    pram_(&basicdata1.mater[i__ - 1], &tt, &basicdata1.tf, &basicdata1.ddry[i__ - 1]
	    , &basicdata1.water[i__ - 1], &sph, &cnd);
    cu = cnd;
    rsu = dens * sph;
    rlw = dens * basicdata1.hlat * wat;

    phase_(&i__, &rlw, &cu, &rsu, &c_b644, &c_b644, &tt, &basicdata1.tf, &dxx, &
	    gaold, &ga, &t1, &t2);

/* --- EVALUATE WHETHER GA GOES OUT OF ELEMENT */
    if (v != 0.f) {
	     goto L11;
    }
    if (ga > dxx) {
	     goto L12;
    }
    basicdata1.t[i__ - 1] = t1;
    tabv_(&im1);
    basicdata1.xfa[i__ - 1] = ga;
    goto L200;
/* --- SECOND PHASE PLANE WITHIN ELEMENT EXISTS */
L11:
    basicdata1.t[i__ - 1] = t1;
    tabv_(&im1);
    gv = v + gaold;
    basicdata1.xfb[i__ - 1] = gv - ga;
    basicdata1.xfa[i__ - 1] = ga;
    if (gv > ga) {
	     goto L200;
    }
    basicdata1.xfa[i__ - 1] = -1e10f;
    basicdata1.xfb[i__ - 1] = 0.f;
    if (i__ < *ibzero) {
	      goto L200;
    }
    basicdata1.xfa[i__ - 1] = ga;
    basicdata1.xfb[i__ - 1] = gv - ga;
    goto L200;

L12:
/* --- OUT BELOW */

    if (abs(dxx) < 0.01) { dxx =0.01; } ; // to protect
    gdn = basicdata1.xfa[ip1 - 1];
    basicdata1.capx[i__ - 1] = rsu * dxx;
    basicdata1.conx[i__ - 1] = cu / dxx;
    basicdata1.xfa[i__ - 1] = -1e10f;
    if (i__ == basicdata1.imax1) {
	     return 0;
    }
    ++i__;
    im1 = i__ - 1;
    ip1 = i__ + 1;
    ip2 = i__ + 2;
    basicdata1.dtfaz = dtwas - basicdata1.dtfaz;
    dtwas = basicdata1.dtfaz;
    r__1 = td - basicdata1.tf;
//    basicdata1.t[i__ - 1] = basicdata1.tf - r_sign(&c_b127, &r__1);
    if (r__1 >=0)  basicdata1.t[i__ - 1] = basicdata1.tf - c_b127;
    else basicdata1.t[i__ - 1] = basicdata1.tf + c_b127;

    asmabv_(&im1);
    tabv_(&im1);
    gaold = 0.f;
    if (gdn == -1e10f) {
	     goto L10;
    }
    vdn = basicdata1.xfb[i__ - 1];
    if (vdn != 0.f) {
	     goto L13;
    }
    v = gdn;
    goto L10;
L13:
    gaold = vdn + gdn;
    v = 0.f;
    basicdata1.xfb[i__ - 1] = 0.f;
    goto L10;

L200:
    if (*ibzero < *itzero) {
	     return 0;
    }

/* --- SECTION FOR LOWER PHASE PLANE */
    i__ = *ibzero;
    dtwas = basicdata1.dt;
    basicdata1.dtfaz = basicdata1.dt;
    gold = basicdata1.xfa[i__ - 1];
    v = basicdata1.xfb[i__ - 1];
    im1 = i__ - 1;
    ip1 = i__ + 1;
    ip2 = i__ + 2;

//	int phase_c =0; // changed
L20:

    asmblo_(&ip2);
    dxx = basicdata1.dx[i__ - 1];
    wat = basicdata1.water[i__ - 1];
    dens = basicdata1.ddry[i__ - 1];
    tt = basicdata1.t[i__ - 1];
    td = basicdata1.t[ip1 - 1];
    pram_(&basicdata1.mater[i__ - 1], &td, &basicdata1.tf, &basicdata1.ddry[i__ - 1]
	    , &basicdata1.water[i__ - 1], &sph, &cnd);


    cd = cnd;
    rsd = dens * sph;
    rlw = dens * basicdata1.hlat * wat;

//   if (phase_c < 100)

	phase_(&i__, &rlw, &c_b644, &c_b644, &cd, &rsd, &basicdata1.tf, &td, &dxx, &
	    gbold, &gb, &t1, &t2);


//  phase_c++; }
//   else return 0;
/* --- EVALUATE WHETHER G GOES OUT OF ELEMENT */

	if (v != 0.f) {
	   goto L21;
  }
  if (gb < 0.f) {
	   goto L22;
  }
    basicdata1.t[ip1 - 1] = t2;


    tblo_(&ip2);

    basicdata1.xfa[i__ - 1] = gb;

    return 0;


L21:
    basicdata1.t[ip1 - 1] = t2;

    tblo_(&ip2);

    basicdata1.xfa[i__ - 1] = gold;
    basicdata1.xfb[i__ - 1] = gb - gold;
    if (gb > gold) {
 	     return 0;
    }

    basicdata1.xfa[i__ - 1] = -1e10f;
    basicdata1.xfb[i__ - 1] = 0.f;
    return 0;
L22:
/* --- OUT ABOVE */
    gup = basicdata1.xfa[im1 - 1];
    basicdata1.capx[i__ - 1] = rsd * dxx;

    if (abs(dxx) < 0.01){   dxx =0.01;}; // to protect
    basicdata1.conx[i__ - 1] = cd / dxx;
    basicdata1.xfa[i__ - 1] = -1e10f;
    if (i__ == basicdata1.is) {
	     return 0;
    }
    --i__;
    im1 = i__ - 1;
    ip1 = i__ + 1;
    ip2 = i__ + 2;
    basicdata1.dtfaz = dtwas - basicdata1.dtfaz;
    dtwas = basicdata1.dtfaz;
    r__1 = tt - basicdata1.tf;
//    basicdata1.t[ip1 - 1] = basicdata1.tf - r_sign(&c_b127, &r__1);
    if (r__1 >=0 )  basicdata1.t[ip1 - 1] = basicdata1.tf - c_b127;
    else  basicdata1.t[ip1 - 1] = basicdata1.tf + c_b127;

    asmblo_(&ip2);
    tblo_(&ip2);
    gbold = basicdata1.dx[i__ - 1];
    v = 0.f;
    if (gup == -1e10f) {
       goto L20;
    }
    gold = gup;
    vup = basicdata1.xfb[i__ - 1];
    if (vup != 0.f) {
	     goto L23;
    }
    v = gbold - gup;
    goto L20;
L23:
    gbold = gup;
    v = 0.f;
    basicdata1.xfb[i__ - 1] = 0.f;
    goto L20;


} /* asmtwo_ */

 int asmblo_(integer *i1)
{
    integer i__1 = 0;

    static double denm = 0.0, cond = 0.0, conu = 0.0;
    static integer i__ = 0, j = 0;
    static double r__ = 0.0;
    static integer i2 = 0;
    static double condt1 = 0., conut1 = 0.;
    static integer jj = 0;
    static double rc = 0.0, condth = 0.0, conuth = 0.0;
    static integer im1 = 0, ip1 = 0;
    static double con = 0.0, rhs = 0.0;

    con = basicdata1.conx[basicdata1.imax1 - 1];

    if (abs(basicdata1.dt) < 0.0001) { basicdata1.dt =0.01; }; // to protect

    rc = basicdata1.capx[basicdata1.imax1 - 1] * .5f / basicdata1.dt;
    denm = rc + con + basicdata1.hbot;
    if (abs(denm) < 0.0001) {   denm =30.0;}; // to protected   debugged Q. Z. 12/Dec/2000;

    basicdata1.s[basicdata1.imax - 1] = con / denm;
    if (abs(basicdata1.t[basicdata1.imax - 1]) >= abs(denm)) {basicdata1.t[basicdata1.imax - 1] = -3.0;}  //zbl added
    basicdata1.e[basicdata1.imax - 1] = (rc * basicdata1.t[basicdata1.imax - 1] -
	    basicdata1.gflux + basicdata1.hbot * basicdata1.tbot + basicdata1.ht[
	    basicdata1.imax - 1]) / denm;
    if (*i1 >= basicdata1.imax) {
	     return 0;
    }

    i2 = *i1;
    if (*i1 == basicdata1.is) {
	     i2 = basicdata1.is + 1;
    }
    jj = i2 + basicdata1.imax1;
    i__1 = basicdata1.imax1;
    for (j = i2; j <= i__1; ++j) {
	      i__ = jj - j;
	      im1 = i__ - 1;
	      ip1 = i__ + 1;

   if ( abs(basicdata1.dt) < 0.00010) { basicdata1.dt =0.01; }; // to protect
	 rc = (basicdata1.capx[i__ - 1] + basicdata1.capx[im1 - 1]) / basicdata1.dt;
	 conut1 = basicdata1.conx[im1 - 1] * basicdata1.theta1;
	 condt1 = basicdata1.conx[i__ - 1] * basicdata1.theta1;
	 rhs = (rc - conut1 - condt1) * basicdata1.t[i__ - 1] + conut1 *
		basicdata1.t[im1 - 1] + condt1 * basicdata1.t[ip1 - 1];
	 rhs = rhs + basicdata1.theta * basicdata1.ht[i__ - 1] + basicdata1.theta1 *
		basicdata1.htold[i__ - 1];
	 conuth = basicdata1.conx[im1 - 1] * basicdata1.theta;
	 condth = basicdata1.conx[i__ - 1] * basicdata1.theta;
	 r__ = conuth + rc + condth;
	 denm = r__ - condth * basicdata1.s[ip1 - 1];

   if (abs(denm) < 0.0001)  { denm =30.0;}; // to protected   debugged Q. Z. 12/Dec/2000;
	 basicdata1.s[i__ - 1] = conuth / denm;
	 basicdata1.e[i__ - 1] = (rhs + condth * basicdata1.e[ip1 - 1]) / denm;
   if(basicdata1.e[i__ - 1] > 50.0) {basicdata1.e[i__ - 1] = 50.0;}
	 if(basicdata1.e[i__ - 1] < -50.0) {basicdata1.e[i__ - 1] = -50.0;}
    }
    if (*i1 > basicdata1.is) {
	     return 0;
    }

    if (abs(basicdata1.dt) < 0.0001) { basicdata1.dt =0.01; }; // to protect
    rc = basicdata1.capx[basicdata1.is - 1] * .5f / basicdata1.dt;
    conu = basicdata1.conx[basicdata1.ism1 - 1];
    cond = basicdata1.conx[basicdata1.is - 1];
    rhs = rc * basicdata1.t[basicdata1.is - 1] + basicdata1.sflux + basicdata1.ht[i__
	    - 1];
    r__ = conu + rc + cond;
    denm = r__ - cond * basicdata1.s[i2 - 1];

    if (abs(denm) < 0.00010) {denm =30.0;}; // to protected   debugged Q. Z. 12/Dec/2000;
    basicdata1.s[basicdata1.is - 1] = conu / denm;
    basicdata1.e[basicdata1.is - 1] = (rhs + cond * basicdata1.e[i2 - 1]) / denm;
    if(basicdata1.e[i__ - 1] > 50.0) {basicdata1.e[i__ - 1] = 50.0;}
		if(basicdata1.e[i__ - 1] < -50.0) {basicdata1.e[i__ - 1] = -50.0;}
    return 0;
} /* asmblo_ */

 int tblo_(integer *i1)
{
    integer i__1 = 0;

    static integer i__ = 0;

    if (*i1 > basicdata1.imax) {
	     return 0;
    }
    i__1 = basicdata1.imax;
    for (i__ = *i1; i__ <= i__1; ++i__) {
	      basicdata1.t[i__ - 1] = basicdata1.s[i__ - 1] * basicdata1.t[i__ - 2] +
		basicdata1.e[i__ - 1];
    }
    return 0;
} /* tblo_ */

 int asmabv_(integer *i1)
{
    integer i__1 = 0;

    static double cond = 0.0, denm = 0.0, conu = 0.0;
    static integer i__ = 0;
    static double r__ = 0.0, condt1 = 0.0, conut1 = 0.0, rc = 0.0, condth = 0.0, conuth = 0.0;
    static integer im1 = 0, ip1 = 0;
    static double rhs = 0.0;
    static integer isp1 = 0;

    if (*i1 < basicdata1.is) {
	     return 0;
    }

   if ( abs(basicdata1.dt) < 0.0001) { basicdata1.dt =0.01; }; // to protect

    rc = basicdata1.capx[basicdata1.is - 1] * .5f / basicdata1.dt;

    conu = basicdata1.conx[basicdata1.ism1 - 1];
    cond = basicdata1.conx[basicdata1.is - 1];
    rhs = rc * basicdata1.t[basicdata1.is - 1] + basicdata1.sflux + basicdata1.ht[
	    basicdata1.is - 1];
    r__ = conu + rc + cond;
    denm = r__ - conu * basicdata1.s[basicdata1.ism1 - 1];

    if (abs(denm) < 0.0001)  {  denm =30.0;}; // to protected   debugged Q. Z. 12/Dec/2000;
    basicdata1.s[basicdata1.is - 1] = cond / denm;

    basicdata1.e[basicdata1.is - 1] = (rhs + conu * basicdata1.e[basicdata1.ism1 - 1])
	     / denm;
    if (*i1 == basicdata1.is) {
	     return 0;
    }

    isp1 = basicdata1.is + 1;
    i__1 = *i1;

   for (i__ = isp1; i__ <= i__1; ++i__) {
  	   im1 = i__ - 1;
	     ip1 = i__ + 1;

   if ( basicdata1.dt ==0.0) { basicdata1.dt =0.01; }; // to protect
	rc = (basicdata1.capx[i__ - 1] + basicdata1.capx[im1 - 1]) / basicdata1.dt;
	conu = basicdata1.conx[im1 - 1];
	cond = basicdata1.conx[i__ - 1];
	conut1 = basicdata1.theta1 * conu;
	condt1 = basicdata1.theta1 * cond;
	rhs = (rc - conut1 - condt1) * basicdata1.t[i__ - 1] + conut1 *
		basicdata1.t[im1 - 1] + condt1 * basicdata1.t[ip1 - 1];
	rhs = rhs + basicdata1.theta * basicdata1.ht[i__ - 1] + basicdata1.theta1 *
		basicdata1.htold[i__ - 1];
	conuth = basicdata1.theta * conu;
	condth = basicdata1.theta * cond;
	r__ = conuth + rc + condth;
	denm = r__ - conuth * basicdata1.s[im1 - 1];

   if (abs(denm) < 0.0001)  {  denm =30.0;}; // to protected   debugged Q. Z. 12/Dec/2000;
  	basicdata1.s[i__ - 1] = condth / denm;
	  basicdata1.e[i__ - 1] = (rhs + conuth * basicdata1.e[im1 - 1]) / denm;
		if(basicdata1.e[i__ - 1] > 50.0) {basicdata1.e[i__ - 1] = 50.0;}
		if(basicdata1.e[i__ - 1] < -50.0) {basicdata1.e[i__ - 1] = -50.0;}
    }

    return 0;
} /* asmabv_ */

 int tabv_(integer *i1)
{
    integer i__1 = 0;

    static integer i__ = 0, j = 0, jj = 0;

    if (*i1 < basicdata1.is) {
	     return 0;
    }
    jj = basicdata1.is + *i1;
    i__1 = *i1;
    for (j = basicdata1.is; j <= i__1; ++j) {
	      i__ = jj - j;
	      basicdata1.t[i__ - 1] = basicdata1.s[i__ - 1] * basicdata1.t[i__] + basicdata1.e[i__ - 1];
    }
    return 0;
} /* tabv_ */




 int phase_(integer *i__, double *rlw, double *cu, double *rsu,
	double *cd, double *rsd, double *tt, double *td, double *dxx, double *gold, double *
	g, double *t1, double *t2)
{
    double r__1 = 0.0;

    static double capd = 0.0, cdtf = 0.0, cond = 0.0, hold = 0.0, capu = 0.0, conu = 0.0, cutf = 0.0, star = 0.0, fwas = 0.0, gwas = 0.0,
	    gnew = 0.0, rlwt = 0.0, f = 0.0, h__ = 0.0;
    static integer n = 0;
    /*static double b1, c1, d1, a1, d2, c2, b2, a2, df, dg, ft, a1p, b1p, a2p,
	    b2p;*/
    static double  b1 = 0.0, c1 = 0.0, d1 = 0.0, a1 = 0.0, d2 = 0.0, c2 = 0.0, b2 = 0.0, a2 = 0.0, df = 0.0, dg = 0.0, ft = 0.0, a1p = 0.0, b1p = 0.0, a2p = 0.0, b2p = 0.0;
    static integer im1 = 0, ip1 = 0, ip2 = 0;
    //static double fld, flu, err, alf1, alf2, gam1, gam2, bet1, bet2, rsd4, rsu4;
    static double fld = 0.0, flu = 0.0, err = 0.0, alf1 = 0.0, alf2 = 0.0, gam1 = 0.0, gam2 = 0.0, bet1 = 0.0, bet2 = 0.0, rsd4 = 0.0, rsu4 = 0.0;
    star = *dxx * .05f;
    hold = *dxx - *gold;
    im1 = *i__ - 1;
    ip1 = *i__ + 1;
    ip2 = *i__ + 2;


    if (abs(basicdata1.dtfaz) < 1.0)  {  basicdata1.dtfaz = basicdata1.dtday * 86400.f;}; // debugging by Q.Z. Nov/28/2000

    conu = basicdata1.conx[im1 - 1];
    capu = basicdata1.capx[im1 - 1] * .5f / basicdata1.dtfaz;
    cond = basicdata1.conx[ip1 - 1];
    capd = basicdata1.capx[ip1 - 1] * .5f / basicdata1.dtfaz;
    rsu4 = *rsu * .25f;
    d1 = rsu4 / basicdata1.dtfaz;
    c1 = conu * (1.f - basicdata1.s[im1 - 1]) + capu + d1 * *gold;
    b1 = d1 * *tt;
    a1 = (capu + d1 * *gold) * *tt + conu * 1.00001f * basicdata1.e[im1 - 1] +
	    basicdata1.ht[*i__ - 1];


    if (*i__ == basicdata1.is) {
	     a1 += basicdata1.sflux;
    }
    cutf = *cu * basicdata1.tf;
    a1p = a1 - basicdata1.tf * c1;
    b1p = b1 - basicdata1.tf * d1;
    rsd4 = *rsd * .25f;
    d2 = rsd4 / basicdata1.dtfaz;
    c2 = cond * (1.f - basicdata1.s[ip2 - 1]) + capd + d2 * hold;
    b2 = d2 * *td;
    a2 = (capd + d2 * hold) * *td + cond * basicdata1.e[ip2 - 1] + basicdata1.ht[ip1 - 1];
    cdtf = *cd * basicdata1.tf;
    a2p = a2 - basicdata1.tf * c2;
    b2p = b2 - basicdata1.tf * d2;

    if ((r__1 = a1p * *cu) < 0.f) {
	     goto L1;
    } else if ((r__1= a1p * *cu) == 0) {
	     goto L2;
    } else { if ((r__1= a1p * *cu) > 0)
	     goto L5;
    }
L1:
    *rlw = -(*rlw);
    goto L5;
L2:
    if ((r__1 = a2p * *cd) < 0.f) {
	     goto L5;
    } else if ((r__1 = a2p * *cd) == 0) {
	     goto L3;
    } else { if ((r__1 = a2p * *cd) > 0)
	     goto L4;
    }
L3:
    return 0;
L4:
    *rlw = -(*rlw);
L5:

    if (abs(basicdata1.dtfaz) < 1.0) { basicdata1.dtfaz = basicdata1.dtday * 86400.f;}; // debugging by Q.Z. Nov/28/2000
    rlwt = *rlw / basicdata1.dtfaz;
    alf1 = *cu * (a1p + b1p * *gold);
    gam1 = c1 + d1 * *gold;
    bet1 = *cu + gam1 * *gold;

  //  printf("%5.3f ", *cd); // ??? problem 4/29/2001

//    printf("ok1");

    alf2 = *cd * (a2p + b2p * hold);

 //   printf("ok2, %4.2f", alf2);

    gam2 = c2 + d2 * hold;
    bet2 = *cd + gam2 * hold;

    double dgdgdg= 0.0;                            // to protect
    if ( abs(rlwt * bet1 * bet2 + alf1 * gam2 -alf2 * gam1) < 0.0001)  {  dgdgdg =0.01;}
    else dgdgdg = rlwt * bet1 * bet2 + alf1 * gam2 -alf2 * gam1;
//    dg = (alf2 * bet1 + alf1 * bet2) / (rlwt * bet1 * bet2 + alf1 * gam2 -   alf2 * gam1);
    dg = (alf2 * bet1 + alf1 * bet2) / dgdgdg;


    if (dg >= *dxx) {
	     dg = *dxx * .99f;
    }
    if (dg <= 0.f) {
	     dg = *dxx * .1f;
    }
    *g = *gold + dg * .5f;
    h__ = *dxx - *g;
    ft = rlwt * (*g - *gold);

    if ( abs(*cu + (c1 + d1 * *g) * *g) < 0.01) {   *cu = (c1 + d1 * *g) * *g + 0.01;}; // to protect Q. Z. Dec/12/2000
    if (abs(*cd + (c2 + d2 * h__) * h__) < 0.01) { *cd = -(c2 + d2 * h__) * h__ + 0.01;}; // to protect

    flu = *cu * (a1p + b1p * *g) / (*g * (c1 + d1 * *g) + *cu);
    fld = *cd * (a2p + b2p * h__) / (h__ * (c2 + d2 * h__) + *cd);

    fwas = ft - flu - fld;
    gwas = *g;
    *g = *gold + dg;

    if (*rlw == 0.f) {
	     goto L30;
    }

  err = 1e-5f;

  for (n = 1; n <= 50; ++n) {

	  h__ = *dxx - *g;
	  ft = rlwt * (*g - *gold);

    // the following modification is made on 28/NOv/2000 QZ
   double ttt1;
   double ttt2;
   ttt1 = *g * (c1 + d1 * *g) + *cu;
   ttt2= h__ * (c2 + d2 * h__) + *cd ;
   if ( ttt1 ==0.0) ttt1 = 3.01;
   if (ttt2 ==0.0) ttt2  =3.001;

//	flu = *cu * (a1p + b1p * *g) / (*g * (c1 + d1 * *g) + *cu);
//   fld = *cd * (a2p + b2p * h__) / (h__ * (c2 + d2 * h__) + *cd);

	 flu = *cu * (a1p + b1p * *g) / ttt1;
   fld = *cd * (a2p + b2p * h__) / ttt2;

	 if (*gold <= star || hold <= star) {
       goto L11;
	}

   if (abs(*gold) < 0.01) { /* printf(" =0.0 ");*/ *gold =0.01;}; // to protect
   if (abs(hold) < 0.001) {    hold =0.001;}; // 04/20/2001


	flu = (basicdata1.theta * flu + basicdata1.theta1 * *cu * (*tt - basicdata1.tf) / *gold) * .5f;
	fld = (basicdata1.theta * fld + basicdata1.theta1 * *cd * (*td - basicdata1.tf) / hold) * .5f;



L11:
	f = ft - flu - fld;
	if (dabs(f) <= err) {
       goto L50;
	}

	df = f - fwas;
	if (df == 0.f) {
	    goto L50;
	}
  if(abs(df) < 0.01)  {df = 0.01;}  // to protect, zbl edit
	gnew = *g - (*g - gwas) * f / df;
	gwas = *g;
	fwas = f;
	*g = gnew;
   }

return 0;
/* --- NO LATENT HEAT */
L30:
    err = 1e-5f;
    for (n = 1; n <= 50; ++n) {
	h__ = *dxx - *g;
	flu = 0.f;
	if (*cu > 0.f) {
     // added 11/Dec/2000

    if ( (*cu + (c1 + d1 * *g) * *g) ==0.0) {   *cu = -(c1 + d1 * *g) * *g + 0.01;}; // to protect
	  flu = *cu * (a1p + b1p * *g) / (*cu + (c1 + d1 * *g) * *g);
	}
	fld = 0.f;
	if (*cd > 0.f) {
    if ((*cd + (c2 + d2 * h__) * h__)==0.0) {  *cd = -(c2 + d2 * h__) * h__ + 0.001;}; // to protect 04/20/2001

	    fld = *cd * (a2p + b2p * h__) / (*cd + (c2 + d2 * h__) * h__);
	}
	f = -flu - fld;
	if (dabs(f) <= err) {
	    goto L50;
	}
	df = f - fwas;

	if (df == 0.f) {
	    goto L50;
	}

	gnew = *g - (*g - gwas) * f / df;
	gwas = *g;
	fwas = f;
	*g = gnew;
    }

  return 0;

L50:


    if (abs(*g * (c1 + d1 * *g) + *cu) < 0.01) { *cu = -*g * (c1 + d1 * *g) + 0.01;}; // to proctect
    if (abs(h__ * (c2 + d2 * h__) + *cd) <0.01) {  *cd = -h__ * (c2 + d2 * h__) + 0.001;}; // to protect 04/20/2001
    *t1 = (*g * (a1 + b1 * *g) + cutf) / (*g * (c1 + d1 * *g) + *cu);
    *t2 = (h__ * (a2 + b2 * h__) + cdtf) / (h__ * (c2 + d2 * h__) + *cd);

    if (*g >= 0.f) {
       goto L51;
    }

    if (abs(*gold - *g) <0.01)  {
      basicdata1.dtfaz = basicdata1.dtday * 86400.f; // debugging by Q.Z. 12/Dec/2000
       }
     else  basicdata1.dtfaz = *gold * basicdata1.dtfaz / (*gold - *g);

//     if ((*gold - *g) ==0.0)  { *gold = *g +0.01; printf(" bad");};  // 04/20/2001

//     basicdata1.dtfaz = *gold * basicdata1.dtfaz / (*gold - *g);

    *t1 = *tt;
    return 0;
L51:

    if (*g <= *dxx) {
	     return 0;
    }
//  if ((*g - *gold) > 10.0)   { *g = *gold +1.0;}

    if ((*g - *gold) ==0.0) {
    //     *g = *gold +0.01;
     basicdata1.dtfaz = basicdata1.dtday * 86400.f; // debugging by Q.Z. 12/Dec/2000
     } // to protect
    else
    basicdata1.dtfaz = hold * basicdata1.dtfaz / (*g - *gold);
    *t2 = *td;

    return 0;
 } /* phase_ */


/* **************************************************************
************************************************************** */
//zxd: find two neighbor nodes containing a depth (unit: m)
double Soilthermal::cal_out_soilt(double soilt[211], double profile[211], const double& depth) {
	int idx = 0;
	double Tsoil = 0.0;
	for (int i=9; i<210; i++)
	{
		if (depth>=profile[i] && depth<=profile[i+1]) {idx=i; break;}
	} //get index of upper node
	
	if (abs(profile[idx+1]-profile[idx]) < 0.01) {Tsoil = soilt[idx+1];}
	else {Tsoil=(soilt[idx+1]*(depth-profile[idx])+soilt[idx]*(profile[idx+1]-depth))/(profile[idx+1]-profile[idx]);}
	//printf("idx = %i, tsoil = %.3f, profile[idx] = %.3f, profile[idx + 1] = %.3f, tidx = %.3f\n", idx, Tsoil, profile[idx], profile[idx + 1], soilt[idx]);
	return Tsoil;
};

//zbl : reset values
void Soilthermal::reset() 
{
	//printf("3331: reset\n");
	flag = 0;
	profileflag = 0;
	//printflag = 1;
	ism19 = 0;
	i__ = 0;
	airt19 = 0.0; airt29 = 0.0; airt39 = 0.0; tsoil = 0.0; frontd = 0.0; 
  thawbegin = 0.0; thawend = 0.0; smass9 = 0.0; zz = 0.0; hsnow19 = 0.0; hsnow29 = 0.0; hsnow39 = 0.0; 
	int i;
	int j;
	for (i = 0; i < 211; i++) {
	t9[i] = 0.0; //printf("reset: i = %i, t9 = %.3f\n", i, t9[i]);
	x9[i] = 0.0;
	water9[i] = 0.0;
	dx9[i] = 0.0;
	xfa9[i] = 0.0;
	xfb9[i] = 0.0;
	mytstart[i] = 0.0;
	myxstart[i] = 0.0;
	}
	for (i = 0; i < 10; i++) {
	weight9[i] = 0.0;
	diffsoilt[i] = 0.0;
	mytiso[i] = 0.0;
	}
	for (i = 0; i < 20; i++) {
	mydeptem[i] = 0.0;
	mydepflx[i] = 0.0;
	mythick[i] = 0.0;
	mydxa[i] = 0.0;
	mydxb[i] = 0.0;
	mydense[i] = 0.0;
	mywat[i] = 0.0;
	myvcond[i] = 0.0;
	myvsph[i] = 0.0;
	mymat[i] = 0.0;
	}
	for (i = 0; i < 30; i++) {
		for (j = 0; i < 20; i++) {
		mycondt[i][j] = 0.0;
		myspht[i][j] = 0.0;
		mycondf[i][j] = 0.0;
		mysphf[i][j] = 0.0;
		}
	} 
	is9 = 0;
	calcu_cdsnow = 0.0; // to calculate snow thermal conductivity
	water1 = 0.0;	// moss layer water
	water2 = 0.0;	// orgainc layer water
	water3 = 0.0;	// mineral layer water
	tmax = 0.0;
	tmin = 0.0;
	tave = 0.0;
	my_first = 0.0; my_final = 0.0; my_per = 0.0; my_dtday = 0.0; my_theta = 0.0; my_hlat = 0.0; my_tf = 0.0; my_gflux = 0.0;
	my_ig = 0; my_nst = 0;
	my1_cdsnow = 0.0;
	mykiso = 0; mykallyr = 0; myknodes = 0; myktemp = 0; myksnow = 0; mykenvl = 0; mykflux = 0;
	myliso = 0; mylmax = 0;mynanmax = 0;myldepf = 0;
	myvdep = 0.0;myvdep1 = 0.0;
	myvdens = 0.0;
	myvspace = 0.0;mytop = 0.0;
	myepsmin = 0.0;
	mykinti = 0;
	myvdep2 = 0.0;mydephet = 0.0;
	myindex = 0; mynp = 0;
	myvdepth = 0.0;
	myindex2 = 0;
	myepssno = 0.0;myconvrt = 0.0;myeta0 = 0.0;mydenfac = 0.0;myfcmelt = 0.0;mydenmax = 0.0;
	//printf("3435: airt19 = %.3f, airt29 = %.3f, airt39 = %.3f\n", airt19, airt29, airt39);

};

void Soilthermal::getsnowecd(ofstream& rflog1) {

  const int NUMVAR = 52;
  char ecd[80], dummy[NUMVAR][10];
  ifstream infile;
  int i;

  int vegid[NUMVEG], update[NUMVEG];
  char vegname[NUMVEG][31];

  cout << "Enter name of the snow (.ECD) data file with the parameter values:" << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the snow (.ECD) data file with the parameter values:" << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < NUMVEG+1; i++) {
//    infile >> vegid[i] >> vegname[i] >> kc[i] >> ki[i] >> gva[i] >> tmin[i] >> toptmin[i] >> toptmax[i] >> tmax[i] >> raq10a0[i] >> raq10a1[i] >> raq10a2[i] >> raq10a3[i] >> kn1[i] >> labncon[i]  >> leafmxc[i] >> kleafc[i] >> sla[i] >> cov[i] >> fpcmax[i] >> update[i];
    infile >>vegid[i] >> vegname[i]>> basicf.MAX[i] >> basicf.NST[i]>> basicf.KALLYR[i]>>
     basicf.KNODES[i] >> basicf.KISO[i]>> basicf.KTEMP[i] >> basicf.KSNOW[i] >> basicf.KENVL[i]>> basicf.KFLUX[i]
     >> basicf.LISO[i]>> basicf.TISO[i]>> basicf.LMAX[i]>> basicf.VDEPTH[i]>> basicf.DEPTEM1[i]>>
     basicf.DEPTEM2[i]>> basicf.DEPTEM3[i]>> basicf.DEPTEM4[i]>> basicf.DEPTEM5[i]>> basicf.NDEPF[i]>>
     basicf.VDEP[i]>> basicf.DEPFLX1[i]>> basicf.DEPFLX2[i]>> basicf.DEPFLX3[i]>> basicf.DEPFLX4[i]>>
     basicf.DEPFLX5[i]>> basicf.HLAT[i]>>basicf.TF[i]>> basicf.gflux[i]>> basicf.cdsnow[i]>> basicf.FIRST[i]>>
     basicf.FINAL[i]>> basicf.PER[i]>> basicf.DTDAY[i]>> basicf.THETA[i]>> basicf.TOP[i]>> basicf.IG[i]>>basicf.EPSMIN[i]
     >> basicf.VSPACE[i]>> basicf.VDEN[i]>>basicf.kint[i]>>basicf.VDEP1[i]>>basicf.DEPHET[i]>> basicf.SNOFAL[i]>>
     basicf.EPSSNO[i]>> basicf.CONVRT[i]>> basicf.ETAO[i]>>basicf.DENFAC[i]>> basicf.FCMELT[i]>>
     basicf.DENMAX[i]>>update[i];
    }


  infile.close();
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Soilthermal::getsnowecd(char ecd[80]) {

  const int NUMVAR = 52;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i;

  int vegid[NUMVEG], update[NUMVEG];
  char vegname[NUMVEG][31];

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << endl << "Cannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < NUMVEG+1; i++) {
//    infile >> vegid[i] >> vegname[i] >> kc[i] >> ki[i] >> gva[i] >> tmin[i] >> toptmin[i] >> toptmax[i] >> tmax[i] >> raq10a0[i] >> raq10a1[i] >> raq10a2[i] >> raq10a3[i] >> kn1[i] >> labncon[i] >> leafmxc[i] >> kleafc[i] >> sla[i] >> cov[i] >> fpcmax[i] >> update[i];
    infile >>vegid[i] >> vegname[i]>> basicf.MAX[i] >> basicf.NST[i]>> basicf.KALLYR[i]>>
     basicf.KNODES[i] >> basicf.KISO[i]>> basicf.KTEMP[i] >> basicf.KSNOW[i] >> basicf.KENVL[i]>> basicf.KFLUX[i]
     >> basicf.LISO[i]>> basicf.TISO[i]>> basicf.LMAX[i]>> basicf.VDEPTH[i]>> basicf.DEPTEM1[i]>>
     basicf.DEPTEM2[i]>> basicf.DEPTEM3[i]>> basicf.DEPTEM4[i]>> basicf.DEPTEM5[i]>> basicf.NDEPF[i]>>
     basicf.VDEP[i]>> basicf.DEPFLX1[i]>> basicf.DEPFLX2[i]>> basicf.DEPFLX3[i]>> basicf.DEPFLX4[i]>>
     basicf.DEPFLX5[i]>> basicf.HLAT[i]>>basicf.TF[i]>> basicf.gflux[i]>> basicf.cdsnow[i]>> basicf.FIRST[i]>>
     basicf.FINAL[i]>> basicf.PER[i]>> basicf.DTDAY[i]>> basicf.THETA[i]>> basicf.TOP[i]>> basicf.IG[i]>>basicf.EPSMIN[i]
     >> basicf.VSPACE[i]>> basicf.VDEN[i]>>basicf.kint[i]>>basicf.VDEP1[i]>>basicf.DEPHET[i]>> basicf.SNOFAL[i]>>
     basicf.EPSSNO[i]>> basicf.CONVRT[i]>> basicf.ETAO[i]>>basicf.DENFAC[i]>> basicf.FCMELT[i]>>
     basicf.DENMAX[i]>>update[i];
  }
  infile.close();

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Soilthermal::getsoillecd(ofstream& rflog1) {

  const int NUMVAR = 69;
  char ecd[80], dummy[NUMVAR][8];
  ifstream infile;
  int i;
	string comments;
 
  int lfvegid[NUMVEG], update[NUMVEG];
  char lfvegname[NUMVEG][31];

  cout << "Enter name of soil layer (.ECD) data file with parameter values:" << endl;
  fpara >> ecd;getline(fpara, comments);

  rflog1 << "Enter name of soil layer (.ECD) data file with parameter values:" << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < NUMVEG+1; i++) {
//    infile >> lfvegid[i] >> lfvegname[i] >> minleaf[i] >> aleaf[i] >> bleaf[i] >> cleaf[i] >> update[i];
    infile >> lfvegid[i] >> lfvegname[i]>> basicf.THICK1[i]>> basicf.DXA1[i]>>
     basicf.DXB1[i]>> basicf.MAT1[i]>> basicf.DENSE1[i]>> basicf.WATER1[i]>>
     basicf.vcond1[i]>> basicf.vsph1[i]>> basicf.cond1[i]>> basicf.spht1[i]>>
     basicf.condf1[i]>> basicf.sphf1[i]>> basicf.THICK2[i]>>basicf.DXA2[i]>>
     basicf.DXB2[i]>> basicf.MAT2[i]>>basicf.DENSE2[i]>>basicf.WATER2[i]>>
     basicf.vcond2[i]>>basicf.vsph2[i]>>basicf.cond2[i]>>basicf.spht2[i]>>
     basicf.condf2[i]>> basicf.sphf2[i]>> basicf.THICK3[i]>> basicf.DXA3[i]>>
     basicf.DXB3[i]>> basicf.MAT3[i]>> basicf.DENSE3[i]>>basicf.WATER3[i]>>
     basicf.vcond3[i]>> basicf.vsph3[i]>> basicf.cond3[i]>> basicf.spht3[i]>>
     basicf.condf3[i]>> basicf.sphf3[i]>>basicf.THICK4[i]>>basicf.DXA4[i]>>
     basicf.DXB4[i]>> basicf.MAT4[i]>> basicf.DENSE4[i]>>basicf.WATER4[i]>>
     basicf.vcond4[i]>> basicf.vsph4[i]>>basicf.cond4[i]>> basicf.spht4[i]>>
     basicf.condf4[i]>> basicf.sphf4[i]>> basicf.THICK5[i]>>basicf.DXA5[i]>>
     basicf.DXB5[i]>>basicf.MAT5[i]>>basicf.DENSE5[i]>>basicf.WATER5[i]>>
     basicf.vcond5[i]>>basicf.vsph5[i]>> basicf.cond5[i]>>basicf.spht5[i]>>
     basicf.condf5[i]>>basicf.sphf5[i]>>basicf.THICK6[i]>>basicf.DXA6[i]>>
     basicf.DXB6[i]>>basicf.MAT6[i]>>basicf.DENSE6[i]>>basicf.WATER6[i]>>update[i];
   }
  infile.close();

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Soilthermal::getsoillecd(char ecd[80]) {

  const int NUMVAR = 69;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;

  int lfvegid[NUMVEG], update[NUMVEG];
  char lfvegname[NUMVEG][31];

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < NUMVEG+1; i++) {
//    infile >> lfvegid[i] >> lfvegname[i] >> minleaf[i] >> aleaf[i] >> bleaf[i] >> cleaf[i] >> update[i];
    infile >> lfvegid[i] >> lfvegname[i]>> basicf.THICK1[i]>> basicf.DXA1[i]>>
     basicf.DXB1[i]>> basicf.MAT1[i]>> basicf.DENSE1[i]>> basicf.WATER1[i]>>
     basicf.vcond1[i]>> basicf.vsph1[i]>> basicf.cond1[i]>> basicf.spht1[i]>>
     basicf.condf1[i]>> basicf.sphf1[i]>> basicf.THICK2[i]>>basicf.DXA2[i]>>
     basicf.DXB2[i]>> basicf.MAT2[i]>>basicf.DENSE2[i]>>basicf.WATER2[i]>>
     basicf.vcond2[i]>>basicf.vsph2[i]>>basicf.cond2[i]>>basicf.spht2[i]>>
     basicf.condf2[i]>> basicf.sphf2[i]>> basicf.THICK3[i]>> basicf.DXA3[i]>>
     basicf.DXB3[i]>> basicf.MAT3[i]>> basicf.DENSE3[i]>>basicf.WATER3[i]>>
     basicf.vcond3[i]>> basicf.vsph3[i]>> basicf.cond3[i]>> basicf.spht3[i]>>
     basicf.condf3[i]>> basicf.sphf3[i]>>basicf.THICK4[i]>>basicf.DXA4[i]>>
     basicf.DXB4[i]>> basicf.MAT4[i]>> basicf.DENSE4[i]>>basicf.WATER4[i]>>
     basicf.vcond4[i]>> basicf.vsph4[i]>>basicf.cond4[i]>> basicf.spht4[i]>>
     basicf.condf4[i]>> basicf.sphf4[i]>> basicf.THICK5[i]>>basicf.DXA5[i]>>
     basicf.DXB5[i]>>basicf.MAT5[i]>>basicf.DENSE5[i]>>basicf.WATER5[i]>>
     basicf.vcond5[i]>>basicf.vsph5[i]>> basicf.cond5[i]>>basicf.spht5[i]>>
     basicf.condf5[i]>>basicf.sphf5[i]>>basicf.THICK6[i]>>basicf.DXA6[i]>>
     basicf.DXB6[i]>>basicf.MAT6[i]>>basicf.DENSE6[i]>>basicf.WATER6[i]>>update[i];
  }
  infile.close();
};

/* **************************************************************
************************************************************** */

void Soilthermal::showsnowecd(int& cmnt) {

 // clrscr();
  cout << endl << "             SNOW PARAMETERS FOR VEGETATION TYPE STANDS " << endl << endl;
  cout <<endl;
  cout <<endl;

  printf("MAX     = %d     NST    = %d     KALLYR = %d     KNODES = %d     KISO   = %d\n", basicf.MAX[cmnt], basicf.NST[cmnt], basicf.KALLYR[cmnt], basicf.KNODES[cmnt],basicf.KISO[cmnt]);
  printf("KTEMP   = %5.1lf KSNOW  = %5.1lf KENVL  = %5.1lf KFLUX  = %5.1lf\n", basicf.KTEMP[cmnt], basicf.KSNOW[cmnt], basicf.KENVL[cmnt], basicf.KFLUX[cmnt]);
  printf("LISO    = %5.1lf TISO   = %5.1lf LMAX   = %d     VDEPTH = %5.1lf\n", basicf.LISO[cmnt], basicf.TISO[cmnt], basicf.LMAX[cmnt], basicf.VDEPTH[cmnt]);
  printf("DEPTEM1 = %5.1lf DEPTEM2= %5.1lf DEPTEM3= %5.1lf DEPTEM4= %5.1lf DEPTEM5= %5.1lf\n", basicf.DEPTEM1[cmnt], basicf.DEPTEM2[cmnt], basicf.DEPTEM3[cmnt],basicf.DEPTEM4[cmnt],basicf.DEPTEM5[cmnt]);
  printf("NDEPF   = %5.1lf VDEP   = %5.1lf DEPFLX1= %5.1lf DEPFLX2= %5.1lf DEPFLX3= %5.1lf\n", basicf.NDEPF[cmnt], basicf.VDEP[cmnt], basicf.DEPFLX1[cmnt],basicf.DEPFLX2[cmnt],basicf.DEPFLX3[cmnt]);
  printf("DEPFLX4 = %5.1lf DEPFLX5= %5.1lf HLAT   = %5.1lf TF     = %5.1lf GFLUX  = %5.1lf\n",basicf.DEPFLX4[cmnt],basicf.DEPFLX5[cmnt],basicf.HLAT[cmnt], basicf.TF[cmnt], basicf.gflux[cmnt]);
  printf("CDSNOW  = %5.1lf FIRST  = %5.1lf FINAL  = %5.1lf PER    = %5.1lf DTDAY  = %5.1lf\n", basicf.cdsnow[cmnt],basicf.FIRST[cmnt],basicf.FINAL[cmnt],basicf.PER[cmnt], basicf.DTDAY[cmnt]);
  printf("THETA   = %5.1lf TOP    = %5.1lf IG     = %5.1lf EPSMIN = %5.1lf VSPACE = %5.1lf\n", basicf.THETA[cmnt],basicf.TOP[cmnt],basicf.IG[cmnt],basicf.EPSMIN[cmnt],basicf.VSPACE[cmnt]);
  printf("VDEN    = %5.1lf KINT   = %5.1lf VDEP1  = %5.1lf DEPHET = %5.1lf SNOFAL = %5.1lf\n", basicf.VDEN[cmnt], basicf.kint[cmnt],basicf.VDEP1[cmnt],basicf.DEPHET[cmnt],basicf.SNOFAL[cmnt]);
  printf("EPSSNO  = %5.1lf CONVRT = %5.1lf ETAO   = %5.1lf DENFAC = %5.1lf FCMELT = %5.1lf\n", basicf.EPSSNO[cmnt], basicf.CONVRT[cmnt], basicf.ETAO[cmnt],basicf.DENFAC[cmnt],basicf.FCMELT[cmnt]);
  printf("DENMAX  = %5.1lf \n",  basicf.DENMAX[cmnt]);
  cout << endl << "Press any key to continue . . ." << endl;
 // while (getch() == '\0');

};

/* **************************************************************
************************************************************** */

void Soilthermal::showsoillecd(int& cmnt) {
  // clrscr();
  cout << endl << "         PARAMETERS FOR SOIL LAYERS " << endl << endl;
  printf("THICK1 =%5.1f  DXA1 =%5.1f DXB1 =%5.1f MAT1 =%d    DENSE1=%5.1f WATER1= %5.3f\n",basicf.THICK1[cmnt],basicf.DXA1[cmnt],basicf.DXB1[cmnt],basicf.MAT1[cmnt],basicf.DENSE1[cmnt],basicf.WATER1[cmnt]);
  printf("VCOND1 =%5.1f  VSPH1=%5.1f COND1=%5.1f SPHT1=%5.1f CONDF1=%5.1f SPHF1 = %5.1f\n",basicf.vcond1[cmnt],basicf.vsph1[cmnt],basicf.cond1[cmnt],basicf.spht1[cmnt],basicf.condf1[cmnt],basicf.sphf1[cmnt]);
  printf("THICK2 =%5.1f  DXA2 =%5.1f DXB2 =%5.1f MAT2 =%d    DENSE2=%5.1f WATER2= %5.3f\n",basicf.THICK2[cmnt],basicf.DXA2[cmnt],basicf.DXB2[cmnt],basicf.MAT2[cmnt],basicf.DENSE2[cmnt],basicf.WATER2[cmnt]);
  printf("VCOND2 =%5.1f  VSPH2=%5.1f COND2=%5.1f SPHT2=%5.1f CONDF2=%5.1f SPHF2 = %5.1f\n",basicf.vcond2[cmnt],basicf.vsph2[cmnt],basicf.cond2[cmnt],basicf.spht2[cmnt],basicf.condf2[cmnt],basicf.sphf2[cmnt]);
  printf("THICK3 =%5.1f  DXA3 =%5.1f DXB3 =%5.1f MAT3 =%d    DENSE3=%5.1f WATER3= %5.3f\n",basicf.THICK3[cmnt],basicf.DXA3[cmnt],basicf.DXB3[cmnt],basicf.MAT3[cmnt],basicf.DENSE3[cmnt],basicf.WATER3[cmnt]);
  printf("VCOND3 =%5.1f  VSPH3=%5.1f COND3=%5.1f SPHT3=%5.1f CONDF3=%5.1f SPHF3 = %5.1f\n",basicf.vcond3[cmnt],basicf.vsph3[cmnt],basicf.cond3[cmnt],basicf.spht3[cmnt],basicf.condf3[cmnt],basicf.sphf3[cmnt]);
  printf("THICK4 =%5.1f  DXA4 =%5.1f DXB4 =%5.1f MAT4 =%d    DENSE4=%5.1f WATER4= %5.3f\n",basicf.THICK4[cmnt],basicf.DXA4[cmnt],basicf.DXB4[cmnt],basicf.MAT4[cmnt],basicf.DENSE4[cmnt],basicf.WATER4[cmnt]);
  printf("VCOND4 =%5.1f  VSPH4=%5.1f COND4=%5.1f SPHT4=%5.1f CONDF4=%5.1f SPHF4 = %5.1f\n",basicf.vcond4[cmnt],basicf.vsph4[cmnt],basicf.cond4[cmnt],basicf.spht4[cmnt],basicf.condf4[cmnt],basicf.sphf4[cmnt]);
  printf("THICK5 =%5.1f  DXA5 =%5.1f DXB5 =%5.1f MAT5 =%d    DENSE5=%5.1f WATER5= %5.3f\n",basicf.THICK2[cmnt],basicf.DXA2[cmnt],basicf.DXB2[cmnt],basicf.MAT2[cmnt],basicf.DENSE2[cmnt],basicf.WATER2[cmnt]);
  printf("VCOND5 =%5.1f  VSPH5=%5.1f COND5=%5.1f SPHT5=%5.1f CONDF5=%5.1f SPHF5 = %5.1f\n",basicf.vcond2[cmnt],basicf.vsph2[cmnt],basicf.cond2[cmnt],basicf.spht2[cmnt],basicf.condf2[cmnt],basicf.sphf2[cmnt]);
  printf("THICK6 =%5.1f  DXA6 =%5.1f DXB6 =%5.1f MAT6 =%d    DENSE6=%5.1f WATER6= %5.3f\n",basicf.THICK6[cmnt],basicf.DXA6[cmnt],basicf.DXB6[cmnt],basicf.MAT6[cmnt],basicf.DENSE6[cmnt],basicf.WATER6[cmnt]);
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');

};


void Soilthermal::getsoiltecd(ofstream& rflog1) {

  const int NUMVAR = 56;
  char ecd[80], dummy[NUMVAR][10];
  ifstream infile;
  int i,j;
	string comments;
  int vegid[NUMVEG], update[NUMVEG];
  char vegname[NUMVEG][31];

  cout << "Enter name of the soil temperature initial (.ECD) data file with the parameter values:" << endl;
  fpara >> ecd;

  rflog1 << "Enter name of the soil temperature initial (.ECD) data file with the parameter values:" << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < MAXCMNT; i++) {
    infile >>vegid[i] >> vegname[i]>> basicf.INDEX[i]>> basicf.VDEPP[i]>> basicf.NP[i];
      //infile >>dummy[0] >> dummy[0]>> dummy[0]>> dummy[0]>> dummy[0];
      //printf("tsoil 3947: i = %i, idx = %i, vdepp = %.3f, np = %.3f\n", i, basicf.INDEX[i], basicf.VDEPP[i], basicf.NP[i]);
    for (j=0;j<25; j++)
    infile >> basicf.DEPTH[i][j] >> basicf.TEMP[i][j];
    infile >>update[i];
    //printf("tsoil 3950: i = %i, j = %i, depth = %.3f, temp = %.3f\n", i, j, basicf.DEPTH[i][j], basicf.TEMP[i][j]);
    }
  infile.close();
};
 
void Soilthermal::getsoiltecd(char ecd[80]) {

  const int NUMVAR = 56;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i,j;

  int vegid[NUMVEG], update[NUMVEG];
  char vegname[NUMVEG][31];

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << endl << "Cannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 1; i < NUMVEG+1; i++) {
    infile >>vegid[i] >> vegname[i]>> basicf.INDEX[i]>> basicf.VDEPP[i]>> basicf.NP[i];
    for (j=0;j<25; j++)
    infile >> basicf.DEPTH[i][j] >> basicf.TEMP[i][j];
    infile >>update[i];
   }
  infile.close();
 };

void Soilthermal::showsoiltecd(int& cmnt) {
  int j;
//  clrscr();
  cout << "INITIALIZATION FOR SOIL TEMPERATURES:" << endl;
  printf(" INDEX = %5.1f  VDEP = %5.1f  NP = %d \n", basicf.INDEX[cmnt], basicf.VDEPP[cmnt],basicf.NP[cmnt]);
  printf("DEPTH     TEMP \n");
    for (j=0;j<15; j++)
     { cout << basicf.DEPTH[cmnt][j] << "    " << basicf.TEMP[cmnt][j]<< endl;
     }
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');

  printf("DEPTH     TEMP \n");
    for (j=15;j<25; j++)
     { cout << basicf.DEPTH[cmnt][j] << "    " << basicf.TEMP[cmnt][j]<< endl;
    }
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');


};


