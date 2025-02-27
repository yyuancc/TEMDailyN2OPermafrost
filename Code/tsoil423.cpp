/* **************************************************************
*****************************************************************
TSOIL423.CPP - object describing general characteristics of soil
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */

#if !defined(TSOIL423_H)
  #include "tsoil423.hpp"
#endif

/* **************************************************************
************************************************************** */

Tsoil4::Tsoil4(void)
{

  text  = -99;
  wsoil = -99;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the soil (.ECD) data file with parameter values: ";
  std::cout << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;
  
  rflog1 << "Enter name of the soil (.ECD) data file with parameter values: ";
  rflog1 << ecd << std::endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd (char ecd[80])
{

  char dummy[12];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    std::cerr << "\nCannot open " << ecd << " for data input" << std::endl;
    exit(-1);
  }

  infile >> dummy >> dummy >> dummy;
  infile >> dummy >> pctpora >> update;
  infile >> dummy >> pctporb >> update;
  infile >> dummy >> fldcapa >> update;
  infile >> dummy >> fldcapb >> update;
  infile >> dummy >> wiltpta >> update;
  infile >> dummy >> wiltptb >> update;

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the data file containing the rooting depths:";
  std::cout << std::endl;
  std::cout << "               (e.g., ROOTZVEG.ECD)" << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the data file containing the rooting depths:";
  rflog1 << std::endl;
  rflog1 << "               (e.g., ROOTZVEG.ECD)" << std::endl;
  rflog1 << ecd << std::endl;

  getrootz(ecd);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(char ecd[80])
{

  const int NUMVAR = 7;
  char dummy[NUMVAR][10];
  ifstream infile;

  int i;
  int dcmnt;
  int  rootveg[MAXCMNT];
  long update[MAXCMNT];
  char vegname[MAXCMNT][31];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    std::cerr << "\nCannot open " << ecd << " for data input" << std::endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> rootveg[dcmnt] >> vegname[dcmnt];
    infile >> rootza[dcmnt] >> rootzb[dcmnt] >> rootzc[dcmnt];
    infile >> minrootz[dcmnt] >> update[dcmnt];
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::lake(double& tair,double& prec,double& rain,double& snowfall,
       		  double& pet, double& eet, int& dm)
{

  rgrndh2o[dm] = 0.0;
  sperc[dm] = 0.0;
  snowpack[dm] = 0.0;
  sgrndh2o[dm] = 0.0;
  moist[dm] = 0.0;

  if (tair >= -1.0)
  {
   rain = prec;
    snowfall = 0.0;
  }
  else
  {
    rain = 0.0;
    snowfall = prec;
  }

  eet = pet;
  h2oyld[dm] = prec - pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::percol(double& rain, double& snowinf, double& eet,
                    double& avlh2o, const int& dm)
{

  double extra;
  double recharge;
  sperc[dm] = 0.0;
  rperc[dm] = 0.0;

  recharge = rain + snowinf;
  if (recharge <= 0.0) { recharge = 0.001; }
  if ((avlh2o + rain + snowinf - eet) > awcapmm)
  {
    extra = rain + snowinf + avlh2o - awcapmm - eet;
    sperc[dm] = snowinf * extra / recharge;
    rperc [dm] = rain * extra / recharge;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Tsoil4::rrunoff(const double& rgrndh2o, const double& rperc)
{

  double rrunof;

  rrunof = 0.5 * (rgrndh2o + rperc);

  return rrunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tsoil4::showecd(void)
{

  std::cout << std::endl << "                   SOIL CHARACTERISTICS OF SITE";
  std::cout << std::endl << std::endl;
  printf("PSAND    = %5.2lf      PSILT = %5.2lf      PCLAY = %5.2lf\n",
         pctsand, pctsilt, pctclay);

  printf("POROSITY = %5.2lf   PCFLDCAP = %5.2lf   PCWILTPT = %5.2lf\n",
         pctpor, pcfldcap, pcwiltpt);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::snowmelt(const double& elev, const double& tair,
                        const double& prevtair, const double& snowpack)
{

  double snowflux = 0.0;

  if (tair >= -1.0)
  {
    if (elev <= 500.0) { snowflux = snowpack;}
    else
    {
      if (prevtair < -1) { snowflux = 0.5 * snowpack; }
      else { snowflux = snowpack; }
    }
  }

  return snowflux;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::srunoff(const double& elev, const double& tair,
                       const double& prevtair, const double& prev2tair,
                       const double& sgrndh2o, const double& sperc)
{

  double srunof = 0.0;

  if (tair >= -1.0)
  {
    if (prevtair < -1.0) { srunof = 0.1 * (sgrndh2o + sperc); }
    else
    {
      if (prev2tair < -1)
      {
	if (elev <= 500.0) { srunof = 0.5 * (sgrndh2o + sperc); }
	else { srunof = 0.25 * (sgrndh2o + sperc); }
      }
      else { srunof = 0.5 * (sgrndh2o + sperc); }
    }
  }

  return srunof;

};

/* *************************************************************
************************************************************* */
// To determine the parameters in DOC modelling, refer to (Neff et al. 2001) in Ecosystems 2001, 29-48.
void Tsoil4::docmb(const double& soilorgc, const double& density)
{
	double soilc; //soil c content in terms of %
	soilc = soilorgc / (density * 10000)*10; //transfer soilorgc(g/m2) to soil c content(%), with the assumption that soil c pools are based on 1m depth  added *10 by YYuan 011122 (100 in paper, *10 to match the value in paper)
	//printf("tsoil 311: soilc=%.3f soilorgc=%.3f, density = %.3f\n", soilc, soilorgc, density);
  docm = 0.15 * log(soilc) + 0.51;
  if (docm < 0.0) {docm = 0.0;}
	if (docm > 1.0) docm = 1.0;
	docb = 0.05 * soilc + 0.09;
	if (docb > 1.0) docb = 1.0;	
	docb *= soilc / (soilc + 1.0);
};

/* *************************************************************
************************************************************* */
// To determine the parfrost thawing extra C 
double Tsoil4::permac(const double& frontddiff,const double& frontdmx,const double&frontdymx)
{
/*   double permacd1 = 5;
   double permacd2 = 4;
   double permacd3 = 3;
   double permacd4 = 2;   
   double permacd5 = 1; */
  
   //printf("tsoil333: frontdiff=%.3f,frontdmx=%.3f,frontdymx=%.3f\n", frontddiff,frontdmx, frontdymx);
   //printf("tsoil334: permacd1=%.3f,permacd2=%.3f,permacd3=%.3f,permacd4=%.3f,permacd5=%.3f\n", permacd1,permacd2,permacd3,permacd4,permacd5);
   double permadc;
    if (frontddiff > 0){
       if (frontdmx < 0 &&frontdmx > -0.3 ) {permadc = frontddiff * permacd1 * 1000;}
       else if (frontdmx <= -0.3 && frontdmx > - 0.5 ) {
         if (frontdymx > -0.3) {permadc = ((-0.3-frontdmx) * permacd2  + (frontdymx + 0.3) * permacd1) * 1000;}
         else {permadc = frontddiff * permacd2 * 1000;}
       }      
       else if (frontdmx <= -0.5 && frontdmx > - 1.0 ) { 
         if (frontdymx > -0.3) {permadc = ((-0.5 - frontdmx) * permacd3 + 0.2 * permacd2 + (frontdymx + 0.3) * permacd1) * 1000;}
         else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadc = ((-0.5 - frontdmx) * permacd3 + (frontdymx + 0.5) * permacd2) * 1000;}
         else {permadc =frontddiff * permacd3 * 1000;}
       }
       else if (frontdmx <= -1.0 && frontdmx > - 2.0 ) {  
         if (frontdymx > -0.3) {permadc = ((-1.0 - frontdmx) * permacd4 + 0.5 * permacd3 + 0.2 * permacd2 + (frontdymx + 0.3) * permacd1) * 1000;}
         else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadc = ((-1.0 - frontdmx) * permacd4 + 0.5 * permacd3 + (frontdymx + 0.5) * permacd2) * 1000;}
         else if (frontdymx > -1.0 && frontdymx <= -0.5 ) {permadc = ((-1.0 - frontdmx) * permacd4 + (frontdymx + 1.0) * permacd3) * 1000;}                                              
         else {permadc = frontddiff * permacd4 * 1000;}
       }
       else if (frontdmx <= -2.0) { 
          if (frontdymx > -0.3) {permadc = ((-2.0 - frontdmx) * permacd5 + 1.0 * permacd4 + 0.5 * permacd3 + 0.2 * permacd2 + (frontdymx + 0.3) * permacd1) * 1000;}
          else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadc = ((-2.0 - frontdmx) * permacd5 + 1.0 * permacd4 + 0.5 * permacd3 + (frontdymx + 0.5) * permacd2) * 1000;}
          else if (frontdymx > -1.0 && frontdymx <= -0.5 ) {permadc = ((-2.0 - frontdmx) * permacd5 + 1.0 * permacd4 + (frontdymx + 1.0) * permacd3) * 1000;} 
          else if (frontdymx > -2.0 && frontdymx <= -1.0 ) {permadc = ((-2.0 - frontdmx) * permacd5 + (frontdymx + 2.0) * permacd4) * 1000;}                                                       else {permadc = frontddiff * permacd5 * 1000;}   
        }               
    } 
    return permadc;


};
/* *************************************************************
************************************************************* */
// To determine the parfrost thawing extra N 
double Tsoil4::perman(const double& frontddiff,const double& frontdmx,const double&frontdymx)
{
   /*double permand1 = 5;
   double permand2 = 4;
   double permand3 = 3;
   double permand4 = 2;   
   double permand5 = 1; */
  double permadn;
   //printf("tsoil373: frontdiff=%.3f,frontdmx=%.3f,frontdymx=%.3f\n", frontddiff,frontdmx, frontdymx);
   //printf("tsoil375: permand1=%.3f,permand2=%.3f,permand3=%.3f,permand4=%.3f,permand5=%.3f\n", permand1,permand2,permand3,permand4,permand5);
  if (frontddiff > 0){
       if (frontdmx < 0 &&frontdmx > -0.3 ) {permadn = frontddiff * permand1 * 1000;}
       else if (frontdmx <= -0.3 && frontdmx > - 0.5 ) {
         if (frontdymx > -0.3) {permadn = ((-0.3-frontdmx) * permand2  + (frontdymx + 0.3) * permand1) * 1000;}
         else {permadn = frontddiff * permand2 * 1000;}
       }      
       else if (frontdmx <= -0.5 && frontdmx > - 1.0 ) { 
         if (frontdymx > -0.3) {permadn = ((-0.5 - frontdmx) * permand3 + 0.2 * permand2 + (frontdymx + 0.3) * permand1) * 1000;}
         else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadn = ((-0.5 - frontdmx) * permand3 + (frontdymx + 0.5) * permand2) * 1000;}
         else {permadn =frontddiff * permand3 * 1000;}
       }
       else if (frontdmx <= -1.0 && frontdmx > - 2.0 ) {  
         if (frontdymx > -0.3) {permadn = ((-1.0 - frontdmx) * permand4 + 0.5 * permand3 + 0.2 * permand2 + (frontdymx + 0.3) * permand1) * 1000;}
         else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadn = ((-1.0 - frontdmx) * permand4 + 0.5 * permand3 + (frontdymx + 0.5) * permand2) * 1000;}
         else if (frontdymx > -1.0 && frontdymx <= -0.5 ) {permadn = ((-1.0 - frontdmx) * permand4 + (frontdymx + 1.0) * permand3) * 1000;}                                       
         else {permadn = frontddiff * permand4 * 1000;}
       }
       else if(frontdmx <= -2.0) { 
          if (frontdymx > -0.3) {permadn = ((-2.0 - frontdmx) * permand5 + 1.0 * permand4 + 0.5 * permand3 + 0.2 * permand2 + (frontdymx + 0.3) * permand1) * 1000;}
          else if (frontdymx > -0.5 && frontdymx <= -0.3 ) {permadn = ((-2.0 - frontdmx) * permand5 + 1.0 * permand4 + 0.5 * permand3 + (frontdymx + 0.5) * permand2) * 1000;}
          else if (frontdymx > -1.0 && frontdymx <= -0.5 ) {permadn = ((-2.0 - frontdmx) * permand5 + 1.0 * permand4 + (frontdymx + 1.0) * permand3) * 1000;} 
          else if (frontdymx > -2.0 && frontdymx <= -1.0 ) {permadn = ((-2.0 - frontdmx) * permand5 + (frontdymx + 2.0) * permand4) * 1000;}                                       
          else {permadn = frontddiff * permand5 * 1000;}   
        }                
    } 
    return permadn;

};

/* *************************************************************
************************************************************** */

void Tsoil4::xtext(int& cmnt, double& pctsilt, double& pctclay)
{

  totpor = fldcap = wiltpt = MISSING;
  awcapmm =  MISSING;

  psiplusc = (pctsilt + pctclay) * 0.01;
  if (psiplusc < 0.01) { psiplusc = 0.01; }

  rootz = (rootza[cmnt] * pow(psiplusc, 2.0)) + (rootzb[cmnt] * psiplusc)
          + rootzc[cmnt];
  if (rootz < minrootz[cmnt]) { rootz = minrootz[cmnt]; }

  pctpor = (pctpora * psiplusc) + pctporb;
  pcfldcap = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor  = rootz * pctpor * 10.0;
  fldcap  = rootz * pcfldcap * 10.0;
  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm = fldcap - wiltpt;

};

/* *************************************************************
************************************************************* */
// added for 3-layer hydrological model, Q. Zhuang, 05/Jan/2003

void Tsoil4::hydm_xtext(int& ez, double& pctsilt, double& pctclay)
 {
 double rootmin;
 double rootmx;

  totpor = fldcap = wiltpt = MISSING;
  awcapmm =  MISSING;

  rootmin = 99.9;
  rootmx = -999.9;

  totpor1 = fldcap1 = wiltpt1 = -999.9;
  totpor2 = fldcap2 = wiltpt2 = -999.9;
  totpor3 = fldcap3 = wiltpt3 = -999.9;

  awcapmm1 =  -999.9;
  awcapmm2 =  -999.9;
  awcapmm3 =  -999.9;

  psiplusc = (pctsilt + pctclay) * 0.01;

  if (psiplusc < 0.01) { psiplusc = 0.01; }

// unit is meter

  rootz = (rootza[ez] * pow(psiplusc, 2.0)) + (rootzb[ez] * psiplusc) + rootzc[ez];

  if (rootz < minrootz[ez]) { rootz = minrootz[ez]; }

  if (rootz< rootmin) {rootmin= rootz;}
  if (rootz> rootmx) { rootmx=rootz;}

  dpwbox1 = rootz * 0.1;  // a thin layer (surface laler 0.1m for the first root zone, maybe the moss, Desborough, 1997
  dpwbox2 = rootz * 0.7; // needed to reconsider it
  dpwbox3 = rootz * 0.2; // mineral account for 20%

  pctpor = (pctpora * psiplusc) + pctporb;
  pcfldcap = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor1  = rootz * 0.1 * pctpor * 10.0 ;
  fldcap1  = rootz * 0.1 * pcfldcap * 10.0;
  wiltpt1  = rootz * 0.1 * pcwiltpt * 10.0;

  totpor2  = rootz * 0.8 * pctpor * 10.0 ;
  fldcap2  = rootz * 0.8 * pcfldcap * 10.0;
  wiltpt2  = rootz * 0.8 * pcwiltpt * 10.0;

  totpor3  = rootz * 0.2 * pctpor * 10.0 * 5;
  fldcap3  = rootz * 0.2 * pcfldcap * 10.0;
  wiltpt3  = rootz * 0.2 * pcwiltpt * 10.0;

  totpor  = rootz * pctpor * 10.0;
  fldcap  = rootz * pcfldcap * 10.0;
  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm1 = fldcap1 - wiltpt1;
  awcapmm2 = fldcap2 - wiltpt2;
  awcapmm3 = fldcap3 - wiltpt3;

  awcapmm = fldcap - wiltpt;  // estimate capacity for daily hydrological model

};

/* *************************************************************
************************************************************** */

