/* **************************************************************
datPERCN.CPP - object to read and write the structure of permafrost cn
               data from/to files
************************************************************** */
using namespace std;
#if !defined(DATPERCN_H)
  #include "datpercn.hpp"
#endif

/* *********************************************************** */

Percndata::Percndata(void)
{

  elvend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Percndata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> permacd1;
  infile >> permacd2;
  infile >> permacd3;
  infile >> permacd4;
  infile >> permacd5; 
  infile >> permand1;
  infile >> permand2;
  infile >> permand3;
  infile >> permand4;
  infile >> permand5; 
  infile  >> source;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { elvend = -1;}

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Percndata::getdel(FILE* infile)
{

  elvend = fscanf(infile,"%f,%f, %s ,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf, %s , %s ",
                  &col,&row,varname,&carea,&permacd1,&permacd2,&permacd3,&permacd4,&permacd5,&permand1,&permand2,&permand3,&permand4,&permand5,source,contnent);
                 
  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Percndata::out(ofstream& ofile, float col, float row, char varname[9], int carea, 
                                     double permacd1, double permacd2, double permacd3, double permacd4, double permacd5, 
                                     double permand1, double permand2, double permand3, double permand4, double permand5,
                                     char source[9],char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << permacd1 << ' ';
  ofile << permacd2 << ' ';
  ofile << permacd3 << ' ';
  ofile << permacd4 << ' '; 
  ofile << permacd5 << ' ';
  ofile << permand1 << ' ';
  ofile << permand2 << ' ';
  ofile << permand3 << ' ';
  ofile << permand4 << ' ';
  ofile << permand5 << ' ';
  ofile << source << ' ';  
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Percndata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea,double permacd1, double permacd2, double permacd3, double permacd4, double permacd5, 
                      double permand1, double permand2, double permand3, double permand4, double permand5, char source[9], char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << permacd1 << ",";
  ofile << permacd2 << ",";
  ofile << permacd3 << ",";
  ofile << permacd4 << ",";
  ofile << permacd5 << ","; 
  ofile << permand1 << ",";
  ofile << permand2 << ",";
  ofile << permand3 << ",";
  ofile << permand4 << ",";
  ofile << permand5 << ", ";
  ofile << source << " , ";
  ofile << contnent;
  ofile << endl;

};

