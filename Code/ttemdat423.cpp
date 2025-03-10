/* *************************************************************
TTEMDAT423.CPP - object to read and write the structure of the
               output data from files generated by the Terrestrial 
               Ecosystem Model (TEM)
************************************************************* */
using namespace std;
#if !defined(TTEMDAT423_H)
  #include "ttemdat423.hpp"
#endif

/* ********************************************************** */

Temdata::Temdata(void)
{

  temend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Temdata::get(ifstream& ifile)
{

  int i;

  lagpos = ifile.tellg();

  ifile >> col >> row;
  ifile >> varname;
  ifile >> tveg;
  ifile >> subtveg;
  ifile >> psiplusc;
  ifile >> qlcon;
  ifile >> carea;
  ifile >> year;

  ifile >> total;
  ifile >> max;
  ifile >> ave;
  ifile >> min;

  for (i = 0; i < CYCLE; i++)  { ifile >> mon[i]; }
  ifile >> contnent;

  ifile.seekg(0, ios::cur);
  curpos = ifile.tellg();

  if (curpos < (lagpos + 10)) { temend = -1; }

  return temend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Temdata::getdel(FILE* ifile)
{

  temend = fscanf(ifile, "%f,%f, %s ,%d,%d,%lf,%d,%d,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf, %s",
                  &col, &row, varname, &tveg, &subtveg, &psiplusc, &qlcon,
                  &carea, &year, &total, &max, &ave, &min,
                  mon+0, mon+1, mon+2, mon+3, mon+4, mon+5,
                  mon+6, mon+7, mon+8, mon+9, mon+10, mon+11, contnent);
  return temend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Temdata::getyrsum(double mon[CYCLE])
{

  int j;

  total = 0.0;
  max   = -900000.0;
  min   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > max) { max = mon[j]; }
    if (mon[j] < min) { min = mon[j]; }
    total += mon[j];
  }

  ave = total / (double) CYCLE;
  if (total <= ((MISSING*CYCLE)+CYCLE)) { total = MISSING; }


};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Temdata::out(ofstream& ofile, float col, float row, char varname[9],
                  int tveg, int subtveg, double psiplusc, int qlcon,
                  int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= ((MISSING*CYCLE)+CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << tveg << ' ';
  ofile << subtveg << ' ';
  ofile << setprecision(2) << psiplusc << ' ';
  ofile << setprecision(0) << qlcon << ' ';
  ofile << carea << ' ';
  ofile << year << ' ';

  ofile << setprecision(1) << predtotl << ' ';
  ofile << predmax << ' ';
  ofile << setprecision(2) << predave << ' ';
  ofile << setprecision(1) << predmin << ' ';

  for (int i = 0; i < CYCLE; i++) { ofile << mon[i] << ' '; }

  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Temdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                     int tveg, int subtveg, double psiplusc, int qlcon,
                     int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= ((MISSING*CYCLE)+CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);
  //ofile.precision(2);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << tveg << ",";
  ofile << subtveg << ",";
  ofile << setprecision(2) << psiplusc << ",";
  ofile << setprecision(0) << qlcon << ",";
  ofile << carea << ",";
  ofile << year << ",";

  ofile << setprecision(1) << predtotl << ",";
  ofile << predmax << ",";
  ofile << setprecision(2) << predave << ",";
  ofile << setprecision(1) << predmin << ",";

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ","; }

  ofile << " " << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Temdata::pctout(ofstream& ofile, float col, float row, char varname[9],
                     int tveg, int subtveg, double psiplusc, int qlcon,
                     int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= ((MISSING*CYCLE)+CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << tveg << ' ';
  ofile << subtveg << ' ';
  ofile << setprecision(2) << ' ';
  ofile << psiplusc << ' ';
  ofile << setprecision(0) << qlcon << ' ';
  ofile << carea << ' ';
  ofile << year << ' ';

  ofile << setprecision(2) << predtotl << ' ';
  ofile << predmax << ' ';
  ofile << predave << ' ';
  ofile << predmin << ' ';

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ' '; }

  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Temdata::poutdel(ofstream& ofile, float col, float row, char varname[9],
                      int tveg, int subtveg, double psiplusc, int qlcon,
                      int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= ((MISSING*CYCLE)+CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << tveg << ",";
  ofile << subtveg << ",";
  ofile << setprecision(2) << psiplusc << ",";
  ofile << setprecision(0) << qlcon << ",";
  ofile << carea << ",";
  ofile << year << ",";

  ofile << setprecision(2) << predtotl << ",";
  ofile << predmax << ",";
  ofile << predave << ",";
  ofile << predmin << ",";

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ","; }

  ofile << " " << contnent;
  ofile << endl;

};

void Temdata::dayoutdel(ofstream& ofile, float col, float row, char varname[9],
                     int tveg, int subtveg, double psiplusc, int qlcon,
                     int carea, long year, double dat[CYCLE][31], char contnent[9], int days[CYCLE])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  int m;

  for (j = 0; j < CYCLE; j++) {

    predtotl = 0.0;
//    predmax   = -900000.0;
    predmax   = -999999.9;
    predmin   =  900000.0;

//     days = daysnumber(j, year);
     for (m=0; m < days[j]; m++)
      {
      if (dat[j][m] > predmax) { predmax = dat[j][m]; }
      if (dat[j][m] < predmin) { predmin = dat[j][m]; }
      predtotl += dat[j][m];
     }
     if (predtotl <= ((MISSING*days[j])+days[j])) { predave = predtotl = MISSING; }
     else predave = predtotl / (double) days[j];

   ofile.setf(ios::fixed,ios::floatfield);
   ofile.setf(ios::showpoint);
   //ofile.precision(1);
   ofile.precision(2);

   ofile << col << "," << row << ", ";
   ofile << varname << " ,";
   ofile << setprecision(0) << tveg << ",";
   ofile << subtveg << ",";
   ofile << setprecision(2) << psiplusc << ",";
   ofile << setprecision(0) << qlcon << ",";
   ofile << carea << ",";
   ofile << year << ",";
   ofile << (j+1) << ",";

   ofile << setprecision(1) << predtotl << ",";
   ofile << predmax << ",";
   ofile << setprecision(2) << predave << ",";
   ofile << setprecision(1) << predmin << ",";
    for (m=0; m < days[j]; m++)
      {
      ofile << setprecision(2) << dat[j][m] << ",";
      }  // end of days
   ofile << contnent;
   ofile << endl;

  }
};

