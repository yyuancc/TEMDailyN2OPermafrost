/* **************************************************************
DATPERMA.HPP - object to read and write the structure of permafrost C,N
               data from/to files
************************************************************** */
using namespace std;
#if !defined(DATPERCN_H)
#define DATPERCN_H

class Percndata {

  public:

     Percndata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[9], int carea, double permacd1, double permacd2, double permacd3, double permacd4, double permacd5, 
                                                                                 double permand1, double permand2, double permand3, double permand4, double permnad5, char source[9], char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], int carea, double permacd1, double permacd2, double permacd3, double permacd4, double permacd5, 
                                                                                    double permand1, double permand2, double permand3, double permand4, double permand5, char source[9], char contnent[9]);


/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // "perma"
     int carea;           // area covered by grid cell (sq. km)
     double permacd1;      // C content in layer 1 0-30
     double permacd2;      // C content in layer 2 30-50
     double permacd3;      // C content in layer 3 50-100
     double permacd4;      // C content in layer 4 100-200
     double permacd5;      // C content in layer 5 200+     
     double permand1;      // N content in layer 1 0-30
     double permand2;      // N content in layer 2 30-50
     double permand3;      // N content in layer 3 50-100
     double permand4;      // N content in layer 4 100-200
     double permand5;      // N content in layer 5 200+     
     char source[9];      // reference to data source
     char contnent[9];    // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int elvend;
     long curpos;
     long lagpos;

};

#endif

