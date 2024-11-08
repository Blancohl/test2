#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
  BOR_CORD VERSION 2.2
  
  This program (developed for BOREAS) performs the geographic coordinate 
  conversions over the BOREAS study area of 51 N --> 60 N latitude and
  111 W --> 93 W longitude.
      
  The BOREAS grid is based on the ellipsoidal version of the Albers Equal-Area
  Conic (AEAC) projection as defined within the North American Datum 1983
  (NAD83).  The origin of the grid is at 111 W, 51 N and the standard
  parallels are set to 52.5 N and 58.5 N as prescribed in 'Map Projections -
  A Working Manual', USGS Professional Paper 1395, John P. Snyder, 1987.
  The projection equations used in this program were all taken from
  this manual.
  
  The software will take input coordinates as Universal Transverse Mercator
  (UTM) northing and easting; decimal degrees of latitude and longitude;
  or BOREAS (x,y) grid.
  
  Output coordinates for the BOREAS grid are always given under the NAD83
  datum.  Output coordinates for UTM and lat, lon coordinates are given under
  both NAD27 and NAD83 datums.

  The software handles conversion between NAD27 and NAD83 by using the
  function F to convert lat, lon coordinates between the two datums.
  The c_gridint function accesses datum shift information contained in the
  file datmshft.dat.  This file must be locatable for the software to work as
  designed.
  
  The BOREAS project thanks the Geodetic Survey of Canada for the FORTRAN-77
  source code and datum shift information files that made this possible.
===============================================================================
  The original development of this program was the result of collaboration
  between
              Jeffrey Newcomer, Scott Goetz, and Joe Lewthwaite
              Hughes STX Corporation    September 1992
-------------------------------------------------------------------------------
   The current version 2.0 is the result of efforts by:
        Jeffrey Newcomer retrofitted the version 1.3 code to procedurally
        handle the change over between datums.
        
        Fred Irani adapted the software in the function c_gridint that handles
        conversion between NAD27 and NAD83 datums from the National
        Transformation Program.
        
 -------------------------------------------------------------------------------
 Update: 10/94 F. Irani, HSTX: Added precision to longitude and latitude output.
 Update: 12/94 F. Irani, HSTX: X,Y Output now fits input format.
===============================================================================
*/

/*
  Make some definitions for machine portability.
*/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
  If you are running on a PC or VAX system, you will want to have
  the following definition activated so that reading of the file
  datmshft.dat by the c_gridint function works properly:

          #define  NRECBYTES 22


  If you are running on a Silicon Graphics or SUN work-station under 
  UNIX or on a Macintosh system, you will want to have the following 
  definition activated so that reading of the file datmshft.dat by 
  the c_gridint function works properly:

          #define  NRECBYTES 21
*/

#define  NRECBYTES 22 /* VAX, PC, etc */


/*
    North American Datum 27 (NAD27) (also Clarke 1866) parameters
*/
#define   A_27      6378206.4      /* NAD27 equatorial radius */
#define   E_27      0.0822719      /* NAD27 eccentricity */
#define   E2_27     (E_27*E_27)    /* NAD27 eccentricity squared */

/*
    North American Datum 83 (NAD83) parameters
*/
#define   A_83      6378137.0      /* NAD83 equatorial radius */
#define   E_83      0.0818192      /* NAD83 eccentricity */
#define   E2_83     (E_83*E_83)    /* NAD83 eccentricity squared */

/*
    Other needed map parameters
*/
#define  RAD2DEG  57.29577951308  /* Radian to degree conversion */
#define  DEG2RAD   0.01745329252  /* Degree to radian conversion */
#define  K0        0.9996         /* Scale factor on UTM Central Meridian */
#define  M0        0.0000         /* Distance along meridian */
                                  /* from equator to point */
/*
  Function to output software version banner
*/
void version_num (void);

/*
  Functions that control coordinate transformations
*/
void process_coords (FILE *, FILE *, FILE *, long int, char *, char *, char *, char *);
long int process_boreas (double, double, double, double,
                         FILE *, FILE *, FILE *);
long int process_utm27 (double, double, double, double,
                        FILE *, FILE *, FILE *);
long int process_utm83 (double, double, double, double,
                        FILE *, FILE *, FILE *);
long int process_latlon27 (double, double, double, double,
                           FILE *, FILE *, FILE *);
long int process_latlon83 (double, double, double, double,
                           FILE *, FILE *, FILE *);

/*
  Functions to read input coordinate information
*/
long int read_boreas (FILE *, FILE *, double *, double *);
long int read_ll (FILE *, FILE *, double *, double *);
long int read_utm (FILE *, FILE *, double *, double *, long int *);

/*
  Functions to perform basic coordinate transformations
*/
long int albers_2_boreas (double, double, double *, double *);
long int albers_2_ll (double, double, double, double, long int, double *,
                      double *, double *, double *);
long int boreas_2_albers (double, double, double *, double *);
void ll_2_albers (double, double, double, double, long int, double *,
                      double *, double *, double *);
long int ll_2_utm (double *, double *, double, double, double, long int,
                   double *, double *, long int *);
long int utm_2_ll (double *, double *, long int *, double, double, double,
                   long int, double *, double *);
              
/*
  Function to perform datum coordinate conversions
*/
long int c_gridint (FILE *, char *, long int, double *, double *,
                    double *, double *);

/*
  Function to output coordinate information
*/
long int write_coords (FILE *, double, double, double, double, double, double,
                       double, double, long int, double, double, long int);

/*
 Function to check user response for quit 
*/ 
long int check_response (char *);

/*
===============================================================================
   S T A R T   O F   M A I N   F U N C T I O N
===============================================================================
*/ 
main ()
{  

  char checkfile[100], filename[100], input_str[20];
  char title0[200], title1[200], title2[200], title3[200];
 
  FILE *datum_fp, *in_fp = NULL, *out_fp = NULL;

  long int choice, cnt, len;

  int  stat;
  
  /*
    Output banner and determine the coordinate conversions to be performed
  */
  version_num ();

  /*
    As long as the user wants to perform coordinate transformations
    keep prompting them for what transformations to perform.
  */
  for ( ; ; )
  {
    /*
      Keep presenting the menu and trying to get a good response from the user.
      Note: All files should be closed at this point.
    */
    do
    {
      puts ("\n\n  Input Coordinate Selection Menu.");
      puts ("      1) BOREAS X, Y.");
      puts ("      2) UTM  Easting, Northing");
      puts ("      3) Longitude, Latitude");
      puts ("      4) Exit.");
      printf (" Enter choice < 1 - 4 >: ");
      gets (input_str);
      len = strlen (input_str);
      if (len > 0)
      {
        cnt = sscanf (input_str, "%ld", &choice);
        if (cnt != 1)
        {
          puts ("*** Problem interpreting input ***");
          puts ("*** Please re-enter response ***");
          len = 0;
        }
      }
    } while (len == 0);

    if (choice == 4)
    {
      puts ("\n\n< Boreas_Coords Completed >");
      exit (0);
    }
	
    /*
      Do until a useable filename has been entered, <cr> or q(uit).
    */
    stat = TRUE;
    for (; ;)
    {
     /*
        Determine how the coordinate values are to be entered
        and if needed attempt to open the input file
      */
      puts ("\n\nEnter Input Filename or RETURN for keyboard input (q = quit)");
      printf (" --> ");
      gets (filename);
      if (!check_response(filename))
      {
	stat = FALSE;
	break;
     }
      else if (strlen (filename) == 0)
      {
	in_fp = (FILE *) NULL;
        break;
      }
      else if ((in_fp = fopen (filename, "r")) == NULL)
      {
	printf ("\n*** Unable to open %s for input ***\n", filename);
	continue;
      }
      break;
    }
    if (!stat) continue;

    /*
      Do forever until able to 
      determine where to send the converted values and
      attempt to open the output file if needed
    */
    stat = TRUE;
    for (; ;)
    {
      puts ("\n\nEnter Output filename or RETURN for screen output (q = quit)");
      printf (" --> ");
      gets (filename);
      if (!check_response(filename))
      {
	stat = FALSE;
	break;
      }
      else if (strlen (filename) == 0)
      {
       out_fp = (FILE *) NULL;
       break;
      }
      else if ((out_fp = fopen (filename, "w")) == NULL)
      {
	    printf ("\n*** Unable to  open %s for output ***\n", filename);
        continue;
      }
      /*
      Save valid user response for output file title.
      ----------------------------------------------- */
      sprintf(title0, "\nFILE: %s\n\n", filename);
      break;
    }
    if (!stat) continue;
    
    /*
      Open the file that contains the datum shift information
      needed by the c_gridint function for conversions between
      NAD27 and NAD83.
    */
    datum_fp = (FILE *) NULL;
    stat = TRUE;
    while (TRUE)
    {
      strcpy(filename, "datmshft.dat");
      printf ("\n\nEnter datum shift file name (Default = %s; q = quit)\n", filename);
      printf (" --> ");
      gets (filename);
      cnt = sscanf (filename, "%s", checkfile);
      if (cnt < 1) strcpy(filename, "datmshft.dat");
      else if ( !check_response(filename) )
      {
	stat = FALSE;
	break;
      }

      if ((datum_fp = fopen (filename, "r")) == NULL)
      {
	printf ("\n*** Unable to open %s for input ***\n", filename);
        continue;
      }
      else break;
    }
    if (!stat) continue;

    /*
      If an output file was specified, then set up the column titles to
      the top three records of the output file for later use 
   */
    if (out_fp != (FILE *) NULL)
    {
      title1[0] = title2[0] = title3[0] = '\0';
      sprintf(title1, 
      			"\n      * = INPUT\n\n   BOREAS GRID       NAD27 COORDINATES       ");
      strcat(title1,  "NAD83 COORDINATES      NAD27 COORDINATES");
      strcat(title1,  "           NAD83 COORDINATES     \n");

      sprintf(title2, " X (KM)    Y (KM)    LONGITUDE  LATITUDE    ");
      sprintf(title3, " --------  -------   ---------- --------    ");

      strcat(title2,  " LONGITUDE   LATITUDE   EASTING   NORTHING   ZONE");
      strcat(title3,  " ----------  --------   --------  ---------  ----");

      strcat(title2,  "   EASTING   NORTHING   ZONE\n");
      strcat(title3,  "   --------  ---------  ----\n");
    }
    
    /*
      Perform the specified coordinate conversion operations
    */
    process_coords (in_fp, out_fp, datum_fp, choice, title0, title1, title2, title3);

    /*
      Before continuing loop, close up the files
    */
    if (in_fp    != (FILE *) NULL) fclose (in_fp);
    if (out_fp   != (FILE *) NULL) fclose (out_fp);
    if (datum_fp != (FILE *) NULL) fclose (datum_fp);
  	 
  }  /* end of forever for loop */

} /* END OF MAIN */

/*=============================================================================
  Function to control processing of coordinates

  Input parameters:
      in_fp == FILE pointer to file of input coordinate values
      out_fp == FILE pointer to file for output coordinate information
      datum_fp == FILE pointer to file containing datum shift information
      choice == The number that designates the user's
                selected input coordinate system

  Output parameters:
      Returned status value.
=============================================================================*/
void process_coords (in_fp, out_fp, datum_fp, choice, title0, title1, title2, title3)
FILE *in_fp, *out_fp, *datum_fp;
long int choice;
char *title0, *title1, *title2, *title3;
{
  double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
  
  long int cnt, in_datum, stat;
  
  char input_str[80];
  
/*
    Set values for the BOREAS Grid and Albers Equal Area Conic
    (AEAC) projection parameters
*/
  aeac_phi0 = 51.00;        /* Origin latitude for AEAC projection */
                            /*  (southern most) */
  aeac_lambda0 = -111.00;   /* Origin longitude for AEAC projection */
                            /*  (western most) */
  aeac_phi1 = 52.50;        /* Latitude of southern standard parallel */
  aeac_phi2 = 58.50;        /* Latitude of northern standard parallel */

  switch (choice)
  {
    case 1 :  /* BOREAS (x,y) */
    {
      puts ("\n\n\nNote that the BOREAS X and Y grid values are assumed");
      puts ("to have been calculated under the NAD83 datum");
      if (out_fp != (FILE *) NULL)
      {
        title1[20] = '*';
        stat = fprintf(out_fp, title0);
        stat = fprintf(out_fp, title1);
        stat = fprintf(out_fp, title2);
        stat = fprintf(out_fp, title3);
        if (stat < 0)
        {
          puts ("\n*** Error writing header records to output file ***");
          fclose (out_fp);
          if (in_fp != (FILE *) NULL) fclose (in_fp);
          fclose (datum_fp);
          datum_fp = in_fp = out_fp = (FILE *) NULL;
       	  return;
        }
      }
      process_boreas (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                      in_fp, datum_fp, out_fp);
      break;
    }
    case 2 :  /* UTM */
    {
      do
      {
        puts ("\n\nEnter the datum for the input UTM coordinates");
        puts ("as 27 for NAD27 or 83 for NAD83 (q = quit)");
        puts (" --> ");
        if (gets (input_str) == NULL)
        {
          puts ("*** Error reading from terminal ***");
          continue;
        }
	if ( !check_response (input_str) ) return;
    	cnt = sscanf (input_str, "%ld", &in_datum);
        if (cnt != 1)
        {
          puts ("*** Error extracting input datum from input record ***");
          continue;
        }
        if (in_datum != 27 && in_datum != 83)
            puts ("*** Datum values need to be given as 27 or 83 ***");
      } while (in_datum != 27 && in_datum !=83);
      if (in_datum == 27) 
      {
        if (out_fp != (FILE *) NULL)
        {
          title1[89] = '*';
          stat = fprintf(out_fp, title0);
          stat = fprintf(out_fp, title1);
          stat = fprintf(out_fp, title2);
          stat = fprintf(out_fp, title3);
          if (stat < 0)
          {
            puts ("\n*** Error writing header records to output file ***");
            fclose (out_fp);
            if (in_fp != (FILE *) NULL) fclose (in_fp);
            fclose (datum_fp);
            datum_fp = in_fp = out_fp = (FILE *) NULL;
       	    return;
          }
        }
      	process_utm27 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                       in_fp, datum_fp, out_fp);
      }
      else if (in_datum == 83)
      {
        if (out_fp != (FILE *) NULL)
        {
          title1[117] = '*';
          stat = fprintf(out_fp, title0);
          stat = fprintf(out_fp, title1);
          stat = fprintf(out_fp, title2);
          stat = fprintf(out_fp, title3);
          if (stat < 0)
          {
            puts ("\n*** Error writing header records to output file ***");
            fclose (out_fp);
            if (in_fp != (FILE *) NULL) fclose (in_fp);
            fclose (datum_fp);
            datum_fp = in_fp = out_fp = (FILE *) NULL;
       	    return;
          }
        }
      	process_utm83 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                       in_fp, datum_fp, out_fp);
      }
      break;
    }
    case 3 :  /* Latitude, longitude */
    {
      do
      {
        puts ("\n\nEnter the datum for the input latitude and longitude");
        puts ("coordinates as 27 for NAD27 or 83 for NAD83. (q = quit)");
        puts (" --> ");
        if (gets (input_str) == NULL)
        {
          puts ("*** Error reading from terminal ***");
          continue;
        }
        if ( !check_response (input_str) ) return;
        cnt = sscanf (input_str, "%ld", &in_datum);
        if (cnt != 1)
        {
          puts ("*** Error extracting input datum from input record ***");
          continue;
        } 
      } while (in_datum != 27 && in_datum !=83);
      if (in_datum == 27) 
      {
        if (out_fp != (FILE *) NULL)
        {
          title1[39] = '*';
          stat = fprintf(out_fp, title0);
          stat = fprintf(out_fp, title1);
          stat = fprintf(out_fp, title2);
          stat = fprintf(out_fp, title3);
          if (stat < 0)
          {
            puts ("\n*** Error writing header records to output file ***");
            fclose (out_fp);
            if (in_fp != (FILE *) NULL) fclose (in_fp);
            fclose (datum_fp);
            datum_fp = in_fp = out_fp = (FILE *) NULL;
       	    return;
          }
        }
     	process_latlon27 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                          in_fp, datum_fp, out_fp);
      }
      else if (in_datum == 83)
      {
        if (out_fp != (FILE *) NULL)
        {
          title1[63] = '*';
          stat = fprintf(out_fp, title0);
          stat = fprintf(out_fp, title1);
          stat = fprintf(out_fp, title2);
          stat = fprintf(out_fp, title3);
          if (stat < 0)
          {
            puts ("\n*** Error writing header records to output file ***");
            fclose (out_fp);
            if (in_fp != (FILE *) NULL) fclose (in_fp);
            fclose (datum_fp);
            datum_fp = in_fp = out_fp = (FILE *) NULL;
            return;
          }
        }
      	process_latlon83 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                          in_fp, datum_fp, out_fp);
      }        
      break;
    }
  }  /* end of switch */

}  /* END OF PROCESS_COORDS */

/*=============================================================================
   Function to control conversion from BOREAS grid
   locations to latitude and longitude and UTM coordinates.
   
   Note that this function sets the in_datum parameter to a constant value
   of 83 since the BOREAS grid system was based on the NAD83 datum.
   
  Input parameters:
      aeac_phi0 == origin latitude for the projection (southern most)
      aeac_phi1 == latitude of southern standard parallel
      aeac_phi2 == latitude of northern standard parallel
      aeac_lambda0 == origin longitude for the projection (western most)
      in_fp == FILE pointer to file of input coordinate values
      datum_fp == FILE pointer to file containing datum shift information
      out_fp == FILE pointer to file for output coordinate information

  Output parameters:
      Returned status value.
=============================================================================*/
long int process_boreas (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                    in_fp, datum_fp, out_fp)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
FILE *in_fp, *datum_fp, *out_fp;
{

  double lat_27, lat_83, lon_27, lon_83;
  double easting_27, easting_83, northing_27, northing_83;
  double x, xgrid, y, ygrid;

  long int cnt, zone_27, zone_83;

  char datum_str[8], input_str[80];

  /*
    Set the input coordinate datum to NAD83 since the
    BOREAS grid is based on the NAD83 datum.
  */
  strcpy (datum_str, "NAD83");
  
  /*
    While there are coordinates to process, handle them accordingly
  */
  while ( read_boreas (in_fp, out_fp, &xgrid, &ygrid) )
  {
    /*
      Convert BOREAS x and y components into x and y AEAC meter
      values under NAD83
    */
    boreas_2_albers (xgrid, ygrid, &x, &y);

    /*
      Convert AEAC meters to latitude and longitude under NAD83
    */
    albers_2_ll (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, (long int) 1, &x, &y,
                 &lat_83, &lon_83);

    /*
      Convert NAD83 latitude and longitude to NAD83 UTM northing and easting
    */
    ll_2_utm (&lat_83, &lon_83, (double) A_83, (double) E_83, (double) E2_83,
              (long int) 1, &northing_83, &easting_83, &zone_83);
    
    /*
      Convert NAD83 latitude and longitude to NAD27 latitude and longitude    
    */
    if ( !c_gridint (datum_fp, datum_str, (long int) 1, &lat_83, &lon_83,
                     &lat_27, &lon_27) )  return (0);
    
    /*
      Convert NAD27 latitude and longitude to NAD27 UTM northing and easting
    */
    ll_2_utm (&lat_27, &lon_27, (double) A_27, (double) E_27, (double) E2_27,
              (long int) 1, &northing_27, &easting_27, &zone_27);
    
    /*
      Output the information as indicated
    */
    if (!write_coords (out_fp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                       northing_27, easting_27, zone_27, northing_83,
                       easting_83, zone_83) )
    {
      puts ("*** Error writing coordinates to output ***");
      return (0);
    }
  } /* end of while loop */
  
  return (1);

}  /* END OF PROCESS_BOREAS */

/*=============================================================================
   Function to control conversion from Universal Transverse
   Mercator (UTM) under NAD27 to latitude, longitude and BOREAS grid

  Input parameters:
      aeac_phi0 == origin latitude for AEAC projection (southern most)
      aeac_phi1 == latitude of southern standard parallel for AEAC projection
      aeac_phi2 == latitude of northern standard parallel for AEAC projection
      aeac_lambda0 == origin longitude for the AEAC projection (western most)
      in_fp == FILE pointer to file of input coordinate values
      datum_fp == FILE pointer to file containing datum shift information
      out_fp == FILE pointer to file for output coordinate information

  Output parameters:
      Returned status value.
=============================================================================*/
long int process_utm27 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                   in_fp, datum_fp, out_fp)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
FILE *in_fp, *datum_fp, *out_fp;
{
  double lat_27, lat_83, lon_27, lon_83;
  double easting_27, easting_83, northing_27, northing_83;
  double temp_north, temp_east, x, xgrid, y, ygrid;

  long int cnt, zone_27, zone_83;

  char datum_str[8], input_str[80];

  /*
    Set the input coordinate datum to NAD27.
  */
  strcpy (datum_str, "NAD27");
  
  /*
    While there are coordinates to process, handle them accordingly
  */
  while ( read_utm (in_fp, out_fp, &northing_27, &easting_27, &zone_27) )
  {
    /*
      Convert UTM northing and easting to latitude and longitude
      under NAD27 datum
    */
    utm_2_ll ( &northing_27, &easting_27, &zone_27, (double) A_27,
               (double) E_27, (double) E2_27, (long int) 1, &lat_27, &lon_27 );

    /*
      Convert NAD27 latitude and longitude to NAD83 latitude and longitude
    */
    if ( !c_gridint (datum_fp, datum_str, (long int) 1, &lat_27, &lon_27,
                     &lat_83, &lon_83) ) return (0);
        
    /*
      Convert NAD83 latitude and longitude to NAD83 UTM northing and easting
    */
    ll_2_utm (&lat_83, &lon_83, (double) A_83, (double) E_83, (double) E2_83,
              (long int) 1, &northing_83, &easting_83, &zone_83);
    
    /*
      Convert NAD83 latitude and longitude values to AEAC meter values
    */
    ll_2_albers (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, (long int) 1,
                 &lat_83, &lon_83, &x, &y);

    /*
      Convert AEAC meters to BOREAS grid coordinates
    */
    albers_2_boreas (x, y, &xgrid, &ygrid);

    /*
      Output the information as indicated
    */
    if (!write_coords (out_fp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                       northing_27, easting_27, zone_27, northing_83,
                       easting_83, zone_83) )
    {
      puts ("*** Error writing coordinates to output ***");
      return (0);
    }
  } /* end of while loop */
  
  return (1);

}  /* END OF PROCESS_UTM27 */

/*=============================================================================
   Function to control conversion from Universal Transverse
   Mercator (UTM) under NAD83 to latitude, longitude  and BOREAS grid

  Input parameters:
      aeac_phi0 == origin latitude for AEAC projection (southern most)
      aeac_phi1 == latitude of southern standard parallel for AEAC projection
      aeac_phi2 == latitude of northern standard parallel for AEAC projection
      aeac_lambda0 == origin longitude for the AEAC projection (western most)
      in_fp == FILE pointer to file of input coordinate values
      datum_fp == FILE pointer to file containing datum shift information
      out_fp == FILE pointer to file for output coordinate information

  Output parameters:
      Returned status value.
=============================================================================*/
long int process_utm83 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                   in_fp, datum_fp, out_fp)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
FILE *in_fp, *datum_fp, *out_fp;
{
  double lat_27, lat_83, lon_27, lon_83;
  double easting_27, easting_83, northing_27, northing_83;
  double temp_north, temp_east, x, xgrid, y, ygrid;

  long int cnt, zone_27, zone_83;

  char datum_str[8], input_str[80];

  /*
    Set the input coordinate datum to NAD83.
  */
  strcpy (datum_str, "NAD83");
  
  /*
    While there are coordinates to process, handle them accordingly
  */
  while ( read_utm (in_fp, out_fp, &northing_83, &easting_83, &zone_83) )
  {
    /*
      Convert UTM northing and easting to latitude and longitude
      under NAD83 datum
    */
    utm_2_ll ( &northing_83, &easting_83, &zone_83, (double) A_83,
               (double) E_83, (double) E2_83, (long int) 1, &lat_83, &lon_83 );

    /*
      Convert NAD83 latitude and longitude to NAD27 latitude and longitude
    */
    if ( !c_gridint (datum_fp, datum_str, (long int) 1, &lat_83, &lon_83,
                     &lat_27, &lon_27) )  return (0);
        
    /*
      Convert NAD27 latitude and longitude to NAD27 UTM northing and easting
    */
    ll_2_utm (&lat_27, &lon_27, (double) A_27, (double) E_27, (double) E2_27,
              (long int) 1, &northing_27, &easting_27, &zone_27);
    
    /*
      Convert NAD83 latitude and longitude values to AEAC meter values
    */
    ll_2_albers (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, (long int) 1,
                 &lat_83, &lon_83, &x, &y);

    /*
      Convert AEAC meters to BOREAS grid components
    */
    albers_2_boreas (x, y, &xgrid, &ygrid);

    /*
      Output the information as indicated
    */
    if (!write_coords (out_fp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                       northing_27, easting_27, zone_27, northing_83,
                       easting_83, zone_83) )
    {
      puts ("*** Error writing coordinates to output ***");
      return (0);
    }
  } /* end of while loop */
  
  return (1);

}  /* END OF PROCESS_UTM83 */

/*=============================================================================
   Function to control conversion from latitude and
   longitude under NAD27 to BOREAS grid and UTM coordinates.

  Input parameters:
      aeac_phi0 == origin latitude for AEAC projection (southern most)
      aeac_phi1 == latitude of southern standard parallel for AEAC projection
      aeac_phi2 == latitude of northern standard parallel for AEAC projection
      aeac_lambda0 == origin longitude for the AEAC projection (western most)
      in_fp == FILE pointer to file of input coordinate values
      datum_fp == FILE pointer to file containing datum shift information
      out_fp == FILE pointer to file for output coordinate information

  Output parameters:
      Returned status value.
=============================================================================*/
long int process_latlon27 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                      in_fp, datum_fp, out_fp)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
FILE *in_fp, *datum_fp, *out_fp;
{
  double lat_27, lat_83, lon_27, lon_83;
  double easting_27, easting_83, northing_27, northing_83;
  double x, xgrid, y, ygrid;

  long int cnt, zone_27, zone_83;

  char datum_str[8], input_str[80];

  /*
    Set the input coordinate datum to NAD27.
  */
  strcpy (datum_str, "NAD27");
  
  /*
    While there are coordinates to process, handle them accordingly
  */
  while ( read_ll (in_fp, out_fp, &lat_27, &lon_27) )
  {
    /*
      Convert NAD27 latitude and longitude to NAD83 latitude and longitude
    */
    if ( !c_gridint (datum_fp, datum_str, (long int) 1, &lat_27, &lon_27,
                     &lat_83, &lon_83) )  return (0);
        
    /*
      Convert NAD83 latitude and longitude values to AEAC meter values
    */
    ll_2_albers (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, (long int) 1,
                 &lat_83, &lon_83, &x, &y);

    /*
      Convert AEAC meters to BOREAS grid components
    */
    albers_2_boreas (x, y, &xgrid, &ygrid);

    /*
      Convert NAD27 latitude and longitude to NAD27 UTM northing and easting
    */
    ll_2_utm (&lat_27, &lon_27, (double) A_27, (double) E_27, (double) E2_27,
              (long int) 1, &northing_27, &easting_27, &zone_27);
    
    /*
      Convert NAD83 latitude and longitude to NAD83 UTM northing and easting
    */
    ll_2_utm (&lat_83, &lon_83, (double) A_83, (double) E_83, (double) E2_83,
              (long int) 1, &northing_83, &easting_83, &zone_83);
    
    /*
      Output the information as indicated
    */
    if (!write_coords (out_fp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                       northing_27, easting_27, zone_27, northing_83,
                       easting_83, zone_83) )
    {
      puts ("*** Error writing coordinates to output ***");
      return (0);
    }
  } /* end of while loop */
  
  return (1);

}  /* END OF PROCESS_LATLON27 */

/*=============================================================================
   Function to control conversion from and
   longitude under NAD83 to BOREAS grid and UTM coordinates.

  Input parameters:
      aeac_phi0 == origin latitude for AEAC projection (southern most)
      aeac_phi1 == latitude of southern standard parallel for AEAC projection
      aeac_phi2 == latitude of northern standard parallel for AEAC projection
      aeac_lambda0 == origin longitude for the AEAC projection (western most)
      in_fp == FILE pointer to file of input coordinate values
      datum_fp == FILE pointer to file containing datum shift information
      out_fp == FILE pointer to file for output coordinate information

  Output parameters:
      Returned status value.
=============================================================================*/
long int process_latlon83 (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0,
                      in_fp, datum_fp, out_fp)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
FILE *in_fp, *datum_fp, *out_fp;
{
  double lat_27, lat_83, lon_27, lon_83;
  double easting_27, easting_83, northing_27, northing_83;
  double x, xgrid, y, ygrid;

  long int cnt, zone_27, zone_83;

  char datum_str[8], input_str[80];

  /*
    Set the input coordinate datum to NAD83.
  */
  strcpy (datum_str, "NAD83");
  
  /*
    While there are coordinates to process, handle them accordingly
  */
  while ( read_ll (in_fp, out_fp, &lat_83, &lon_83) )
  {
    /*
      Convert NAD83 latitude and longitude to NAD27 latitude and longitude
    */
    if ( !c_gridint (datum_fp, datum_str, (long int) 1, &lat_83, &lon_83,
                     &lat_27, &lon_27) )  return (0);
        
    /*
      Convert NAD83 latitude and longitude values to AEAC meter values
    */
    ll_2_albers (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, (long int) 1,
                 &lat_83, &lon_83, &x, &y);

    /*
      Convert AEAC meters to BOREAS grid components
    */
    albers_2_boreas (x, y, &xgrid, &ygrid);

    /*
      Convert NAD27 latitude and longitude to NAD27 UTM northing and easting
    */
    ll_2_utm (&lat_27, &lon_27, (double) A_27, (double) E_27, (double) E2_27,
              (long int) 1, &northing_27, &easting_27, &zone_27);
    
    /*
      Convert NAD83 latitude and longitude to NAD83 UTM northing and easting
    */
    ll_2_utm (&lat_83, &
lon_83, (double) A_83, (double) E_83, (double) E2_83,
              (long int) 1, &northing_83, &easting_83, &zone_83);
    
    /*
      Output the information as indicated
    */
    if (!write_coords (out_fp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                       northing_27, easting_27, zone_27, northing_83,
                       easting_83, zone_83) )
    {
      puts ("*** Error writing coordinates to output ***");
      return (0);
    }
  } /* end of while loop */
  
  return (1);

}  /* END OF PROCESS_LATLON83 */

/*=============================================================================
  Function to convert Albers Equal-Area Conic (x,y) coordinates
  to latitude and longitude coordinates in NAD83 datum.
  
  All input latitude and longitide coordinates are assumed to be
  in decimal degree units.

  Input parameters:
      aeac_phi0 == origin latitude for the projection (southern most)
      aeac_phi1 == latitude of southern standard parallel
      aeac_phi2 == latitude of northern standard parallel
      aeac_lambda0 == origin longitude for the project (western most)
      ncoords == the number of latitude and longitude coordinate pairs
                 contained in the lat and lon arrays
      x == array of AEAC x values in meters
      y == array of AEAC y values in meters

  Output parameters:
      lat == array of latitude coordinates in degrees
      lon == array of longitude coordinates in degrees
=============================================================================*/
long int albers_2_ll (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, ncoords,
                 x, y, lat, lon)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
long int ncoords;
double *x, *y, *lat, *lon;
{

  double c, lambda, m, m_1, m_2, n, phi, q, q_0, q_1, q_2;
  double rho, rho_0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, theta;

  long int i, j;

  /*
    Calculate the needed constants
  */
  temp1 = sin(aeac_phi0*DEG2RAD);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_0 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));
  aeac_lambda0 = aeac_lambda0 * DEG2RAD;

  temp1 = sin(aeac_phi1*DEG2RAD);
  m_1 = cos(aeac_phi1*DEG2RAD)/sqrt(1.0-E2_83*temp1*temp1);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_1 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));

  temp1 = sin(aeac_phi2*DEG2RAD);
  m_2 = cos(aeac_phi2*DEG2RAD)/sqrt(1.0-E2_83*temp1*temp1);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_2 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));

  n = (m_1*m_1 - m_2*m_2)/(q_2 - q_1);
  c = m_1*m_1 + n*q_1;
  rho_0 = A_83*sqrt(c - n*q_0)/n;

  /*
    Calculate lat and lon coordinates from each Albers x, y pair
  */
  for (i = 0; i < ncoords; i++)
  {
    temp1 = x[i];
    temp2 = rho_0 - y[i];
    rho = sqrt( temp1*temp1 + temp2*temp2 );

    theta = atan(temp1/temp2);
  
    temp1 = rho*n/A_83;
    q = (c - temp1*temp1)/n;

    /*
      Set initial value of phi (i.e., lat) and derive iterative solution
    */
    phi = asin (q/2.0);
    lat[i] = phi * RAD2DEG;
    for (j = 0; j < 6; j++)
    {
      temp1 = sin(phi);
      temp2 = E_83*temp1;
      temp3 = temp2*temp2;
      temp4 = 1.0 - temp3;
      temp4 = (temp4*temp4)/(2.0*cos(phi));
      temp5 = q/(1.0 - E2_83);
      temp6 = temp1/(1.0 - temp3);
      temp7 = log( (1.0-temp2)/(1.0+temp2) );
      temp7 = temp7/(2.0*E_83);
      phi = phi + temp4*(temp5 - temp6 + temp7);
      lat[i] = phi * RAD2DEG;
    }

    lambda = aeac_lambda0 + theta/n;
    lon[i] = lambda * RAD2DEG;
  }
  return (0);

}  /* END OF ALBERS_2_LL */

/*=============================================================================
  Function to convert from latitude and longitude under NAD83
  to Albers Equal-Area Conic (x,y) coordinates.
  
  Note that all input latitude and longitide coordinates are assumed to be
  in decimal degree units.  Also the earth radius and eccentricity 
  parameters (A_83, E_83, and E2_83) are all hardcoded here since for the
  BOREAS project, the datum is assumed to have been NAD83.

  Input parameters:
      aeac_phi0 == origin latitude for the projection (southern most)
      aeac_phi1 == latitude of southern standard parallel
      aeac_phi2 == latitude of northern standard parallel
      aeac_lambda0 == origin longitude for the project (western most)
      ncoords == the number of latitude and longitude coordinate pairs
                 contained in the lat and lon arrays
      lat == array of latitude coordinates
      lon == array of longitude coordinates

  Output parameters:
      x == array of AEAC x values in meters
      y == array of AEAC y values in meters
=============================================================================*/
void ll_2_albers (aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0, ncoords,
                 lat, lon, x, y)
double aeac_phi0, aeac_phi1, aeac_phi2, aeac_lambda0;
long int ncoords;
double *lat, *lon, *x, *y;
{
  double c, lambda, m_1, m_2, n, phi, q, q_0, q_1, q_2;
  double rho, rho_0, temp1, temp2, theta;

  long int i;

  /*
    Calculate the needed constants
  */
  temp1 = sin(aeac_phi0*DEG2RAD);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_0 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));
  aeac_lambda0 = aeac_lambda0 * DEG2RAD;

  temp1 = sin(aeac_phi1*DEG2RAD);
  m_1 = cos(aeac_phi1*DEG2RAD)/sqrt(1.0-E2_83*temp1*temp1);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_1 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));

  temp1 = sin(aeac_phi2*DEG2RAD);
  m_2 = cos(aeac_phi2*DEG2RAD)/sqrt(1.0-E2_83*temp1*temp1);
  temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
  q_2 = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));

  n = (m_1*m_1 - m_2*m_2)/(q_2 - q_1);
  c = m_1*m_1 + n*q_1;

  rho_0 = A_83*sqrt(c - n*q_0)/n;

  /*
    Calculate the Albers x and y coordinates for each lat, lon pair
  */
  for (i = 0; i < ncoords; i++)
  {
    phi = lat[i] * DEG2RAD;
    lambda = lon[i] * DEG2RAD;
    theta = n*(lambda - aeac_lambda0);
    temp1 = sin(phi);
    temp2 = (1.0-E_83*temp1)/(1.0+E_83*temp1);
    q = (1.0-E2_83) * (temp1/(1.0-E2_83*temp1*temp1) - log(temp2)/(2.0*E_83));
    rho = A_83*sqrt(c - n*q)/n;
    x[i] = rho * sin(theta) + 0.5;
    y[i] = rho_0 - rho*cos(theta) + 0.5;
  }
}  /* END OF LL_2_ALBERS */

/*=============================================================================
  Function to convert from latitude and longitude
  to UTM northing and easting and zone
  See pp.48-65;269-270, USGS Map Projections Working Manual, J.P.Snyder
  
  Note that input latitude and longitide coordinates are assumed to be
  in decimal degree units.

  Input parameters:
      lat == array of latitude coordinates to be transformed
      lon == array of longitude coordinates to be transformed
      A == semimajor axis of the earth ellipsoid in meters
      E == eccentricity of the assumed ellipsoid
      E2 == eccentricity squared of the assumed ellipsoid
      ncoords == the number of latitude and longitude coordinate pairs
                 contained in the lat and lon arrays

  Output parameters:
      northing == array of northing coordinates
      easting == array of easting coordinates
      zone == array of UTM zones for the coordinates
=============================================================================*/
long int ll_2_utm (lat, lon, A, E, E2, ncoords, northing, easting, zone)
double *lat, *lon, A, E, E2;
long int ncoords;
double *northing, *easting;
long int *zone;
{
  /*
    The flow is from Snyder, solving in order equations:
    8.12 - 8.15, 3.22, 8.9, 8.10
  */

  double lat_rad;
  double e_prime, lambda0;
  double a, c, m, n, t;
  double coef1, coef2, coef3, coef4;

  long int i;

  /*
    Calculate some needed values outside the loop
  */
  e_prime = E2 / (1.0 - E2);  /* Equation 8.12 */
  
  /* The four following values are from Eq 3.21 */
  coef1 = 1.0 - pow(E,2.0)/4.0 - 3.0*pow(E,4.0)/64.0 - 5.0*pow(E,6.0)/256.0;
  coef1 = A * coef1 * DEG2RAD;
  coef2 = A * (3.0*pow(E,2.0)/8.0 + 3.0*pow(E,4.0)/32.0 + 45.0*pow(E,6.0)/1024.0);
  coef3 = A * (15.0*pow(E,4.0)/256.0 + 45.0*pow(E,6.0)/1024.0);
  coef4 = A * (35.0*pow(E,6.0)/3072.0);

  /*
    Calculate the Northing and Easting for each
    pair of input latitude and longitude values
  */
  for (i = 0; i < ncoords; i++)
  {
    /* 
      Assign appropriate central meridian to UTM zone based on input longitude 
    */
    zone[i] = (long int) ( (lon[i]+177.0)/6.0 + 1.5);
    lambda0 = (double) (zone[i]-1)*6 - 177;

    /*
      Equations 8.12 - 8.15
    */
    lat_rad = lat[i] * DEG2RAD;
    n = A / sqrt(1.0 - E2 * sin(lat_rad) * sin(lat_rad));
    t = tan(lat_rad) * tan(lat_rad);
    c = e_prime * cos(lat_rad) * cos(lat_rad);
    a = cos(lat_rad) * ((lon[i] - lambda0) * DEG2RAD);

    /*
      Equation 3.22, The distance along central meridian to given latitude
    m = (111132.0894 * lat[i]) - (16216.94 * sin(2.0 * lat_rad)) 
         + (17.21 * sin(4.0 * lat_rad)) - (0.02 * sin(6.0 * lat_rad));	
    */
    m = (coef1 * lat[i]) - (coef2 * sin(2.0 * lat_rad)) 
         + (coef3 * sin(4.0 * lat_rad)) - (coef4 * sin(6.0 * lat_rad));	

    /*
      Equations 8.9 and 8.10
    */
    easting[i] = K0 * n * (a + (1.0 - t + c) * pow(a,3.0) / 6.0
                    + (5.0 - 18.0 * t + pow(t,2.0) + 72.0 * c - 58.0 * e_prime)
                    * pow(a,5.0) / 120.0) + 500000;
	           
    northing[i] = K0 * (m + n * tan(lat_rad) 
                     * (pow(a,2.0) / 2.0 + (5.0 - t + 9.0 * c + 4.0 * pow(c,2.0))
                     * pow(a,4.0) / 24.0 + (61.0 - 58.0 * t + pow(t,2.0) 
                     + 600.0 * c - 330.0 * e_prime) * pow(a,6.0) / 720.0));
	            
  }  /* end of for */
  
  return(0);

}  /* END OF LL_2_UTM */

/*=============================================================================
  Function to calculate Latitude and Longitude coordinates from UTM coordinates
  See pp.48-65;270-271, USGS Map Projections Working Manual, J.P.Snyder

  Note that output latitude and longitide coordinates are in decimal
  degree units.

  Input parameters:
      northing == array of northing coordinates in decimal meters
      easting == array of easting coordinates in decimal meters
      zone == array of UTM zones for the coordinates
      A == semimajor axis of the earth ellipsoid in meters
      E == eccentricity of the assumed ellipsoid
      E2 == eccentricity squared of the assumed ellipsoid
      ncoords == the number of latitude and longitude coordinate pairs
                 contained in the lat and lon arrays

  Output parameters:
      lat == array of latitude coordinates to be transformed
      lon == array of longitude coordinates to be transformed
=============================================================================*/
long int utm_2_ll (northing, easting, zone, A, E, E2, ncoords, lat, lon)
double *northing, *easting;
long int *zone;
double A, E, E2;
long int ncoords;
double *lat, *lon;
{
  /*
    The flow is from Snyder, solving in order equations:
    8.12, 8.20, 3.24, 7.19, 3.26, 8.21-8.25, 8.17, 8.18 
  */

  double e_prime, east_temp, lambda0, m, mu;
  double phi1, phi1_deg, temp, temp2;
  double c1, d, e1, n1, t1, r1;
	
  long int i;
	
  /*
    Calculate some needed values outside the loop
  */
  e_prime = E2 / (1.0 - E2);
  e1 = (1.0 - sqrt(1.0 - E2)) / (1.0 + sqrt(1.0 - E2));
  temp = (A * (1.0 - E2/4.0 - 3.0*pow(E,4.0)/64.0 - 5.0*pow(E,6.0)/256.0));
  
  /*
    Calculate the latitude and longitude values for each pair of
    input Northing and Easting values
  */
  for (i = 0; i < ncoords; i++)
  {
    /*
      Calculate longitude at Central Meridian for given UTM zone
    */
    lambda0 = (double) (zone[i]-1)*6 - 177;
  
    /*
      Subtract "false" easting (Central Meridian value) from entered values
    */
    east_temp = easting[i] - 500000.0;

    /*
      Calculate values for Eqns 8.12, 8.20, 3.24, 7.19
    */
    m = M0 + northing[i] / K0;
    mu = m / temp;

    /*
      Calculate value for Eqn 3.26
    */
    phi1 = mu + (3.0*e1/2.0 - 27.0*pow(e1,3.0)/32.0) * sin(2.0*mu)
              + (21.0*e1*e1/16.0 - 55.0*pow(e1,4.0)/32.0) * sin(4.0*mu)
              + (151.0*pow(e1,3.0)/96.0) * sin(6.0*mu);
    phi1_deg = phi1*RAD2DEG;
  
    /*
      Calculate values for Eqns 8.21 - 8.25
    */
    c1 = e_prime * cos(phi1) * cos(phi1);
    t1 = tan(phi1) * tan(phi1);
    n1 = A / sqrt(1.0 - E2 * sin(phi1) * sin(phi1));
    temp2 = 1.0 - E2 * sin(phi1) * sin(phi1);
    r1 = A * (1.0 - E2) / pow(temp2,1.5);
    d = east_temp / (n1 * K0);

    /*
      Calculate values for Eqns 8.17 and 8.18
    */
    lat[i]=phi1_deg - (n1 * tan(phi1) / r1) * (RAD2DEG * (pow(d,2.0)/2.0 
                    - (5.0 + 3.0*t1 + 10.0*c1 - 4.0*pow(c1,2.0) 
                    - 9.0*pow(e_prime,2.0)) * pow(d,4.0)/24.0 
                    + (61.0 + 90.0*t1 + 298.0*c1 + 45.0*pow(t1,2.0) 
                    - 252.0*pow(e_prime,2.0) - 3.0*pow(c1,2.0)) * pow(d,6.0)/720.0));
			   
    lon[i] = lambda0 + (RAD2DEG * ( d - (1.0 + 2.0*t1 + c1) * pow(d,3.0)/6.0
                     + (5.0 - 2.0*c1 + 28.0*t1 - 3.0*pow(c1,2.0) + 8.0*e_prime
                     + 24.0*pow(t1,2.0) ) * pow(d,5.0) / 120.0) / cos(phi1));
  }

 return(0);

}  /* END OF UTM2LL */

/*=============================================================================
  Function to convert Albers Equal-Area Conic
  coordinates into the BOREAS (x,y) grid components
  
  Input parameters:
      x == X coordinate in units of AEAC decimal meters
      y == Y coordinate in units of AEAC decimal meters

  Output parameters:
      xgrid == BOREAS X grid coordinate
      ygrid == BOREAS Y grid coordinate
=============================================================================*/
long int albers_2_boreas (x, y, xgrid, ygrid)
double x, y, *xgrid, *ygrid;
{

  *xgrid = x/1000.0;
  *ygrid = y/1000.0;

  return (1);

}  /* END OF ALBERS_2_BOREAS */

/*=============================================================================
  Function to convert BOREAS (x,y) grid component values
  into Albers Equal-Area Conic meters
  
  Input parameters:
      xgrid == BOREAS X grid coordinate in decimal kilometers
      ygrid == BOREAS Y grid coordinate in decimal kilometers

  Output parameters:
      x == X coordinate in units of AEAC decimal meters
      y == Y coordinate in units of AEAC decimal meters
=============================================================================*/
long int boreas_2_albers (xgrid, ygrid, x, y)
double xgrid, ygrid, *x, *y;
{

  *x = xgrid*1000.0;
  *y = ygrid*1000.0;

  return (1);

}  /* END OF BOREAS_2_ALBERS */

/*=============================================================================
  Function to get BOREAS (x,y) coordinate components

  Input parameters:
      infp == FILE pointer to input file of coordinates
      outfp == FILE pointer to output file of coordinates

  Output parameters:
      xgrid == BOREAS grid x coordinate in kilometers (east from origin)
      ygrid == BOREAS grid y coordinate in kilometers (north from origin)
=============================================================================*/
long int read_boreas (infp, outfp, xgrid, ygrid)
FILE *infp, *outfp;
double *xgrid, *ygrid;
{
  long int cnt;
  
  char input_str[80];
  
  /*
    The following do while loop allows reprompting of the
    user for erroneous values and handles returns to calling
    function if the input is a disk file
  */
  cnt = 0;
  do
  {
    /*
      Get an input record from the specified file or user
    */
    if (infp != (FILE *) NULL)
    {
      if (fgets (input_str, 80, infp) == NULL)
      {
        if ( ferror (infp) ) puts ("*** Error reading from input file ***");
        else if ( feof (infp) ) puts ("\n\nEnd of input file reached");
        return (0);
      }
    }
    else
    {
      puts ("\nEnter the X and Y BOREAS grid values as xxx.xx yyy.yy\n");
      puts ("(0 <= X <= 1000     0 <= Y <= 1000) ; q = quit");
      printf (" --> ");
      if (gets (input_str) == NULL)
      {
        puts ("*** Error reading from terminal ***");
        continue;
      }
      if ( !check_response (input_str) ) return (0);
    }

    /*
      See if anything was given in the input string
    */
    if ( strlen(input_str) == 0)
    {
      puts ("No coordinates found on input record");
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }
    /*
      Since it seems like something is in the input string
      try to extract the expected coordinates from it
    */
    cnt = sscanf (input_str, "%lf %lf", xgrid, ygrid);
    if (cnt != 2)
    {
      puts ("*** Error extracting BOREAS coordinates from input record ***");
      cnt = 0;
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }
  
    /*
      Check the ygrid value for validity
    */
    if (*ygrid > 1000.0 || *ygrid < 0.0)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified BOREAS Y coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only BOREAS Y values ***\n", outfp);
        fputs ("*** between 0.0 and 1000.00 are accepted     ***\n", outfp);
        puts ("\n*** Specified BOREAS Y coordinate is outside ***");
        puts ("*** anticipated range.  Only BOREAS Y values ***");
        puts ("*** between 0.0 and 1000.00 are accepted     ***");
        return (0);
      }
      puts ("\n*** Specified BOREAS Y coordinate is outside ***");
      puts ("*** anticipated range.  Only BOREAS Y values ***");
      puts ("*** between 0.0 and 1000.00 are accepted     ***");
      cnt = 0;

    }
    /*
      Check the xgrid value for validity
    */
    if (*xgrid > 1000.00 || *xgrid < 0.0)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified BOREAS X coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only BOREAS X values ***\n", outfp);
        fputs ("*** between 0.0 and 1000.00 are accepted     ***\n", outfp);
        puts ("\n*** Specified BOREAS X coordinate is outside ***");
        puts ("*** anticipated range.  Only BOREAS X values ***");
        puts ("*** between 0.0 and 1000.00 are accepted     ***");
        return (0);

      }
      puts ("\n*** Specified BOREAS X coordinate is outside ***");
      puts ("*** anticipated range.  Only BOREAS X values ***");
      puts ("*** between 0.0 and 1000.00 are accepted     ***");
      cnt = 0;
    }
    
  } while (cnt == 0);  /* end do while loop */
  
  return (1);  /* successful return */
 
}  /* END OF READ_BOREAS */

/*=============================================================================
  Function to get latitude and longitude values

  Input parameters:
      infp == FILE pointer to input file of coordinates
      outfp == FILE pointer to output file of coordinates

  Output parameters:
      lat == latitude coordinate in decimal degrees (+ N and - S)
      lon ==longitude coordinate in decimal degrees (+ E and - W)
=============================================================================*/
long int read_ll (infp, outfp, lat, lon)
FILE *infp, *outfp;
double *lat, *lon;
{
  long int cnt;
  
  char input_str[80];
  
  /*
    The following do while loop allows reprompting of the
    user for erroneous values and handles returns to calling
    function if the input is a disk file
  */
  cnt = 0;
  do
  {
    /*
      Get an input record from the specified file or user and
      extract the needed latitude and longitude values
    */
    if (infp != (FILE *) NULL)
    {
      if (fgets (input_str, 80, infp) == NULL)
      {
        if ( ferror (infp) ) puts ("*** Error reading from input file ***");
        else if ( feof (infp) ) puts ("\n\nEnd of input file reached");
        return (0);
      }
    }
    else
    {
      puts ("\n\nEnter longitude latitude in decimal degrees (xx.xx yy.yy)");
      puts ("(-111.00 <=lon<= -93.00 and 50.75 <=lat<= 60.00); q = quit");
      printf ("( --> ");
      if ( gets (input_str) == NULL) 
      {
        puts ("*** Error reading from terminal ***");
   		continue;
      }
      if ( !check_response (input_str) ) return (0);
    }

    /*
      See if anything was given in the input string
    */
    if ( strlen(input_str) == 0)
    {
      puts ("No coordinates found on input record");
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }
    /*
      Since it seems like something is in the input string
      try to extract the expected coordinates from it
    */
    cnt = sscanf (input_str, "%lf %lf", lon, lat);
    if (cnt != 2)
    {
      puts ("*** Error extracting lon, lat from input record ***");
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }
  
    /*
      Check the latitude value for validity
    */
    if (*lat > 60.00 || *lat < 50.75)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified latitude coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only latitude values ***\n", outfp);
        fputs ("*** between 50.75 and 60.00 are accepted     ***\n", outfp);
        puts ("\n*** Specified latitude coordinate is outside     ***");
        puts ("*** the anticipated range.  Only latitude values ***");
        puts ("*** between 50.75 and 60.00 are accepted         ***");
        return (0);
      }
      else
      {
        puts ("\n*** Specified latitude coordinate is outside     ***");
        puts ("*** the anticipated range.  Only latitude values ***");
        puts ("*** between 50.75 and 60.00 are accepted         ***");
        cnt = 0;
      }
    }
    /*    
	  Check the longitude value for validity
    */
  if (*lon > -93.0 || *lon < -111.0)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified longitude coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only longitude values ***\n", outfp);
        fputs ("*** between -93.00 and -111.00 are accepted   ***\n", outfp);
        puts ("\n*** Specified longitude coordinate is outside  ***");
        puts ("***  anticipated range.  Only longitude values ***");
        puts ("*** between -93.00 and -111.00 are accepted    ***");
        return (0);
      }
      else
      {
        puts ("\n*** Specified longitude coordinate is outside  ***");
        puts ("***  anticipated range.  Only longitude values ***");
        puts ("*** between -93.00 and -111.00 are accepted    ***");
        cnt = 0;
      }
    }

  } while (cnt == 0);  /* End of do while loop */

  return (1);  /* successful return */
  
}  /* END OF READ_LL */

/*=============================================================================
  Function to get UTM northing, easting, and zone values

  Input parameters:
      infp == FILE pointer to input file of coordinates
      outfp == FILE pointer to output file of coordinates

  Output parameters:
      northing == UTM northing in decimal meters
      easting == UTM easting in decimal meters
      zone == UTM zone
=============================================================================*/
long int read_utm (infp, outfp, northing, easting, zone)
FILE *infp, *outfp;
double *northing, *easting;
long int *zone;
{
  long int cnt;
  
  char input_str[80];
  
  /*
    The following do while loop allows reprompting of the
    user for erroneous values and handles returns to calling
    function if the input is a disk file
  */
  cnt = 0;
  do
  {
    /*
      Get an input record from the specified file or user and
      extract the needed northing, easting, and zone values
    */
    if (infp != (FILE *) NULL)
    {
      if (fgets (input_str, 80, infp) == NULL)
      {
        if ( ferror (infp) ) puts ("*** Error reading from input file ***");
        else if ( feof (infp) ) puts ("\n\nEnd of input file reached");
        return (0);
      }
    }
    else
    {
      puts ("\n\nEnter: Easting  Northing  Zone  (q = quit)");
      printf ("\r --> ");
      if (gets (input_str) == NULL)
      {
        puts ("*** Error reading from terminal ***");
        continue;
      }
      if ( !check_response (input_str) ) return (0);
    }

    /*
      See if anything was given in the input string
    */
    if ( strlen(input_str) == 0)
    {
      puts ("No coordinates found on input record");
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }
    /*
      Since it seems like something is in the input string
      try to extract the expected coordinates from it
    */
    cnt = sscanf (input_str, "%lf %lf %ld", easting, northing, zone);
    if (cnt != 3)
    {
      puts("*** Error extracting  easting, northing, zone from input ***");
      if (infp != (FILE *) NULL) return (0);
      else continue;
    }

  
    /*
      Check the easting value for validity
    */
    if (*easting > 710496 || *easting < 289504)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified easting coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only easting values ***\n", outfp);
        fputs ("*** between 289504 and 710496 are accepted  ***\n", outfp);
        puts ("\n*** Specified easting coordinate is outside ***");
        puts ("*** anticipated range.  Only easting values ***");
        puts ("*** between 289504 and 710496 are accepted  ***");
      }
      else
      {
        puts ("\n*** Specified easting coordinate is outside ***");
        puts ("*** anticipated range.  Only easting values ***");
        puts ("*** between 289504 and 710496 are accepted  ***");
        cnt = 0;
      }
    }
    /*
      Check the northing value for validity
    */
    if (*northing > 6655000 || *northing < 5649600)
    {
      if (outfp != (FILE *) NULL)
      {
        fputs ("*** Specified northing coordinate is outside ***\n", outfp);
        fputs ("*** anticipated range.  Only northing values ***\n", outfp);
        fputs ("*** between 5649600 and 6655000 are accepted ***\n", outfp);
        puts ("\n*** Specified northing coordinate is outside of ***");
        puts ("*** anticipated range.  Only northing values    ***");
        puts ("*** between 5649600 and 6655000 are accepted    ***");
      }
      else
      {
        puts ("\n*** Specified northing coordinate is outside of ***");
        puts ("*** anticipated range.  Only northing values    ***");
        puts ("*** between 5649600 and 6655000 are accepted    ***");
        cnt = 0;
      }
    }

} while (cnt == 0);  /* End do while loop */
  
  return (1);  /* Successful return */
  
}  /* END OF READ_UTM */

/*=============================================================================
Function to compute the shifts in latitude and longitude
between the NAD27 and NAD83 datums.

  Input parameters:
      fp == FILE pointer to datum shift information file
      indatum == string giving input coordinate datum
      ncoors == number of coordinate pairs for which to calculate shifts
      inlat == array of input latitudes to be transformed
      inlong == array of input longitudes to be transformed 

  Output parameters:
      olat == array of shifted output latitudes
      olong == array of shifted output longitudes

Original program written  by Marc Veronneau          January 5, 1989
                          Modified by Mario Berube   May 9, 1989
                          Modified by Steve Farley   January 18 1990
                          Converted to C and modified
                          for BOREAS by Fred Irani, Hughes STX,   Feb 1994     

  This function is based on the FORTRAN-77 subroutine developed for the:
  National Transformation (for converting between NAD27 and NAD83 in Canada)
  by:
  
  		Geodetic Survey of Canada
  		615 Booth Street
  		Ottawa, Canada
  		K1A0E9
  		Phone: (613) 995-44120

=============================================================================*/
long int c_gridint(fp, indatum, ncoors, inlat, inlong, olat, olong)
FILE *fp;
char *indatum;
long int	ncoors;
double *inlat, *inlong, *olat, *olong;
{
  
  char astring[81];	/* Character buffer to read from transform shift file */
  
                         /* limits of current transformation grid: */
  double elong =  89.0,	 /* East longitute */
         llat  =  48.0,  /* Lower Latitude */
         ulat  =  60.0,  /* Upper Latitude */
         wlong = 112.0;  /* West longitude (updated from 110 on 3/23 */

                  /* Transformation grid shift values read from file */
  double alat,    /* First nearest Latitude shift */
         along,   /* First nearest longitude shift */
         blat,    /* Next latitude shift on same line of latitude */
         blong,   /* Next longitude shift on same line of latitude */
         clat,    /* First Nearest Latitude shift on next latitude line */
         clong,   /* First Nearest Longitude shift on next latitude line */
         dlat,    /* Next Latitude shift on upper latitude line */
         dlong,   /* Next Longitude shift on upper latitude line */
         diflat,  /* Calculated Latitude shift for specified lat/long */
         diflong, /* Calculated Longitude shift for specified lat/long */
         dir;     /* Direction of calculation: NAD27->83 or reverse (1, -1) */

  double n,      /* Number of grid shift values to desired latitude  */
         ngrid,  /* Number of grid shift values from south to north  */
         w,      /* Number of grid shift values to desired longitude */
         wgrid;  /* Number of grid shift values from east to west */

  double itemp, jtemp, lon_temp;  /* Temporary calculation values */

  long int i,        /* Target latitude index */
           iread,    /* target record index */
           irow,     /* Target row index */
           iscan,    /* sscanf status */
           j,        /* Target longitude index*/
           jread,    /*  Second read record index */
           num_long, /* Number of values to next longitude */
           pair,     /* Current user-specified lat long pair */
           ibytes,   /* Number of bytes to seek to desired record */
           jbytes;   /* Number of bytes to seek to desired record */
  	
  /*
    Initalize grid increment variables.
  */ 
  ngrid  =  (double) 5.0/60.0;
  wgrid  =  (double) 5.0/60.0;
    
  /*
    Calculate the dimensions of the grid point file
  */
  num_long = (long int) ((wlong - elong) / wgrid) + 1;	/* grid width  */
  
  /*
    Set dir as per datum shift direction.
  */
  if ( (strcmp(indatum, "NAD27") == 0) )
  {
    dir = 1.0;
  }
  else if ( (strcmp(indatum, "NAD83") == 0) )
  {
    dir = -1.0;
  }
  else
  {
    printf ("\n*** Invalid input datum (%s) ***\n", indatum);
    return (0);
  }

  /*
    Process the number of coordinates passed in
    
    Since the software requires a positive longitude value in its
    interpolation, assign the incoming negative longitude value
    to temp place holder variable
  */ 
  for ( pair = 0; pair < ncoors; ++pair)
  { 
 
    lon_temp = -inlong[pair];
    
    /*
    Check that the input lat and long are useable.  Report and return if not.
    ------------------------------------------------------------------------- */
    if (lon_temp < elong)
    {
    	printf("Input or calculated longitude of %lf is east of datum shift boundary.\n",
    	        lon_temp);
    	puts("This condition should not have occurred under the original design of\n");
    	puts("BOR_CORD for BOREAS purposes.\n");
    	return(0);
    }
    if (lon_temp > wlong)
    {
    	printf("Input or calculated longitude of %lf is west of datum shift boundary.\n",
    	        lon_temp);
    	puts("This condition should not have occurred under the original design of\n");
    	puts("BOR_CORD for BOREAS purposes.\n");
       	return(0);
    }
    if (inlat[pair] < llat)
    {
    	printf("Input or calculated latitude of %lf is south of datum shift boundary.\n",
    	        lon_temp);
    	puts("This condition should not have occurred under the original design of\n");
    	puts("BOR_CORD for BOREAS purposes.\n");
       	return(0);
    }
    /*
    For BOREAS, calculated latitude should be only slightly above 60 degrees north.
    If this is the case, treat the input latitude as if it were at 60 degrees 
    and report this to users on final output.
    ------------------------------------------------------------------------------- */
    if (inlat[pair] > ulat)
    {
     printf("Input or calculated latitude of %lf is north of 60.\n", inlat[pair]);
     puts("Latitude shifts will be handled as if at 60 degrees\n");
     puts("Note that datum shift file edge coordinates are interpolated");
     puts("differently than interior coordintes\n");
     it
emp = (double)  60.0 - llat;
    }
    else if (inlat[pair] == ulat)
    {
     printf(
     "Input or calculated latitude of %lf is at 60 degrees north.\n", inlat[pair]);
     puts("Note that datum shift file edge coordinates are interpolated");
     puts("differently than interior coordintes.\n");
     itemp = (double)  inlat[pair] - llat;

    }
    else
    { 
        itemp = (double)  inlat[pair] - llat;
    }
    
    itemp /= ngrid;
    i = (long int) itemp; /* Truncate result */

    jtemp = (double) lon_temp - elong;
    jtemp /= wgrid;
    j = (long int) jtemp; /* Like i, truncate result */
 
    /*
      Calculate where to seek to in the input file of datum shift values
      in order to start searching for the closest 4 grid points to
      compute needed latitude and longitude adjustments
    */
    irow = (long int) i * num_long;
    
   
    iread = irow  + j;  /* iread is: current row + # longs needed */
    jread = iread + 1;  /* jread is: get the next lat/long pair   */

    /*
      Read upper nearest lat long coords for interpolation.
    */ 
    ibytes =  iread * NRECBYTES;
    jbytes =  jread * NRECBYTES;

    fseek(fp, ibytes, 0);
    if (fgets(astring, 80, fp) == NULL)
    {
      if (ferror (fp)) puts ("*** Error reading from datum shift file ***");
      else if (feof (fp)) puts ("\n\nEnd of datum shift file reached");
      return (0);
    }
    iscan = sscanf(astring, "%lf %lf",&alat, &along);
    if ( (iscan != 2) && (inlat[pair] != ulat) )
    {
      puts ("*** Error extracting datum shift values from input 1 ***");
      return (0);
    }
    fseek(fp, jbytes, 0);
    if (fgets(astring, 80, fp) == NULL)
    {
      if (ferror (fp))
      {
       puts ("*** Error reading from datum shift file ***");
      }
      else if ( feof (fp) )
      {  
      	puts ("\n\nEnd of datum shift file reached");
      }
      return (0);
    }
    iscan = sscanf(astring, "%10lf %10lf", &blat, &blong);
    if (iscan != 2)
    {
      puts ("*** Error extracting datum shift values from input 2 ***");
      return (0);
    }
      
    /*
    Handle edge coordinates differently than interior coordinates.
    -------------------------------------------------------------- */
    if (inlat[pair] < 60.0)
    {
      /*
        Read upper nearest lat long coords for interpolation.
      */
      iread += num_long;
      jread += num_long;
      ibytes = iread * NRECBYTES;
      jbytes = jread * NRECBYTES;
      fseek(fp, ibytes, 0);
      if (fgets(astring, 80, fp) == NULL)
      { 
        if (ferror (fp))
        {
        	puts ("*** Error reading from datum shift file ***");
        }
        else if ( (feof (fp) ) && (inlat[pair] != ulat) )
        {  
        	puts ("\n\nEnd of datum shift file reached");
        }
        return (0);
      }
      iscan = sscanf(astring, "%lf %lf", &clat, &clong);
      if (iscan != 2)
      {
        puts ("*** Error extracting datum shift values from input 3 ***");
        return (0);
      }
      fseek(fp, jbytes, 0);
      if (fgets(astring, 80, fp) == NULL)
      {
        if (ferror (fp))
        {
        	puts ("*** Error reading from datum shift file ***");
        }
        else if ( (feof (fp) ) && (inlat[pair] != ulat) )
        {  
        	puts ("\n\nEnd of datum shift file reached");
        }
        return (0);
      }
      iscan = sscanf(astring, "%10lf %10lf", &dlat, &dlong);
      if (iscan != 2)
      {
        puts ("*** Error extracting datum shift values from input 4 ***");
        return (0);
      }
       n = (double) (( inlat[pair] - (llat  + ((i-1) * ngrid))) / ngrid);

    }  /* inlat[pair] < 60 */
    else 
    {
    	/*
    	Interpolate northern edge coordinates without northern values.
    	-------------------------------------------------------------- */
    	clat  = alat;
    	clong = along;
    	dlat  = blat;
    	dlong = blong;
    	n = (double) (( 60.0 - (llat  + ((i-1) * ngrid))) / ngrid);
    }
    
    w = (double) (( lon_temp - (elong + ((j-1) * wgrid))) / wgrid);

    diflat  = (double) (alat  + (clat  - alat ) * n + (blat - alat ) * w + 
                       (alat - blat  - clat  + dlat ) * n * w);

    diflong = (double) (along + (clong - along) * n + (blong - along) * w + 
                       (along - blong - clong + dlong) * n * w);

    /*
      Apply corrections and convert seconds to decimal degrees.
    */
    olat[pair] = (double) (inlat[pair] + (dir * ( diflat  / 3600.0 )));
    olong[pair] = (double) (lon_temp + (dir * ( diflong / 3600.0 )) );
    olong[pair] = -olong[pair];
  }  /* for each pair */

  return (1);

}  /* END OF C_GRIDINT */

/*=============================================================================
  Function to output the coordinate information

  Input parameters:
      outfp == FILE pointer to output file of coordinates
      xgrid == BOREAS X grid coordinate in decimal kilometers
      ygrid == BOREAS Y grid coordinate in decimal kilometers
      lat_27 == latitude coordinate in decimal degrees for NAD27
      lon_27 ==longitude coordinate in decimal degrees for NAD27
      lat_83 == latitude coordinate in decimal degrees for NAD83
      lon_83 ==longitude coordinate in decimal degrees for NAD83
      northing_27 == UTM northing in decimal meters for NAD27
      easting_27 == UTM easting in decimal meters for NAD27
      zone_27 == UTM zone for NAD27
      northing_83 == UTM northing in decimal meters for NAD83
      easting_83 == UTM easting in decimal meters for NAD83
      zone_83 == UTM zone for NAD83

  Output parameters:
      Returned status code.
=============================================================================*/
long int write_coords (outfp, xgrid, ygrid, lat_27, lon_27, lat_83, lon_83,
                  northing_27, easting_27, zone_27, northing_83, easting_83,
                  zone_83)
FILE *outfp;
double xgrid, ygrid, lat_27, lon_27, lat_83, lon_83;
double northing_27, easting_27;
long int zone_27;
double northing_83, easting_83;
long int zone_83;
{
  long int stat;
  
  /*
    Output the coordinate information to the user specified location
  */
  if (outfp == (FILE *) NULL)
  {
    stat=printf("\n\nNAD83: BOREAS X  =  %07.3f    BOREAS Y = %07.3f\n", xgrid, ygrid);
    stat=printf("NAD83  Longitude = %10.5lf  Latitude = %9.5lf\n", lon_83, lat_83);
    stat=printf("NAD83  Easting   = %9.1lf   Northing = %9.1lf  Zone = %2ld\n",
			    easting_83, northing_83, zone_83);
    stat=printf("NAD27  Longitude = %10.5lf  Latitude = %9.5lf\n", lon_27, lat_27);
    stat=printf("NAD27  Easting   = %9.1lf   Northing = %9.1lf  Zone = %2ld\n",
			    easting_27, northing_27, zone_27);
    if (stat < 0) return (0);
  }
  else
  {
    stat = fprintf(outfp, " %8.3lf %8.3lf   %10.5lf %8.5lf",
                            xgrid, ygrid, lon_27, lat_27);
    stat = fprintf(outfp, "     %10.5lf  %7.5lf %10.1lf  %9.1lf   %2.2ld   ",
                           lon_83, lat_83, easting_27, northing_27,  zone_27);
    stat = fprintf(outfp, "%9.1lf  %9.1lf   %2.2ld  \n",
                           easting_83, northing_83,  zone_83);
    if (stat < 0) return (0);
  }
  
  return (1);
  
}  /* END OF WRITE_COORDS */

/*=============================================================================
  Function to Check a user response for a quit and return FALSE
  if user does not wish to continue; TRUE if they do.
=============================================================================*/ 
long int check_response (response)
char	*response;
{
  if ( (strcmp(response, "q") == 0) || (strcmp(response, "Q") == 0) )
  {
    return(FALSE);
  }
  else
  {
    return(TRUE);
  }
     
} /* END OF CHECK_RESPONSE */

/*=============================================================================
  Function to output the software version number banner
=============================================================================*/ 
void version_num()
{
  puts ("\n *============================================*");
  puts   (" *  Version 2.2      BOREAS Coordinate        *");
  puts   (" *                Transformation Software     *");
  puts   (" *============================================*");
} /* END OF VERSION_NUM */

