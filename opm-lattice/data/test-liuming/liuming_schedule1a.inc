RPTSCHED
'PRES' 'SOIL' 'SWAT' 'RESTART=2' 'FIP=1' 'WELLS=1' 'SUMMARY=1' /
  
RPTSCHED
'RESTART=2' 'FIP=2' 'WELLS=2' 'VFPPROD=1' 'CPU=2' 'WELSPECS' 'FIPTR=2' 
/
  
WELSPECS
--Well specification data
 'PROD1 '  'G'   5  10  12000  'OIL' /
/

COMPDAT
--Define connection data
 'PROD1 '  2*    1  10  'OPEN'  2*  1.0  1*  /
/

WCONPROD
'PROD1'  OPEN  BHP  1*  1*  1*  1*  1*  2612.5  /
/

TUNING
0.1 5 8* / 
11* / 
2* 80 30 6* /

TSTEP
 50 50 50 50          -- to day 200
 50 50 50 50 50 50
 50 50 50 50 50 50
 50 50 50 50 50 50
 50 50 50 50 50 50
 50 50 50 50 50 50
 50 50 50 50 50 50 /  -- to day 2000
 