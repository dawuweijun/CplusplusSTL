--RPTSCHED
--'VWAT=1' /

--RPTSCHED
--'PRES' 'SOIL' 'SWAT' 'RESTART=2' 'FIP=1' 'WELLS=1' 'SUMMARY=1' /

WELSPECS
--'INJE1'   'G'   1   1     1*    'GAS'   /
'INJE1'   'G'   1   1     1*    'WATER'   /
'PROD1'   'G'   10  10    1*    'OIL'   /
/

COMPDAT
'INJE1'   2*   1   1     'OPEN'    2*   0.5   2*   /
'PROD1'   2*   3   3     'OPEN'    2*   0.5   2*   /
/

WCONINJE
--'INJE*'   'GAS'   'OPEN'   'RATE'   100000.0   1*   10000    /
'INJE*'   'WATER'   'OPEN'   'RATE'   100000.0   1*   10000    /
/

WCONPROD
'PROD*'   'OPEN'    'ORAT'   20000.0   4*   1000    /
/
