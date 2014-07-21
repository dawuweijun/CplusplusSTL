SCHEDULE
WELSPECS
--Well specification data
-----name----group--I---J---depth---phase
     'INJE1'   'G'  24 25  9110  'WATER'/
--	 'PROD2'   'G'   5  1  9110  'OIL'  /
	 'PROD3'   'G'   8  2  9110  'OIL'  /
	 'PROD4'   'G'  11  3  9110  'OIL'  /
	 'PROD5'   'G'  10  4  9110  'OIL'  /
	 'PROD6'   'G'  12  5  9110  'OIL'  /
	 'PROD7'   'G'   4  6  9110  'OIL'  /
	 'PROD8'   'G'   8  7  9110  'OIL'  /
	 'PROD9'   'G'  14  8  9110  'OIL'  /
	 'PROD10'  'G'  11  9  9110  'OIL'  /
	 'PROD11'  'G'  12 10  9110  'OIL'  /
	 'PROD12'  'G'  10 11  9110  'OIL'  /
	 'PROD13'  'G'   5 12  9110  'OIL'  /
	 'PROD14'  'G'   8 13  9110  'OIL'  /
	 'PROD15'  'G'  11 14  9110  'OIL'  /
	 'PROD16'  'G'  13 15  9110  'OIL'  /
	 'PROD17'  'G'  15 16  9110  'OIL'  /
	 'PROD18'  'G'  11 17  9110  'OIL'  /
	 'PROD19'  'G'  12 18  9110  'OIL'  /
	 'PROD20'  'G'   5 19  9110  'OIL'  /
	 'PROD21'  'G'   8 20  9110  'OIL'  /
	 'PROD22'  'G'  11 21  9110  'OIL'  /
	 'PROD23'  'G'  15 22  9110  'OIL'  /
	 'PROD24'  'G'  12 23  9110  'OIL'  /
	 'PROD25'  'G'  10 24  9110  'OIL'  /
	 'PROD26'  'G'  17 25  9110  'OIL'  /
/

COMPDAT
--Define connection data
----------I--J---K1--K2---status
 'INJE1'  24 25 11  15  'OPEN'  2*  1.0  1*  /
-- 'PROD2'   5  1  2   4  'OPEN'  2*  1.0  1*  /
 'PROD3'   8  2  2   4  'OPEN'  2*  1.0  1*  /
 'PROD4'  11  3  2   4  'OPEN'  2*  1.0  1*  /
 'PROD5'  10  4  2   4  'OPEN'  2*  1.0  1*  /
 'PROD6'  12  5  2   4  'OPEN'  2*  1.0  1*  /
 'PROD7'   4  6  2   4  'OPEN'  2*  1.0  1*  /
 'PROD8'   8  7  2   4  'OPEN'  2*  1.0  1*  /
 'PROD9'  14  8  2   4  'OPEN'  2*  1.0  1*  /
 'PROD10' 11  9   2   4  'OPEN'  2*  1.0  1*  /
 'PROD11' 12 10   2   4  'OPEN'  2*  1.0  1*  /
 'PROD12' 10 11   2   4  'OPEN'  2*  1.0  1*  /
 'PROD13'  5 12   2   4  'OPEN'  2*  1.0  1*  /
 'PROD14'  8 13   2   4  'OPEN'  2*  1.0  1*  /
 'PROD15' 11 14   2   4  'OPEN'  2*  1.0  1*  /
 'PROD16' 13 15   2   4  'OPEN'  2*  1.0  1*  /
 'PROD17' 15 16   2   4  'OPEN'  2*  1.0  1*  /
 'PROD18' 11 17   2   4  'OPEN'  2*  1.0  1*  /
 'PROD19' 12 18   2   4  'OPEN'  2*  1.0  1*  /
 'PROD20'  5 19   2   4  'OPEN'  2*  1.0  1*  /
 'PROD21'  8 20   2   4  'OPEN'  2*  1.0  1*  /
 'PROD22' 11 21   2   4  'OPEN'  2*  1.0  1*  /
 'PROD23' 15 22   2   4  'OPEN'  2*  1.0  1*  /
 'PROD24' 12 23   2   4  'OPEN'  2*  1.0  1*  /
 'PROD25' 10 24   2   4  'OPEN'  2*  1.0  1*  /
 'PROD26' 17 25   2   4  'OPEN'  2*  1.0  1*  /
/

WCONINJE
--                              Max        BHP
--                              Rate       limit
 'INJE1' 'WATER' 'OPEN' 'RATE'  5000  1*   4000   /
/

WCONPROD
--                     Max       BHP
--                     Rate      limit
-- 'PROD2' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD3' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD4' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD5' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD6' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD7' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD8' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD9' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD10' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD11' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD12' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD13' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD14' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD15' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD16' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD17' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD18' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD19' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD20' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD21' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD22' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD23' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD24' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD25' 'OPEN' 'ORAT' 1500  4*  1000   /
 'PROD26' 'OPEN' 'ORAT' 1500  4*  1000   /
/


TSTEP
-- to time 100
 1 49 50
-- to time 300
 2*100
-- to time 360
/

WELTARG
--	'PROD*'  'O'  100.0 /
-- 'PROD2' 'O'  100.0 /
 'PROD3' 'O'  100.0 /
 'PROD4' 'O'  100.0 /
 'PROD5' 'O'  100.0 /
 'PROD6' 'O'  100.0 /
 'PROD7'  'O'  100.0 /
 'PROD8'  'O'  100.0 /
 'PROD9'  'O'  100.0 /
 'PROD10' 'O'  100.0 /
 'PROD11' 'O'  100.0 /
 'PROD12'  'O'  100.0 /
 'PROD13'  'O'  100.0 /
 'PROD14'  'O'  100.0 /
 'PROD15'  'O'  100.0 /
 'PROD16'  'O'  100.0 /
 'PROD17'  'O'  100.0 /
 'PROD18'  'O'  100.0 /
 'PROD19'  'O'  100.0 /
 'PROD20'  'O'  100.0 /
 'PROD21'  'O'  100.0 /
 'PROD22'  'O'  100.0 /
 'PROD23'  'O'  100.0 /
 'PROD24'  'O'  100.0 /
 'PROD25'  'O'  100.0 /
 'PROD26'  'O'  100.0 /
/

TSTEP
 60
/

WELTARG
--	'PROD*'  'O' 1500.0 /
-- 'PROD2' 'O'  1500.0 /
 'PROD3' 'O'  1500.0 /
 'PROD4' 'O'  1500.0 /
 'PROD5' 'O'  1500.0 /
 'PROD6' 'O'  1500.0 /
 'PROD7'  'O'  1500.0 /
 'PROD8'  'O'  1500.0 /
 'PROD9'  'O'  1500.0 /
 'PROD10' 'O'  1500.0 /
 'PROD11' 'O'  1500.0 /
 'PROD12'  'O'  1500.0 /
 'PROD13'  'O'  1500.0 /
 'PROD14'  'O'  1500.0 /
 'PROD15'  'O'  1500.0 /
 'PROD16'  'O'  1500.0 /
 'PROD17'  'O'  1500.0 /
 'PROD18'  'O'  1500.0 /
 'PROD19'  'O'  1500.0 /
 'PROD20'  'O'  1500.0 /
 'PROD21'  'O'  1500.0 /
 'PROD22'  'O'  1500.0 /
 'PROD23'  'O'  1500.0 /
 'PROD24'  'O'  1500.0 /
 'PROD25'  'O'  1500.0 /
 'PROD26'  'O'  1500.0 /
/
-- to time 420
TSTEP
 5*60
-- to time 720
 1*60 
-- to time 900
 180 
/

