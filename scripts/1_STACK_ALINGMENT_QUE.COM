$!      
$!         _/_/_/  _/_/_/    _/        
$!          _/    _/    _/  _/           Image Processing Language 
$!         _/    _/_/_/    _/ 
$!        _/    _/        _/             (c) Andres Laib, IBT/ETH Zürich
$!     _/_/_/  _/        _/_/_/_/            Scanco Medical
$!        
$! IPL Batch Scanco
$!========================================================================================
$! This script submits scans to stack alignment algorithm
$
$! INPUTS:
$!              a. The grayscale ".AIM" file
$
$! OUTPUTS:
$!              A. Aligned AIM and GOBJ files
$!              B. A logfile in the scan folder
$!
$!====================================================================
$! call as:
$! @1_STACK_ALINGMENT_QUE.COM   C00001.AIM
$!====================================================================
$ FILENAME 	 = F$PARSE("''P1'",,,"NAME") - ""
$ DEV      		= F$PARSE("''P1'",,,"DEVICE")
$ DIR_READ   	= F$PARSE("''P1'",,,"DIRECTORY")
$ DIR1   		:== 'DEV''DIR_READ''FILENAME'
$ ! Identify the user account
$ STRIPPED = F$EXTRACT(1, F$LENGTH(DEV) - 2, DEV)  ! remove the [ and ]
$ USER = F$ELEMENT(0, ".", STRIPPED)
$
$ QUEUE_NAME :== SYS$FAST
$ LOG_FILE :== 'DEV''DIR_READ''FILENAME'_STACK_ALIGNMENT.LOG
$ SCRIPT = "DISK1:[" + USER + "]2_STACK_ALINGMENT_IPL_V1.COM"
$ INPUT1 = P1
$ SUBMIT/QUEUE='QUEUE_NAME'/PARAM=('INPUT1')/LOG=('LOG_FILE') 'SCRIPT'
$
$ EXIT