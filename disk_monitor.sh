#!/bin/bash

LOG_FILE=$1 ; 

while [ 1 -eq 1 ] ; 
	do
		date >> $LOG_FILE
		df >> $LOG_FILE 
		sleep 1s ;
	done ;


