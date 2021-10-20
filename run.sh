#!/bin/bash
nextflow run main.nf\
	-with-tower\
	-c configs/institutional/uppmax.config\
	-c configs/conf/upprun.config
