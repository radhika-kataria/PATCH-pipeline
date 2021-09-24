#!/bin/bash
nextflow run main.nf\
	-with-tower\
	-c /home/theo/postdoc/rhadikaproj/PATCH-pipeline/configs/institutional/laptop.config\
	-c /home/theo/postdoc/rhadikaproj/PATCH-pipeline/configs/conf/laprun.config
