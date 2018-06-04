# SLR
SLR is a tool for scaffolding using long reads.
License
=========

Copyright (C) 2014 Junwei Luo(luojunwei@hpu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Junwei Luo(luojunwei@hpu.edu.cn)
College of Computer Science and Technology
Henan Polytechnic University
Jiaozuo
454000
China


Scaffolder: SLR
=================

1) Introduction

	SLR is an scaffolder which aims to determine the orientations and orders of contigs. 
	The contigs can be produced by any assembler.
	The input data of SLR is the long reads (fasta format) and the contigs.

2) Before installing and running

	First, Please build and install Samtools and Bamtools. And add enviroment vairable BAMTOOLS_HOME which is the path of bamtools.

3) Installing.

	SLR should run on Linux operating sysetm with gcc. We test SLR using gcc4.6.3 on Ubuntu.
	Create a main directory (eg:SLR). Copy all source code to this directory.
	cd SLR
	export BAMTOOLS_HOME=/path_bamtools/
	make

4) Running.

	Run command line: 
  bwa index contigs.fasta
  bwa mem -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y contigs.fasta longreads.fasta > aligning.sam
  samtools view -Sb aligning.sam > aligning.bam
	SLR -c <contigs.fa> -r <aligning.bam> -p <output_directory> -m <minimum_contig_length> -n <minimum_read_length>
	<contigs.fa>: 
		The file includes contigs produced by one assembler.
	<aligning.bam>:
		The aligning result between the contigs and the long reads.
	<output_directory>:
		The directory of scaffolding result.
	<minimum_contig_length>: 
		The contigs whose lengthes are larger than minimum_contig_length will be used for scaffolding.
	<minimum_read_length>: 
		The long reads whose lengthes are larger than minimum_read_length will be used for scaffolding.
