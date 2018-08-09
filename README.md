# SLR
SLR is a scaffolding tool based on long reads and contig classification.
=========
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

Junwei Luo(luojunwei@hpu.edu.cn),
College of Computer Science and Technology,
Henan Polytechnic University,
Jiaozuo,
454000,
China


Scaffolder: SLR
=================

1) Introduction
```
    SLR is an scaffolder which aims to determine the orientations and orders of contigs. 
    The contigs can be produced by any assembler.
    The input data of SLR is the long reads (fasta format) and the contigs (fasta format). SLR can classify the contigs into unique contigs and ambiguous contigs.
```
2) Before installing and running
```
    First, Please build and install BWA, Samtools and Bamtools. 
    And add enviroment vairable "BAMTOOLS_HOME_INCLUDE" which is the path that 
    includes two directories "api" and "shared" of bamtools.
    And add enviroment vairable "BAMTOOLS_HOME_LIB" which is the path that 
    includes the file "libbamtools.a".
```
3) Installing.
```
    SLR should run on Linux operating sysetm with gcc. We test SLR using gcc4.6.3 on Ubuntu.
    Create a main directory (eg:SLR). Copy all source code to this directory.
	cd SLR
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```
4) Running.
```
    Step 1: bwa index contigs.fasta
    Step 2: bwa mem -t8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y contigs.fasta longreads.fasta > aligning.sam
    Step 3: samtools view -Sb aligning.sam > aligning.bam
    Step 4: SLR -c <contigs.fa> -r <aligning.bam> -p <output_directory> [-m <minimum_contig_length>] [-n <minimum_read_length>]
	<contigs.fa>: 
	    The file includes contigs produced by one assembler.
	<aligning.bam>:
	    The aligning result between the contigs and the long reads.
	<output_directory>:
	    The output directory of scaffolding result.
	<minimum_contig_length>: 
	    The contigs whose lengths are larger than minimum_contig_length will be used for scaffolding. Default:2000
	<minimum_read_length>: 
	    The long reads whose lengths are larger than minimum_read_length will be used for scaffolding. Default:3000
```
5) Output.
```
    The output file "scaffold_set.fa" is the scaffolding result. Unique contigs are in "unique-contig-set.fa". Ambigous contigs are in "ambiguous-contig-set.fa"
```
6) Running other scaffolding tools based on contig classification
```
    Firstly, users can run other scaffolding tools using unique contig set and long read set, which generate a scaffolding result. Next users can insert ambiguous contigs to the scaffolds.
    Step 1: Running SLR and Getting the contig set file: unique-contig-set.fa.
    Setp 2: Running other scaffolding tools (SSPACE-LR, LINKS) based on unique-contig-set.fa, and Getting scaffolding result: temp-scaffold-set.fa.
    Step 3: bwa index temp-scaffold-set.fa
    Step 4: bwa mem temp-scaffold-set.fa unique-contig-set.fa > unique-contig.sam
    Step 5: samtools view -Sb unique-contig.sam > unique-contig.bam
    Step 6: bwa mem -t8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y contigs.fasta longreads.fasta > aligning.sam
    Step 7: samtools view -Sb aligning.sam > aligning.bam
    Step 8: SLR-unique-ambiguous -c <contigs.fa> -r <aligning.bam> -u <unique-contig-set.fa> -s <temp-scaffold-set.fa> -b <unique-contig.bam> -p <output_directory> 

