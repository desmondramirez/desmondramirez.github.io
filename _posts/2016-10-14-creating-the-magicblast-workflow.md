---
layout: post
title:  "Scripting a pipeline"
date:   2016-10-10
categories: shell bash programming forloops magicblast BLAST transcriptomics
---

# Working on a script pipeline
This strips the landmark1 sequences from PIA fasta files using **grep_LANDMARK1.sh**:

```bash
#grep_LANDMARK1.sh
#created Sat Oct  8 12:29:27 PDT 2016 by Des Ramirez

#tells the computer what language to expect, basically.
#!/bin/bash

#changing the working directory to the specific folder my files are in.
cd /home/desramirez/PIA_datasets/LIT_1.1/rtrans

#creates a list of files that end in .fas. The '*' means 'wildcard', so anything can be in front of the '.fas'. 'ls' tells the computer to list anything that ends in '.fas'. The quotes tell the computer to make a list of the stuff in between them. I'm not totally clear on what the **'$'** means, but I think it does something like "print out the literal output of 'ls *.fas'," which is the name of each of the files I need. It sort of turns things into variables. See '$i' below
## An important note: it's important to have consistent file names if you're working with lists!
LIST="$(ls *.fas)"

#starts a *for loop*. Here I'm saying that for every item (*i*) in the list I just made above (has to be called with the $LIST here in the for loop)...
for i in $LIST

#... do the things below the word "do".
do

#'grep' grabs any line that contains the string "LANDMARK1", AND the next line after it (specified by '-A1') and shoves both lines into a file named after each item in the list (note that we call 'i' here using the '$'). The files will be consistently named with 'BAITS' as a prefix and '.fasta' as a suffix, with the name of the file in the middle.
grep -A1 "LANDMARK1" $i >> BAITS.$i.fasta

#'mv' moves any file that starts with 'BAITS' to the BAITS folder in another directory. Note that 'mv' **moves** the file. If I just wanted to make a copy of the file in another directory, I would use 'cp'
mv BAITS* /home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS

#here I'm 'closing' the 'do' command that I called earlier in the script, saying that after moving the BAITS file, the for loop is done and the computer can move onto applying the loop to the next item in the list.
done


#Above, the for loop is done, but my script isn't done! I have more folders that I'm extracting lines from, so I repeat the same code module as I used above, just changing the path to reflect the next folder I want to apply this script to.

cd /home/desramirez/PIA_datasets/LIT_1.1/ctrans

LIST="$(ls *.fas)"

for i in $LIST

do

grep -A1 "LANDMARK1" $i >> BAITS.$i.fasta

mv BAITS* /home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS

done

cd /home/desramirez/PIA_datasets/LIT_1.1/prc

LIST="$(ls *.fas)"

for i in $LIST

do

grep -A1 "LANDMARK1" $i >> BAITS.$i.fasta

mv BAITS* /home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS

done


cd /home/desramirez/PIA_datasets/LIT_1.1/rdn

LIST="$(ls *.fas)"

for i in $LIST

do

grep -A1 "LANDMARK1" $i >> BAITS.$i.fasta

mv BAITS* /home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS

done


```
Finally, when the script ends, it has run on all my folders. There are other ways to do this, but this is how I did it. If I had more folders, I might make a script to make a separate script for each folder automatically. I do this later on today for a slightly different purpose!

I need nucleotides to run magicblast, but the fasta files I pulled from only had protein sequences, but:
1. the blast hits should come from molluscs.

2. the hits need to be nucleotide, not amino acid, sequences
So, the next step is killing two birds with one stone. **tblastn** using the landmark proteins for each gene to find nucleotide sequences from (in this case) molluscs. I only need a few blast hits, but more might be better for capturing sequence diversity for mapping. I might need to adjust this, depending on how well magicblast is able to map reads... I might need to get more specific sequences for each of the molluscan classes to do the mapping. We'll see. I'm starting with **r-opsin phototransduction** genes first, since this is the basis of my final chapter. If I have enough time later, I can run the same scripts on other genes of interest.

```bash

#!/bin/bash

cd /home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS

LIST="$(ls *rtrans.fas.fasta)"

for i in $LIST

do

~/.local/bin/tblastn -remote -query $i -out $i.tblastn_results -outfmt 7

awk '{print $2}' $i.tblastn_results > $i.tblastn_seqs.list

done

```

#### Unfortunately, I can't restrict to a specific taxon using the commandline, so I have to just do it manually online :-1:

* I used tblastn online, restricting to mollusc hits (taxid: 6447).
* I allowed 50 hits per bait sequence.
* I downloaded each hit table in csv format, used Excel to remove duplicate hits and saved to a new file
* I used multiple **sed** commands to clean up the gi numbers so I could use batch entrez to download the sequences in fasta format*

```bash
#!/bin/bash

cd /Users/desmondramirez/Box\ Sync/DISSERTATION_WORK/RESEARCH/ancestral_r-opsin_expression_in_mollusc_skin/rtrans

LIST="$(ls *_unique_nt_blast_hits.txt)"

for i in $LIST

do

sed -e 's/\([0-9]\)\|.*/\1/g' $i >> gi.$i

sed -i '' 's/gi\|//g' gi.$i

done

```

Used NCBI batch entrez to download nt fasta files for each gene.
Each file is named 'genename'nt_hits.fasta

#### UPDATE! Batch Entrez was not working properly, so I had to just suck it up and use something else to get the sequences I needed.

I ended up finding and using [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started)

which is a way to use command line to run the perl scripts that are the actual "E-utilities" that NCBI provides. NCBI provides code to install Entrez Direct via command line, and then it runs a script to install any missing dependencies (Perl type stuff), and sets the paths properly so that you can call the commands from anywhere.

```bash

LIST="$(ls gi.*)"

for i in $LIST

do

#I needed to make a variable for each file that was a list of the file contents. Then I could call this variable when I needed it for the -query flag in the esearch command.
MYVAR="$(cat $i)"
esearch -db nucleotide -query "$MYVAR" |
efetch -format fasta >> $i.fasta

done

```

Cluster sequence baits with CD-HIT so that there's no so much to blast using **cd-hit_local.sh**

```bash

cd /Users/desmondramirez/Box\ Sync/DISSERTATION_WORK/RESEARCH/ancestral_r-opsin_expression_in_mollusc_skin/rtrans

LIST="$(ls *unique_nt_blast_hits.txt.fasta)"

for i in $LIST

do

cd-hit -i $i -o ntBAITS.$i -c 0.75 -n 5 -M 16000 -d 0 -T 8

done


```
To get everything properly named without having the names be super long, I can use this for loop to rename all of the final files from all of this manipulation to something simple to use in the code later. Then I move each of the files to a folder

```bash

for i in *.gi.*
do
    mv "$i" "${i/.gi./.}"
done

for i in ntBAITS*_unique_nt_blast_hits.txt.fasta
do
    mv "$i" "${i/_unique_nt_blast_hits.txt.fasta/.fasta}"
done

mv ntBAITS*.fasta rtrans_magicblast_dbs

```


The last step before magicblast is to make a blast database for the nt fasta files I have for each gene, which I did using **makeblastdb.sh**

```bash
#!/bin/bash

#This path will need to change depending on where the script is being run. Since I was having to access batch entrez by hand, I was doing that locally in my computer to avoid having to shuffle too many things between the cluster and my local computer.
#Once the databases are created, I can upload them to the cluster and do the rest of the pipeline there. This stuff is really just getting all the pieces in order before running the magicblast pipeline.
#cd /Users/desmondramirez/Box\ Sync/DISSERTATION_WORK/RESEARCH/ancestral_r-opsin_expression_in_mollusc_skin/rtrans/rtrans_magicblast_dbs

LIST="$(ls *.fasta)"

for i in $LIST

do

makeblastdb -in $i -dbtype nucl -parse_seqids -out $i.blastdb

done

for i in $LIST

do

mkdir "$i"_blastdb
mv $i*.n* "$i"_blastdb

done

```

I uploaded all the baits blast databases using scp

On my remote server, I created a couple of variables for the session, since I need to be bouncing around between folders:

```bash

#where my blast databases are stored
DBS="/home/desramirez/data/ancestral_r-opsin_mollusc_skin/BAITS/rtrans/rtrans_magicblast_dbs"

```

I want to make symbolic links to all the datasets that I've downloaded. They are contained in multiple folders, so I'm going to visit each folder and make the symlinks, but direct all the symlinks to a single folder. In that folder I'll subdivide the links a bit further, as most of the SRAs are paired .fastq files, but a few are older sequencing files that are only singles. The magicblast commands are slightly different depending on whether the data are paired end or not.

```bash
LIST="$(ls *.fastq)"

for i in $LIST

do

ln -s ./$i ~/data/ancestral_r-opsin_mollusc_skin/transcriptomes/SRAs/paired/$i

done

LIST="$(ls *.fastq.gz)"

for i in $LIST

do

ln -s /home/desramirez/data/SRAs_from_EBI/single/$i ~/data/ancestral_r-opsin_mollusc_skin/transcriptomes/SRAs/single/$i

done

```

I want to make a magicblast script for each gene, since each gene has it's own blast database

```bash

```

Now I can run magicblast on each SRA for each r-opsin phototransduction gene blast_db

```bash

LIST="$(ls .fastq.gz)"

for i in $LIST

do

#magicblast -query  -query_mate  -db  -infmt fastq -outfmt tabular -num_threads 10

magicblast -query $i -db  -infmt fastq -outfmt tabular -num_threads 10

I made symlinks for both the Aculifera and Bivalvia assemblies for easy reference

```bash
LIST="$(ls *.fa)"

for i in $LIST

do

ln -s /home/desramirez/data/Bivalvia_Assemblies_Gonzalez/$i /home/desramirez/data/ancestral_r-opsin_mollusc_skin/transcriptomes/assemblies/$i

done

```

There was a long cue with all these jobs I had running on the cluster, plus I had met my user cap for the number of jobs I could have going at once. Luckily, I had both the bivalve and Aculifera assemblies on my local computer, so I made blastdbs locally, and then blasted locally as well.

the script is called "local_makeblastdb.sh"

```bash

#!/bin/sh

LIST="$(ls *.fa)"

for i in $LIST

do

makeblastdb -in $i -dbtype nucl -parse_seqids -out $i.nt.blastdb

done

```


Now I want to blast the gene.fasta files against these blastdbs to pull out good blast hits for each gene from each species.

``` bash

LIST="Arr
Gq_alpha
PLC
TRP
DAGK_rdgA
Gq_beta
rdgB
Gprk1
Gq_gamma
rdgC
Gprk2
PKC
ropsin"

for i in $LIST

do

  sed -e "s/template/$i/g" template_sblastn_only.sh | grep blastn >> all_rtrans_sblastn_only_serial.sh

done


```
This put the **sed** commands in to the file called "all_rtrans_sblastn_only_serial.sh", but I had to add the wrapper of the list of the names of the blast databases for each transcriptome. I also made a job on the cluster to submit this script, but I realized just now that I could have just put the job info into the "all rtrans..." script.


The output from blast is a bit annoying-- it's a table with information about each hit in terms of % identity, etc. What I really need for the next step are the sequence IDs. It's a tab delimited file, so I can use **awk** to grab just them!

```bash

for i in $(ls SRR*); do awk -F'\t' '{print $1}' $i > "$i"_IDs.txt; done &


for i in $(ls *blast_results); do awk -F'\t' '{print $2}' $i > "$i"_IDs.txt; done &
```
##### It takes a while for these scripts to run... it's a lot of searching!! I actually made it executable so that the cluster could just run it, on the head node.
