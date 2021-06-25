# rnas_and_parasites_origin_of_life
Source code for the article: Parasites sustain and enhance RNA-like replicators through spatial self-organisation

###############################################################################################

This is the source code used to produce all results of the article:
Parasites Sustain and Enhance RNA-like Replicators Through Spatial Self-Organisation
by Enrico Sandro Colizzi and Paulien Hogeweg, published in PLoS Computational Biology.

It simulates an RNA-like replicator-parasite system on a lattice (as described in the publication).

For all I care, you can take this code and do with it anything you want (citing us if you do is appreciated, of course), 
but you should probably check PLoS license.

The code is written in c, and makes use of the CASH libraries, 
which implement Cellular Automata dynamics, and are also freely available at the link:
http://bioinformatics.bio.uu.nl/rdb/software.html

I compile it with gcc under Ubuntu, and should compile just as well under any Linux system. 
I don't know if this compiles with other OS... in fact, I have no idea how to compile anything under any other OS :P

The source code contains a parameter file, parameters.c, which has all the parameters that can be tuned to obtain the various results.
- If you want mutations not to affect parasites, you'll have to change it in replication.c
- For the ablation experiments you have to use the function Meteorites()
- For the invasion experiment you'll have to set boundary conditions to FIXED

The program generates three types of outputs (see the file graphical_output.c):
- "dataMut.txt" is a text file which contains periodic dumps of the populations on lattice in random order
  **-- BEWARE --** This file can get VERY big!
- "movieMut" is a directory which contains a collection of numbered png files that can be watched subsequently very fast to produce the illusion of a movie 8O
- "backupMut" is a directory which collects text files each containing a whole dump of the information of the lattice in ordered fashion.

Do write me if you have questions, you discover bugs, or if you just like the pretty videos you can make with it!

Have fun!
Enrico Sandro Colizzi (email: istidina -at- gmail -dot- com)

###############################################################################################
