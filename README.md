AMPA-NF
=======

AMPA-NF is a pipeline for assessing the antimicrobial domains of proteins, 
with a focus on the design on new antimicrobial drugs. The application provides fast discovery of 
antimicrobial patterns in proteins that can be used to develop new peptide-based drugs against pathogens. 


Quick start 
-----------

Clone the git repository on your computer with the following command:

    $ git clone ggit@github.com:cbcrg/ampa-nf.git
    

Make sure you have installed the dependencies required listed below. 


When done, move in the project folder named `ampa-nf`, 
which contains an example dataset file named `example.fa`. 

Launch the pipeline by entering the following command 
on your shell terminal:

    $ ./nextflow ampa.nf
    

By default the pipeline is executed against the provided tutorial dataset and the result is saved into
a file named `bigampa.data` in the execution folder.

Check the *Pipeline parameters*  section below to see how enter your data on the program command line.    


Pipeline parameters
-------------------

#####in

  * The query sequences multi-fasta file (default: example.fa)
  * Example: `nextflow ampa.nf --in /some/path/your-query.fa`
  
  
#####out

  * Path to the file where save the results (default: bigampa.data)
  * Example: `nextflow ampa.nf --out /some/path/your-result.data`
  
  
#####t

  * Threshold value (default: 0.225)
  * Example: `nextflow ampa.nf --t 0.3`  


#####w

  * Window size value (default: 7)
  * Example: `nextflow ampa.nf --w 8`   


Dependencies 
------------

- PERL
- Math::CDF
- Math::Round 
- Bio::SeqIO


Cite
----

[AMPA: an automated web server for prediction of protein antimicrobial regions](http://bioinformatics.oxfordjournals.org/content/28/1/130.long)