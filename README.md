AMPA-NF
=======

AMPA-NF is a pipeline for assessing the antimicrobial domains of proteins, 
with a focus on the design on new antimicrobial drugs. The application provides fast discovery of 
antimicrobial patterns in proteins that can be used to develop new peptide-based drugs against pathogens. 


Quick start 
-----------

Make sure you have installed the dependencies required listed at the end of this document. 

Install the Nextflow runtime by running the following command: 

    $ curl -fsSL get.nextflow.io | bash


When done, you can execute AMPA by entering the command shown below:

    $ ./nextflow run cbcrg/ampa-nf


By default the pipeline is executed by using a small dataset included with the project. Check the *Pipeline parameters* section below to see how specify your input data on the program command line.


Pipeline parameters
-------------------

#####--in

  * The query sequences multi-fasta file (default: example.fa)
  * Example: `nextflow ampa.nf --in /some/path/your-query.fa`
  
  
#####--out

  * Path to the file where save the results (default: bigampa.data)
  * Example: `nextflow ampa.nf --out /some/path/your-result.data`
  
  
#####--t

  * Threshold value (default: 0.225)
  * Example: `nextflow ampa.nf --t 0.3`  


#####--w

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
