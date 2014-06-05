##RDP Primer/Probe Match tool

### Intro
RDP's ProbeMatch tool check the primer/probe coverage from an input sequence file. See http://rdp.cme.msu.edu/probematch/ for more details.

<a name="Tutorial"></a>
The sample input and primer files can be downloaded from RDP tutorial http://rdp.cme.msu.edu/tutorials/init_process/RDPtutorial_INITIAL-PROCESS.html.

### Setup
Build from command line: ant jar
or see RDPTools (https://github.com/rdpstaff/RDPTools) to install.

### Usage

* Run ProbeMatch

		java -jar /path/to/ProbeMatch.jar
        usage: PrimerMatch <primer_list | primer_file> <seq_file> 
            -n,--maxDist <arg>   Give a max distance
            -o,--outFile <arg>   Write output to a file
        
	If multiple primers are used, you can either:

        * Use plain text file with one primer per line
        * Use fasta formatted file; the primer name will be the sequence id and the primer will be the sequence string.
        * Enter the primers as command-line argument, separated by "," without any space between them. 

	An example command using the sequence data and primer from RDP tutorial, [see link above](#Tutorial).

        java -jar /path/to/ProbeMatch.jar AYTGGGYDTAAAGNG 1.TCA.454Reads.fna -n 2 -o probematch_out.txt
	The output is a tab delimited text file containing sequences that match at least one of the primer(s) within the specified distance. 
    Each line contains the seqname, description, the index of closest primer, the name of the closest primer (if applicable), the starting position of the bases after primer, 
    the distance between the closest primer and the sequence, and the seq_primer_region. Example output from the above command:

                #seqname	desc                                            primer_index	primer_name	position	mismatches	seq_primer_region
                HC9DO0P01APXU0	rank=0000156 x=178.0 y=1306.0 length=375	1	NA              25              0               ATTGGGCATAAAGGG
                HC9DO0P01AVGVN	rank=0000166 x=241.0 y=1185.0 length=72         1	NA              24              1               CATTGGGCGTAAGGG
                HC9DO0P01BB1KW	rank=0000271 x=430.0 y=366.5 length=372         1	NA              25              0               ACTGGGCGTAAAGGG