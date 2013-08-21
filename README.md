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
        USAGE: PrimerMatch <primer_list | primer_file> <seq_file> [max_dist]
        
	If multiple primers are used, you can either save the primers in a text file with one primer per line, or enter the primers as command-line argument, separated by "," without any space betwee them. 
	An example command using the sequence data and primer from RDP tutorial, [see link above](#Tutorial).

        java -jar /path/to/ProbeMatch.jar AYTGGGYDTAAAGNG 1.TCA.454Reads.fna 2 > probematch_out.txt
	The output is a tab delimited text file containing sequeneces that match at least one of the primer(s) within the specified distance. 
    Each line contains the seqname, description, the index of closest primer, the starting position of the bases after primer, 
    the distance between the closest primer and the sequence. Example output from the above command:
                
		HC9DO0P01APXU0  rank=0000156 x=178.0 y=1306.0 length=375        0       25      0
		HC9DO0P01AVGVN  rank=0000166 x=241.0 y=1185.0 length=72 0       24      1
		HC9DO0P01BB1KW  rank=0000271 x=430.0 y=366.5 length=372 0       25      0
