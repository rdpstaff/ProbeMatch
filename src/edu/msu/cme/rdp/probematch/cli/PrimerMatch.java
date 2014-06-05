/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.msu.cme.rdp.probematch.cli;

import edu.msu.cme.rdp.probematch.myers99.BitVector64;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result.BitVector64Match;
import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;
import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class PrimerMatch {

    private static final Options options = new Options();

    static {
        options.addOption("n", "maxDist", true, "Give a max distance");
        options.addOption("o", "outFile", true, "Write output to a file");
    }
    
    
    public static void main(String [] args) throws Exception {
 
        PrintStream out = new PrintStream(System.out);
        int maxDist = Integer.MAX_VALUE;
        
        try {
            CommandLine line = new PosixParser().parse(options, args);
            if (line.hasOption("outFile")) {
                out = new PrintStream(new File(line.getOptionValue("outFile")));
            }
            if (line.hasOption("maxDist")){
                maxDist = Integer.valueOf(line.getOptionValue("maxDist"));
            }
            args = line.getArgs();
            
            if(args.length != 2){
                throw new Exception("Unexpected number of command line arguments");
            }
        } catch (Exception e){
            System.err.println("Error: " + e.getMessage());
            new HelpFormatter().printHelp("PrimerMatch <primer_list | primer_file> <seq_file>", options);
            return;
        }
        
	List<PatternBitMask64> primers = new ArrayList();
	if(new File(args[0]).exists()) {
            File primerFile= new File(args[0]);
            SequenceFormat seqformat = SeqUtils.guessFileFormat(primerFile);
	    
            if(seqformat.equals(SequenceFormat.FASTA)){
                SequenceReader reader = new SequenceReader(primerFile);
                Sequence seq;
                
                while ( (seq=reader.readNextSequence()) != null){
                    primers.add(new PatternBitMask64(seq.getSeqString(), true, seq.getSeqName()));
                }
                reader.close();
            } else{
                BufferedReader reader = new BufferedReader(new FileReader(args[0]));
                String line;
                
                while((line = reader.readLine()) != null) {
                    line = line.trim();
                    if(!line.equals("")) {
                        primers.add(new PatternBitMask64(line, true));
                    }
                }
                reader.close();
            } 
	} else {
	    for(String primer : args[0].split(",")) {
		primers.add(new PatternBitMask64(primer, true));
	    }
	}

	SeqReader seqReader = new SequenceReader(new File(args[1]));
	Sequence seq;
        String primerRegion;

        out.println("#seqname\tdesc\tprimer_index\tprimer_name\tposition\tmismatches\tseq_primer_region");
	while((seq = seqReader.readNextSequence()) != null) {
	    for(int index = 0;index < primers.size();index++) {
		PatternBitMask64 primer = primers.get(index);
		BitVector64Result results = BitVector64.process(seq.getSeqString().toCharArray(), primer, maxDist);

		for(BitVector64Match result : results.getResults()) {
                    primerRegion = seq.getSeqString().substring(Math.max(0, result.getPosition() - primer.getPatternLength()), result.getPosition());

                    if(result.getPosition() < primer.getPatternLength()) {
                        for(int pad = result.getPosition(); pad < primer.getPatternLength();pad++) {
                            primerRegion = "x" + primerRegion;
                        }
                    }

		    out.println(seq.getSeqName() + "\t" + seq.getDesc() + "\t" + (index + 1) + "\t" + primer.getPrimerName() + "\t" + result.getPosition() + "\t" + result.getScore() + "\t" + primerRegion);
		}
	    }
	}
        out.close();
        seqReader.close();
    }

}
