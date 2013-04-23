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
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author fishjord
 */
public class PrimerMatch {

    public static void main(String [] args) throws Exception {
        if(args.length != 2 && args.length != 3) {
            System.err.println("USAGE: PrimerMatch <primer_list | primer_file> <seq_file> [max_dist]");
	    return;
        }

	int maxDist = Integer.MAX_VALUE;
	if(args.length == 3) {
	    maxDist = Integer.valueOf(args[2]);
	}

	List<PatternBitMask64> primers = new ArrayList();
	if(new File(args[0]).exists()) {
	    BufferedReader reader = new BufferedReader(new FileReader(args[0]));
	    String line;

	    while((line = reader.readLine()) != null) {
		line = line.trim();
		if(!line.equals("")) {
		    primers.add(new PatternBitMask64(line, true));
		}
	    }
	    reader.close();
	} else {
	    for(String primer : args[0].split(",")) {
		primers.add(new PatternBitMask64(primer, true));
	    }
	}

	SeqReader seqReader = new SequenceReader(new File(args[1]));
	Sequence seq;

	while((seq = seqReader.readNextSequence()) != null) {
	    for(int index = 0;index < primers.size();index++) {
		PatternBitMask64 primer = primers.get(index);
		BitVector64Result results = BitVector64.process(seq.getSeqString().toCharArray(), primer, maxDist);
		
		for(BitVector64Match result : results.getResults()) {
		    System.out.println(seq.getSeqName() + "\t" + seq.getDesc() + "\t" + index + "\t" + result.getPosition() + "\t" + result.getScore());
		}
	    }
	}
    }

}
