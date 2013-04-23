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

package edu.msu.cme.rdp.probematch.oneoff;

import edu.msu.cme.rdp.probematch.myers99.BitVector64;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result;
import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class BondTreeTest {

    private static char[] loadDb() throws IOException {
        final File seqFile = new File("/work/fishjord/other_projects/bond_tree/all_genomes.fasta");
        final int desiredDbSize = 10000000;
        SeqReader seqReader = new SequenceReader(seqFile);
        Sequence seq = null;

        StringBuilder db = new StringBuilder();
        while((seq = seqReader.readNextSequence()) != null) {
            String seqString = seq.getSeqString();
            if(seqString.length() + db.length() >= desiredDbSize) {
                int basesNeeded = desiredDbSize - db.length();

                db.append(seqString.substring(0, basesNeeded));
                break;
            }

            db.append(seqString);
        }

        System.out.println("Desired dbSize= " + db.length() + " desiredSize= " + desiredDbSize);

        return db.toString().toCharArray();
    }

    public static void main(String [] args) throws Exception {
        char[] db = loadDb();

        BufferedReader reader = new BufferedReader(new FileReader("/work/fishjord/other_projects/bond_tree/primers.txt"));
        String line;

        while((line = reader.readLine()) != null) {
            line = line.trim();
            if(line.equals("")) {
                continue;
            }

            PatternBitMask64 primer = new PatternBitMask64(line, true);
            long startTime = System.currentTimeMillis();

            BitVector64Result results = BitVector64.process(db, primer, 0);

            System.out.println(line + "\t" + (System.currentTimeMillis() - startTime) + "\t" + results.getResults().size());
        }
    }

}
