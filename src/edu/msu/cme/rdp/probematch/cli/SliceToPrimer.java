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

import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import java.util.Set;
import java.util.HashSet;
import edu.msu.cme.rdp.probematch.myers99.BitVector64;
import edu.msu.cme.rdp.probematch.core.DPMAligner;
import edu.msu.cme.rdp.probematch.core.DPMAligner.DPMAlignment;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result.BitVector64Match;
import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.PrintStream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class SliceToPrimer {

    private static final Options options = new Options();

    static {
        options.addOption("e", "edit-dist", true, "Maximum edit distance for matches (applies to both fprimer and rprimer, default=3)");
        options.addOption("n", "no-ambiguities", false, "Don't allow ambiguity matches");
        options.addOption("k", "keep-primers", false, "Include primers in sliced sequences");
        options.addOption(null, "min-length", true, "Output filter, min length (only effects sequences output)");
        options.addOption(null, "max-length", true, "Output filter, max length (only effects sequences output)");
        options.addOption(null, "fedit-dist", true, "Maximum edit distance for forward primer (cannot be applied with --edit-dist)");
        options.addOption(null, "redit-dist", true, "Maximum edit distance for reverse primer (cannot be applied with --edit-dist)");
        options.addOption("o", "seq-out", true, "Output sequence file (default=stdout)");
        options.addOption("s", "stats-out", true, "Output stats file (default=stderr)");
    }

    private static PatternBitMask64[] translateStringPrimers(String[] primers, boolean allowAmbiguities, boolean reverse) {
        PatternBitMask64[] ret = new PatternBitMask64[primers.length];
        for (int index = 0; index < primers.length; index++) {
            String primer = (reverse) ? IUBUtilities.reverseComplement(primers[index]) : primers[index];
            ret[index] = new PatternBitMask64(primer, allowAmbiguities);
        }

        return ret;
    }

    private static class PrimerMatch {

        int primerIndex;
        int start;
        int end;
        int score;

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final PrimerMatch other = (PrimerMatch) obj;
            if (this.primerIndex != other.primerIndex) {
                return false;
            }
            if (this.start != other.start) {
                return false;
            }
            if (this.end != other.end) {
                return false;
            }
            if (this.score != other.score) {
                return false;
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 97 * hash + this.primerIndex;
            hash = 97 * hash + this.start;
            hash = 97 * hash + this.end;
            hash = 97 * hash + this.score;
            return hash;
        }
    }

    public static void main(String[] args) throws Exception {
        //args = "--fedit-dist 4 --redit-dist=4 -k --max-length=400 --min-length=280 -o java_sliced_edit4.fasta TGCGAYCCSAARGCBGACTC ATSGCCATCATYTCRCCGGA /scratch/fishjord/tae_kwon_primer_match/all_genomes.fasta".split(" ");
        PatternBitMask64[] fprimers;
        String[] fprimerStrs, rprimerStrs;
        PatternBitMask64[] rprimers;

        FastaWriter seqOut;
        PrintStream statsOut;

        int fEdit = 3;
        int rEdit = 3;
        int minLength = Integer.MIN_VALUE;
        int maxLength = Integer.MAX_VALUE;
        boolean allowAmbiguities = true;
        boolean keepPrimers = false;

        SequenceReader inSeqs;

        try {

            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("edit-dist")) {
                fEdit = rEdit = Integer.parseInt(line.getOptionValue("edit-dist"));

                if (line.hasOption("redit-dist") || line.hasOption("fedit-dist")) {
                    throw new Exception("edit-dist, [fedit-dist, redit-dist] are mutually exclusive");
                }
            }

            if (line.hasOption("fedit-dist")) {
                fEdit = Integer.parseInt(line.getOptionValue("fedit-dist"));
            }

            if (line.hasOption("no-ambiguities")) {
                allowAmbiguities = false;
            }

            if (line.hasOption("keep-primers")) {
                keepPrimers = true;
            }

            if (line.hasOption("redit-dist")) {
                rEdit = Integer.parseInt(line.getOptionValue("redit-dist"));
            }

            if (line.hasOption("seq-out")) {
                seqOut = new FastaWriter(new File(line.getOptionValue("seq-out")));
            } else {
                throw new Exception("Must specify seq-out");
            }

            if (line.hasOption("stats-out")) {
                statsOut = new PrintStream(new File(line.getOptionValue("stats-out")));
            } else {
                statsOut = System.out;
            }

            if(line.hasOption("min-length")) {
                minLength = Integer.parseInt(line.getOptionValue("min-length"));
            }

            if(line.hasOption("max-length")) {
                maxLength = Integer.parseInt(line.getOptionValue("max-length"));
            }

            args = line.getArgs();

            if (args.length != 3) {
                throw new Exception("Unexpected number of command line arguments");
            }

            fprimers = translateStringPrimers(args[0].split(","), allowAmbiguities, false);
            fprimerStrs = args[0].split(",");
            rprimers = translateStringPrimers(args[1].split(","), allowAmbiguities, true);
            rprimerStrs = args[1].split(",");
            inSeqs = new SequenceReader(new File(args[2]));

        } catch (Exception e) {
            new HelpFormatter().printHelp("SliceToPrimer [options] <f,p,r,i,m,e,r> <r,p,r,i,m,e,r> <in_seq_file>", options);
            System.err.println("ERROR: " + e.getMessage());
            return;
        }

        Sequence seq;

        statsOut.println("orig_seqid\tsliced_seqid\tfprimer\tstart\tend\tscore\trprimer\tstart\tend\tscore\tlength");

        ScoringMatrix sccoringMatrix = ScoringMatrix.getDefaultNuclMatrix();

        DPMAligner[] faligners = new DPMAligner[fprimers.length];
        for(int index = 0;index < faligners.length;index++) {
            faligners[index] = new DPMAligner(fprimerStrs[index], Integer.MAX_VALUE);
        }

        try {
            while ((seq = inSeqs.readNextSequence()) != null) {
                Set<PrimerMatch> fprimerMatches = new HashSet();
                Set<PrimerMatch> rprimerMatches = new HashSet();

                for (int index = 0; index < fprimers.length; index++) {
                    PatternBitMask64 primer = fprimers[index];

                    for (BitVector64Match r : BitVector64.process(seq.getSeqString().toCharArray(), primer, fEdit).getResults()) {
                        PrimerMatch match = new PrimerMatch();
                        match.start = r.getPosition() - (primer.getPatternLength() + r.getScore());
                        match.end = r.getPosition();
                        match.score = r.getScore();
                        match.primerIndex = index;
                        fprimerMatches.add(match);
                    }
                }

                for (int index = 0; index < rprimers.length; index++) {
                    PatternBitMask64 primer = rprimers[index];

                    for (BitVector64Match r : BitVector64.process(seq.getSeqString().toCharArray(), primer, rEdit).getResults()) {
                        PrimerMatch match = new PrimerMatch();
                        match.start = r.getPosition() - (primer.getPatternLength() + r.getScore());
                        match.end = r.getPosition();
                        match.score = r.getScore();
                        match.primerIndex = index;
                        rprimerMatches.add(match);
                    }
                }

                if (fprimerMatches.isEmpty() || rprimerMatches.isEmpty()) {
                    statsOut.println(seq.getSeqName() + "\tEither/or no forward/reverse primer hits");
                    continue;
                }
                for (PrimerMatch fmatch : fprimerMatches) {
                    PrimerMatch bestReverse = null;
                    int bestScore = Integer.MAX_VALUE;
                    for (PrimerMatch rmatch : rprimerMatches) {
                        if (rmatch.start > fmatch.end && rmatch.start - fmatch.end < bestScore) {
                            bestReverse = rmatch;
                            bestScore = rmatch.start - fmatch.end;
                        }
                    }

                    if(bestReverse == null) {
                        statsOut.println(seq.getSeqName() + "\tNo reverse primer before " + fmatch.end);
                        continue;
                    }

                    String slicedSeq = null;
                    if (keepPrimers) {
                        slicedSeq = seq.getSeqString().substring(fmatch.start, bestReverse.end);
                    } else {
                        slicedSeq = seq.getSeqString().substring(fmatch.end, bestReverse.start);
                    }

                    String seqid = seq.getSeqName() + "_" + fmatch.primerIndex + "_" + fmatch.start;
                    if (slicedSeq.length() > minLength && slicedSeq.length() < maxLength) {
                        seqOut.writeSeq(seqid, "", slicedSeq);
                    }

                    DPMAlignment seqs = faligners[fmatch.primerIndex].align(seq.getSeqString().substring(fmatch.start, fmatch.end));
                    System.err.println(">" + seqid);
                    System.err.println(fprimerStrs[fmatch.primerIndex]);
                    System.err.println(seq.getSeqString().substring(fmatch.start, fmatch.end));
                    System.err.println();
                    System.err.println(seqs.getAlignedMatchFragment());
                    System.err.println(seqs.getAlignedProbe());

                    System.err.println();

                    statsOut.println(seq.getSeqName() + "\t" + seqid + "\t"
                            + fmatch.primerIndex + "\t" + fmatch.start + "\t" + fmatch.end + "\t" + fmatch.score + "\t"
                            + bestReverse.primerIndex + "\t" + bestReverse.start + "\t" + bestReverse.end + "\t" + bestReverse.score + "\t"
                            + slicedSeq.length());
                }
            }
        } catch(Exception e) {
            e.printStackTrace();
        } finally {
            statsOut.close();
            seqOut.close();
        }
    }
}
