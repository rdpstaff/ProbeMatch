/*
 * DPMAligner.java
 *
 * Created on April 17, 2004, 3:50 PM
 */
package edu.msu.cme.rdp.probematch.core;

import edu.msu.cme.rdp.readseq.utils.IUBUtilities;

/**
 *
 * @author  farrisry
 */
public class DPMAligner {

    private char[] probe;
    private int maxHammingDistance;
    private boolean[] requiresMatch = null;
    private boolean mismatchesOnly;
    private int[][] dp_matrix = null;
    private int[][] dp_path = null;
    private int currentMaxSize = 0;
    
    private static final int NONE = 0;
    private static final int MATCH = 1;
    private static final int INSERT = 2;
    private static final int DELETE = 4;
    private int insertScore;
    private int deleteScore;
    private static final int mismatchScore = 1;
    private final int badScore;

    public static class DPMAlignment {

        private String alignedProbe;
        private String alignedMatchFragment;

        public DPMAlignment(String alignedProbe, String alignedMatchFragment) {
            this.alignedProbe = alignedProbe;
            this.alignedMatchFragment = alignedMatchFragment;
        }

        public String getAlignedMatchFragment() {
            return alignedMatchFragment;
        }

        public String getAlignedProbe() {
            return alignedProbe;
        }
    }

    /** Creates a new instance of DPMAligner */
    public DPMAligner(String probe, int maxHammingDistance) {
        this(probe.toCharArray(), maxHammingDistance);
    }

    public DPMAligner(char[] probe, int maxHammingDistance) {
        this.probe = probe;
        this.maxHammingDistance = maxHammingDistance;
        this.mismatchesOnly = false;
        this.badScore = maxHammingDistance * 2;
        this.insertScore = 1;
        this.deleteScore = 1;
    }

    public DPMAligner(char[] probe, int maxHammingDistance, boolean mismatchesOnly) {
        this(probe, maxHammingDistance);
        this.mismatchesOnly = mismatchesOnly;
        if (mismatchesOnly) {
            this.insertScore = this.badScore;
            this.deleteScore = this.badScore;
        }
    }

    public DPMAligner(char[] probe, int maxHammingDistance, boolean mismatchesOnly, boolean[] exactMatchPositions) {
        this(probe, maxHammingDistance, mismatchesOnly);
        if (exactMatchPositions.length != probe.length) {
            throw new RuntimeException("BUG: exactMatchPositions not the same length as probe.");
        }
        this.requiresMatch = exactMatchPositions;
    }

    public DPMAlignment align(String textString) {
        char[] matchingFragment = textString.toCharArray();

        if (probe.length > currentMaxSize) {
            currentMaxSize = probe.length * 2;
            dp_matrix = new int[currentMaxSize][currentMaxSize];
            dp_path = new int[currentMaxSize][currentMaxSize];
        }
        if (matchingFragment.length > currentMaxSize) {
            currentMaxSize = matchingFragment.length * 2;
            dp_matrix = new int[currentMaxSize][currentMaxSize];
            dp_path = new int[currentMaxSize][currentMaxSize];
        }

        for (int i = 0; i < matchingFragment.length + 1; i++) {
            dp_matrix[i][0] = 0;
            dp_path[i][0] = NONE;
        }
        for (int j = 0; j < probe.length + 1; j++) {
            dp_matrix[0][j] = j;
            dp_path[0][j] = NONE;
        }

        for (int i = 1; i < matchingFragment.length + 1; i++) {
            for (int j = 1; j < probe.length + 1; j++) {
                int match = dp_matrix[i - 1][j - 1];
                int insert = dp_matrix[i - 1][j] + insertScore;
                int delete = dp_matrix[i][j - 1] + deleteScore;

                if (!IUBUtilities.matches(matchingFragment[i - 1], probe[j - 1])) {
                    match += mismatchScore;
                    if (requiresMatch != null && requiresMatch[j - 1]) {
                        match += badScore;
                        insert += badScore;
                        delete += badScore;
                    }
                }

                if (match < insert && match < delete) {
                    dp_matrix[i][j] = match;
                    dp_path[i][j] = MATCH;
                } else if (insert < match && insert < delete) {
                    dp_matrix[i][j] = insert;
                    dp_path[i][j] = INSERT;
                } else if (delete < match && delete < insert) {
                    dp_matrix[i][j] = delete;
                    dp_path[i][j] = DELETE;
                } else if (delete == insert && delete == match) {
                    dp_matrix[i][j] = delete;
                    dp_path[i][j] = DELETE | INSERT | MATCH;
                } else if (delete == insert) {
                    dp_matrix[i][j] = delete;
                    dp_path[i][j] = DELETE | INSERT;
                } else if (delete == match) {
                    dp_matrix[i][j] = delete;
                    dp_path[i][j] = DELETE | MATCH;
                } else if (insert == match) {
                    dp_matrix[i][j] = match;
                    dp_path[i][j] = MATCH | INSERT;
                } else {
                    throw new RuntimeException("Impossible scoring state reached. This cannot happen. You're halucinating.");
                }
            }
        }

        int bestPosition = 0;
        for (int i = 0; i < matchingFragment.length + 1; i++) {
            if (dp_matrix[i][probe.length] <= dp_matrix[bestPosition][probe.length]) {
                bestPosition = i;
            }
        }

        int score = dp_matrix[bestPosition][probe.length];

        /******* ADD PATHING STUFFS HERE ************/
        if (score <= maxHammingDistance) {

            int i = bestPosition;
            int j = probe.length;
            int backwardsIndex = probe.length + matchingFragment.length;
            char[] alignedProbe = new char[backwardsIndex]; // as large as the max possible size
            char[] alignedMatchFragment = new char[backwardsIndex];
            boolean keepGoing = true;
            while (keepGoing) {
                switch (dp_path[i][j]) {
                    case MATCH:
                    case MATCH | INSERT:
                    case MATCH | INSERT | DELETE:
                    case MATCH | DELETE:
                        backwardsIndex--;
                        alignedProbe[backwardsIndex] = probe[j - 1];
                        alignedMatchFragment[backwardsIndex] = matchingFragment[i - 1];
                        
                        i--;
                        j--;
                        break;
                    case INSERT | DELETE:
                    case INSERT:
                        backwardsIndex--;
                        alignedProbe[backwardsIndex] = '-';
                        alignedMatchFragment[backwardsIndex] = matchingFragment[i - 1];
                        
                        i--;
                        break;
                    case DELETE:
                        backwardsIndex--;
                        alignedProbe[backwardsIndex] = probe[j - 1];
                        
                        alignedMatchFragment[backwardsIndex] = '-';
                        j--;
                        break;
                    case NONE:
                        keepGoing = false;
                        break;
                    default:
                        throw new RuntimeException("BUG: Got lost in the Matrix.");
                }
//                System.out.println();
            }

            // copy our buffers to the correct-sized arrays
            int alignedLength = probe.length + matchingFragment.length - backwardsIndex;
            char[] forAlignedProbe = new char[alignedLength];
            char[] forAlignedMatchFragment = new char[alignedLength];
            System.arraycopy(alignedProbe, backwardsIndex, forAlignedProbe, 0, alignedLength);
            System.arraycopy(alignedMatchFragment, backwardsIndex, forAlignedMatchFragment, 0, alignedLength);

            return new DPMAlignment(new String(forAlignedProbe), new String(forAlignedMatchFragment));
        } else {
            return null;
        }

    }
}
