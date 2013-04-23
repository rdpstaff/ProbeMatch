package edu.msu.cme.rdp.probematch.myers99;

import edu.msu.cme.rdp.probematch.myers99.BitVector64Result.BitVector64Match;
import java.util.Set;
import java.util.HashSet;

/*
 * MyersBitVector.java
 *
 * Created on February 21, 2004, 8:42 AM
 */
/**
 *
 * @author Ryan Farris
 */
public class BitVector64 {

    public static BitVector64Result process(String str, PatternBitMask64 pattern) {
        return process(str.toCharArray(), pattern, Integer.MAX_VALUE);
    }

    public static BitVector64Result process(char[] str, PatternBitMask64 pattern, int maxScore) {
        long P = -1, M = 0, Ebit = (long) 1 << (Myers99.WORD_SIZE_64 - 1);
        int score = Myers99.WORD_SIZE_64;
        int position = 0 - pattern.getRemainder();

        BitVector64Match regionMatch = null, bestMatch = null;

        Set<BitVector64Match> ret = new HashSet();

        for (int i = 0; i < str.length + pattern.getRemainder(); i++) {
            position++;

            char c = (i < str.length) ? str[i] : 0;

            long U, X, Y;
            U = pattern.getCharMask(c);
            X = (((U & P) + P) ^ P) | U;
            U |= M;

            Y = P;
            P = M | ~(X | Y);
            M = Y & X;

            if ((P & Ebit) != 0) {
                score += 1;
            } else if ((M & Ebit) != 0) {
                score -= 1;
            }

            Y = P << 1;
            P = (M << 1) | ~(U | Y);
            M = Y & U;

            if (score <= maxScore) {
                if (regionMatch == null || score <= regionMatch.getScore()) {
                    regionMatch = new BitVector64Match(position, score);
                }

                if (bestMatch == null || score <= bestMatch.getScore()) {
                    bestMatch = regionMatch;
                }
            } else if (regionMatch != null) {
                ret.add(regionMatch);
                regionMatch = null;
            }
        }

        if (regionMatch != null) {
            ret.add(regionMatch);
            regionMatch = null;
        }

        return new BitVector64Result(ret, bestMatch);
    }

    public static BitVector64Match processBest(String str, PatternBitMask64 pattern) {
        return processBest(str.toCharArray(), pattern, Integer.MAX_VALUE);
    }

    public static BitVector64Match processBest(char[] str, PatternBitMask64 pattern, int maxScore) {
        long P = -1, M = 0, Ebit = (long) 1 << (Myers99.WORD_SIZE_64 - 1);
        int score = Myers99.WORD_SIZE_64;
        int position = 0 - pattern.getRemainder();

        int bestScore = maxScore;
        int pos = -1;

        for (int i = 0; i < str.length + pattern.getRemainder(); i++) {
            position++;

            char c = (i < str.length) ? str[i] : 0;

            long U, X, Y;
            U = pattern.getCharMask(c);
            X = (((U & P) + P) ^ P) | U;
            U |= M;

            Y = P;
            P = M | ~(X | Y);
            M = Y & X;

            if ((P & Ebit) != 0) {
                score += 1;
            } else if ((M & Ebit) != 0) {
                score -= 1;
            }

            Y = P << 1;
            P = (M << 1) | ~(U | Y);
            M = Y & U;

            if (score <= bestScore) {
                bestScore = score;
                pos = position;
            }
        }

        if (pos == -1) {
            return null;
        } else {
            return new BitVector64Match(pos, bestScore);
        }
    }
}
