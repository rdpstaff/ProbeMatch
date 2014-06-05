/*
 * PatternBitMask.java
 *
 * Created on February 21, 2004, 8:44 AM
 */
package edu.msu.cme.rdp.probematch.myers99;

/**
 *
 * @author  Ryan Farris
 */
public class PatternBitMask64 {

    public static final long one = 1;
    public final long[] bitMask = new long[Myers99.SIGMA];
    public final int remainder;
    private final int patternLength;
    private String primerName= "NA";
    
    /** Creates a new instance of PatternBitMask 
     *  NEED TO ADD IUPAC   
     */
    public PatternBitMask64(String patternString, boolean matchAmbiguity) {
        this(patternString.toLowerCase().toCharArray(), matchAmbiguity);
    }

    public PatternBitMask64(String patternString, boolean matchAmbiguity, String primer) {
        this(patternString.toLowerCase().toCharArray(), matchAmbiguity);
        primerName = primer;
    }
    
    public PatternBitMask64(char[] pattern, boolean matchAmbiguity) {
        if (pattern.length > 64) {
            throw new IllegalStateException("PatternBitMask: pattern greater then 64.");
        }

        // set the masks
        long currentBit = 1;
        for (char c : pattern) {
            char lower = Character.toLowerCase(c);
            char upper = Character.toUpperCase(c);
            bitMask[lower] |= currentBit;
            bitMask[upper] |= currentBit;
            if (matchAmbiguity) {
                switch (lower) {
                    case 'r':
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        break;
                    case 'y':
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        break;
                    case 'm':
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        break;
                    case 'k':
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        break;
                    case 's':
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        break;
                    case 'w':
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        break;
                    case 'h':
                        bitMask['m'] |= currentBit;
                        bitMask['M'] |= currentBit;
                        bitMask['w'] |= currentBit;
                        bitMask['W'] |= currentBit;
                        bitMask['y'] |= currentBit;
                        bitMask['Y'] |= currentBit;
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        break;
                    case 'b':
                        bitMask['k'] |= currentBit;
                        bitMask['K'] |= currentBit;
                        bitMask['s'] |= currentBit;
                        bitMask['S'] |= currentBit;
                        bitMask['y'] |= currentBit;
                        bitMask['Y'] |= currentBit;
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        break;
                    case 'v':
                        bitMask['s'] |= currentBit;
                        bitMask['S'] |= currentBit;
                        bitMask['r'] |= currentBit;
                        bitMask['R'] |= currentBit;
                        bitMask['m'] |= currentBit;
                        bitMask['M'] |= currentBit;
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        break;
                    case 'd':
                    case 'i':
                        bitMask['r'] |= currentBit;
                        bitMask['R'] |= currentBit;
                        bitMask['k'] |= currentBit;
                        bitMask['K'] |= currentBit;
                        bitMask['w'] |= currentBit;
                        bitMask['W'] |= currentBit;
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        break;
                    case 'n':
                    case 'x':
                        bitMask['r'] |= currentBit;
                        bitMask['R'] |= currentBit;
                        bitMask['y'] |= currentBit;
                        bitMask['Y'] |= currentBit;
                        bitMask['m'] |= currentBit;
                        bitMask['M'] |= currentBit;
                        bitMask['k'] |= currentBit;
                        bitMask['K'] |= currentBit;
                        bitMask['s'] |= currentBit;
                        bitMask['S'] |= currentBit;
                        bitMask['w'] |= currentBit;
                        bitMask['W'] |= currentBit;
                        bitMask['h'] |= currentBit;
                        bitMask['H'] |= currentBit;
                        bitMask['b'] |= currentBit;
                        bitMask['B'] |= currentBit;
                        bitMask['v'] |= currentBit;
                        bitMask['V'] |= currentBit;
                        bitMask['d'] |= currentBit;
                        bitMask['D'] |= currentBit;
                        bitMask['i'] |= currentBit;
                        bitMask['I'] |= currentBit;
                        bitMask['g'] |= currentBit;
                        bitMask['G'] |= currentBit;
                        bitMask['a'] |= currentBit;
                        bitMask['A'] |= currentBit;
                        bitMask['c'] |= currentBit;
                        bitMask['C'] |= currentBit;
                        bitMask['t'] |= currentBit;
                        bitMask['T'] |= currentBit;
                        bitMask['u'] |= currentBit;
                        bitMask['U'] |= currentBit;
                        bitMask['x'] |= currentBit;
                        bitMask['X'] |= currentBit;
                        bitMask['n'] |= currentBit;
                        bitMask['N'] |= currentBit;
                        break;
                }
            }
            currentBit <<= 1;
        }

        // fill out the remainder of the mask with 1s
        remainder = Myers99.WORD_SIZE_64 - pattern.length;
        for (int i = pattern.length - 1; i < Myers99.WORD_SIZE_64; i++) {
            for (int j = 0; j < Myers99.SIGMA; j++) {
                bitMask[j] |= currentBit;
            }
            currentBit <<= 1;
        }

        patternLength = pattern.length;
    }
    
    public final long getCharMask(int character) {
        return bitMask[character];
    }

    public final int getRemainder() {
        return remainder;
    }

    public final int getPatternLength() {
        return patternLength;
    }
    
    public final String getPrimerName() {
        return primerName;
    }
}