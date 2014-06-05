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

package edu.msu.cme.rdp.probematch.core;

import java.util.Set;
import java.util.HashSet;

import edu.msu.cme.rdp.probematch.core.DPMAligner.DPMAlignment;
import edu.msu.cme.rdp.probematch.myers99.BitVector64;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result.BitVector64Match;
import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;

/**
 *
 * @author fishjord
 */
public class SeqFragMatcher {


    private String fragment;
    private PatternBitMask64 patternMask;
    private DPMAligner aligner;
    private int maxDist;

    public SeqFragMatcher(String fragment) {
        this(fragment, Integer.MAX_VALUE);
    }

    public SeqFragMatcher(String fragment, int maxDist) {
        this(fragment, maxDist, true);
    }

    public SeqFragMatcher(String fragment, int maxDist, boolean allowAmbiguities) {
        this.fragment = fragment;
        this.maxDist = maxDist;

        patternMask = new PatternBitMask64(fragment, allowAmbiguities);
        aligner = new DPMAligner(fragment, Integer.MAX_VALUE);
    }

    public FragmentMatch getBestMatch(String sequence) {
        return convertBitVector64Match(sequence, BitVector64.process(sequence, patternMask).getBestResult());
    }

    public Set<FragmentMatch> getMatches(String sequence) {
        Set<FragmentMatch> ret = new HashSet();

        for(BitVector64Match match : BitVector64.process(sequence, patternMask).getResults()) {
            ret.add(convertBitVector64Match(sequence, match));
        }

        return ret;
    }

    private FragmentMatch convertBitVector64Match(String sequence, BitVector64Match match) {
        int worstCaseStart = match.getPosition() - (fragment.length() + match.getScore());

        DPMAlignment alignment = aligner.align(sequence.substring(worstCaseStart, match.getPosition()));

        return new FragmentMatch(alignment.getAlignedMatchFragment(), alignment.getAlignedProbe(), match.getPosition() - alignment.getAlignedMatchFragment().replace("-", "").length(), match.getPosition(), match.getScore());
    }
}
