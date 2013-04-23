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
package edu.msu.cme.rdp.probematch.myers99;

import java.util.Set;
import java.util.Collections;

/**
 *
 * @author fishjord
 */
public class BitVector64Result {

    public static class BitVector64Match {

        private int position;
        private int score;

        public BitVector64Match(int position, int score) {
            this.position = position;
            this.score = score;
        }

        public int getPosition() {
            return position;
        }

        public int getScore() {
            return score;
        }
    }

    private Set<BitVector64Match> results;
    private BitVector64Match bestResult;

    public BitVector64Result(Set<BitVector64Match> results, BitVector64Match bestResult) {
        this.results = Collections.unmodifiableSet(results);
        this.bestResult = bestResult;
    }

    public BitVector64Match getBestResult() {
        return bestResult;
    }

    public Set<BitVector64Match> getResults() {
        return results;
    }
}
