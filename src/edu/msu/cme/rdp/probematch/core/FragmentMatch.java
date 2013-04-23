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

/**
 *
 * @author fishjord
 */
public class FragmentMatch {
    private int start;
    private int end;
    private int score;
    private String alignedFragment;
    private String alignedMatchRegion;

    public FragmentMatch(String alignedProbe, String alignedMatchRegion, int start, int end, int score) {
        this.start = start;
        this.end = end;
        this.score = score;
        this.alignedFragment = alignedProbe;
        this.alignedMatchRegion = alignedMatchRegion;
    }

    public String getAlignedMatchRegion() {
        return alignedMatchRegion;
    }

    public String getAlignedFragment() {
        return alignedFragment;
    }

    public int getEnd() {
        return end;
    }

    public int getScore() {
        return score;
    }

    public int getStart() {
        return start;
    }
}
