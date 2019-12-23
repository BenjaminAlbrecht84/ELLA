/*
 * QueryHits.java Copyright (C) 2019. University of Tuebingen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package ella.model.aligner.utils.daaMerger;

import ella.model.aligner.utils.SparseString;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;

public class QueryHits {

    private final ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, ArrayList<AlignmentHit>>> hitMap = new ConcurrentHashMap<>();
    private int queryLength;
    private String queryName;

    public void add(AlignmentHit h, SparseString acc, int frame) {

        if (!hitMap.containsKey(acc))
            hitMap.put(acc, new ConcurrentHashMap<>());
        ConcurrentHashMap<Integer, ArrayList<AlignmentHit>> frameMap = hitMap.get(acc);

        if (!frameMap.containsKey(frame))
            frameMap.put(frame, new ArrayList<>());
        ArrayList<AlignmentHit> alignmentHits = frameMap.get(frame);

        alignmentHits.add(h);

    }

    public ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, ArrayList<AlignmentHit>>> getHitMap() {
        return hitMap;
    }

    public boolean isEmpty() {
        return hitMap.isEmpty();
    }

    public ArrayList<AlignmentHit> getAllHits() {
        ArrayList<AlignmentHit> allAlignmentHits = new ArrayList<>();
        for (SparseString acc : hitMap.keySet()) {
            for (int frame : hitMap.get(acc).keySet()) {
                allAlignmentHits.addAll(hitMap.get(acc).get(frame));
            }
        }
        return allAlignmentHits;
    }

    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (SparseString acc : hitMap.keySet()) {
            buf.append(">Acc: ").append(acc).append("\n");
            for (int frame : hitMap.get(acc).keySet()) {
                buf.append("\t>>Frame: ").append(frame).append("\n");
                for (AlignmentHit h : hitMap.get(acc).get(frame))
                    buf.append(h).append("\n");
            }
        }
        return buf.toString();
    }

    public int getQueryLength() {
        return queryLength;
    }

    public void setQueryLength(int queryLength) {
        this.queryLength = queryLength;
    }

    public String getQueryName() {
        return queryName;
    }

    public void setQueryName(String queryName) {
        this.queryName = queryName;
    }
}
