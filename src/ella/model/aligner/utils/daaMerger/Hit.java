/*
 * Copyright 2017 Benjamin Albrecht
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package ella.model.aligner.utils.daaMerger;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import ella.model.aligner.utils.SparseString;

public class Hit {

    public enum FrameDirection {
        POSITIVE, NEGATIVE
    }

    private String readName;
    private int totalQueryLenth, queryLength;
    private byte[] packedQuerySequence;
    private FrameDirection frame;
    private int rawScore, refStart, queryStart, refLength, refCover;

    private int subjectID = -1;

    public ArrayList<Byte> editOperations;

    public Hit(int frame, int rawScore, int refStart, int queryStart, String subjectName, String readName, ArrayList<Byte> editOperations,
               ArrayList<Object[]> subjectInfo, int queryLength, int totalQueryLength, byte[] packedQuerySequence) {
        this.frame = frame < 0 ? FrameDirection.NEGATIVE : FrameDirection.POSITIVE;
        this.rawScore = rawScore;
        this.refStart = refStart;
        this.queryStart = queryStart;
        Object[] subject = {new SparseString(subjectName), null};
        this.subjectID = Collections.binarySearch(subjectInfo, subject, new InfoComparator());
        this.readName = readName;
        this.editOperations = editOperations;
        this.queryLength = queryLength;
        this.totalQueryLenth = totalQueryLength;
        this.packedQuerySequence = packedQuerySequence;
    }

    public Hit(int refStart, int bitScore, int rawScore, int queryStart, int refLength, int queryLength, int subjectID, int refCover,
               ArrayList<Byte> editOperations) {
        this.refStart = new Integer(refStart);
        this.refLength = new Integer(refLength);
        this.refCover = refCover;
        this.queryStart = new Integer(queryStart);
        this.rawScore = new Integer(rawScore);
        this.subjectID = subjectID;
        this.editOperations = editOperations;
    }

    public double getIdentity() {
        double matches = 0, len = 0;
        for (Byte opByte : editOperations) {
            int op = opByte & 0xFF;
            switch (op >>> 6) {
                case (0): // handling match
                    matches += (op & 63);
                    len += (op & 63);
                    break;
                case (1): // handling insertion
                    len += (op & 63);
                    break;
                case (2): // handling deletion
                    len++;
                    break;
                case (3): // handling substitution
                    len++;
                    break;
            }
        }
        return (matches / len) * 100;
    }

    public int getRefStart() {
        return refStart;
    }

    public int getRawScore() {
        return rawScore;
    }

    public int getQueryStart() {
        return queryStart;
    }

    public FrameDirection getFrame() {
        return frame;
    }

    public void setFrame(int frame) {
        this.frame = frame > 0 ? FrameDirection.POSITIVE : FrameDirection.NEGATIVE;
    }

    public int getSubjectID() {
        return subjectID;
    }

    public ArrayList<Byte> getEditOperations() {
        return editOperations;
    }

    public String getReadName() {
        return readName;
    }

    public int getTotalQueryLenth() {
        return totalQueryLenth;
    }

    public byte[] getPackedQuerySequence() {
        return packedQuerySequence;
    }

    public int getQueryLength() {
        return queryLength;
    }

    public int getRefLength() {
        return refLength;
    }

    public int getRefCover() {
        return refCover;
    }

    public double getRefCoverPercentage() {
        return (double) refCover / (double) refLength;
    }

    private class InfoComparator implements Comparator<Object[]> {
        @Override
        public int compare(Object[] o1, Object[] o2) {
            SparseString s1 = (SparseString) o1[0];
            SparseString s2 = (SparseString) o2[0];
            return s1.toString().compareTo(s2.toString());
        }
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof Hit) {
            Hit h = (Hit) o;
            if (h.getQueryStart() == queryStart && h.getSubjectID() == subjectID && h.getReadName().equals(readName) && h.getRawScore() == rawScore
                    && h.getRefStart() == refStart && h.getEditOperations() != null && editOperations != null
                    && h.getEditOperations().size() == editOperations.size()) {
                for (int i = 0; i < editOperations.size(); i++) {
                    if (editOperations.get(i) != h.getEditOperations().get(i))
                        return false;
                }
                return true;
            }
        }
        return false;
    }

    // FOR DEBUGGING ***********************************

    public void print(String prefix) {
        System.out.println(prefix + " " + "\tQB:[" + queryStart + ", ? ]\tRB:[" + refStart + ", ? ]\tRS: " + rawScore + "\tFR: " + frame);
    }

    public String toString() {
        return ("\tQB:[" + queryStart + ", ? ]\tRB:[" + refStart + ", ? ]\tRS: " + rawScore + "\tFR: " + frame + "\tId: " + getIdentity() + "\tRC: " + getRefCoverPercentage());
    }

}
