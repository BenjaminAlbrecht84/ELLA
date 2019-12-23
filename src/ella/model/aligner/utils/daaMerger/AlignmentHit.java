/*
 * AlignmentHit.java Copyright (C) 2019. University of Tuebingen
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

import ella.model.aligner.utils.Alphabet;

import java.util.ArrayList;
import java.util.Objects;


public class AlignmentHit {

    public enum FrameDirection {
        POSITIVE, NEGATIVE
    }

    private final String readName;
    private String refName;
    private int readId;
    private int queryLength;
    private byte[] packedQuerySequence;
    private FrameDirection frame;
    private final int rawScore;
    private final int refStart;
    private final int queryStart;
    private final int refLength;
    private int refCoveredBases;
    private final int totalQueryLength;
    private int subjectID = -1;
    private final ArrayList<Byte> editOperations;
    private double identity;
    private final double bitScore;
    private double positiveRate;
    private int frameshifts;
    private final int accessPoint;
    private final long filePointer;

    public AlignmentHit(int ref_start, int bitScore, int rawScore, int queryStart, int refLength, int totalQueryLength, int subjectID, ArrayList<Byte> editOperations, int frame, String readName, String refName, long filePointer, int accessPoint) {

        cmpAllAliStatistics(editOperations);

        this.refStart = ref_start;
        this.refLength = refLength;
        this.totalQueryLength = totalQueryLength;
        this.queryStart = frame > 0 ? queryStart + 1 : queryStart - queryLength + 1;
        this.bitScore = bitScore;
        this.rawScore = rawScore;
        this.subjectID = subjectID;
        this.editOperations = editOperations;
        this.readName = readName;
        this.filePointer = filePointer;
        this.accessPoint = accessPoint;
        this.refName = refName;
    }

    public void cmpAllAliStatistics(ArrayList<Byte> editOperations) {
        double matches = 0, len = 0, positives = 0;
        queryLength = refCoveredBases = frameshifts = 0;
        for (Byte opByte : editOperations) {
            int op = opByte & 0xFF;
            switch (op >>> 6) {
                case (0): // handling match
                    queryLength += ((op & 63) * 3);
                    refCoveredBases += (op & 63);
                    matches += (op & 63);
                    positives += (op & 63);
                    len += (op & 63);
                    break;
                case (1): // handling insertion
                    queryLength += ((op & 63) * 3);
                    len += (op & 63);
                    break;
                case (2): // handling deletion
                    refCoveredBases += 1;
                    len++;
                    break;
                case (3): // handling substitution
                    char c = Alphabet.getAaString().charAt(op & 31);
                    if (c == '/') {
                        queryLength -= 1;
                        frameshifts++;
                    } else if (c == '\\') {
                        queryLength += 1;
                        frameshifts++;
                    } else {
                        int posBit = (opByte >> 5) & 1;
                        if (posBit == 1)
                            positives++;
                        queryLength += 3;
                        refCoveredBases += 1;
                    }
                    len++;
                    break;

            }
        }
        positiveRate = (positives / len) * 100;
        identity = (matches / len) * 100;
    }

    public AlignmentHit(int ref_start, int bitScore, int rawScore, int queryStart, int refLength, int queryLength, int totalQueryLength, int subjectID,
                        int refCover, ArrayList<Byte> editOperations, int frame, int frameshifts, String readName, long filePointer, int accessPoint) {
        this.refStart = ref_start;
        this.refLength = refLength;
        this.totalQueryLength = totalQueryLength;
        this.refCoveredBases = refCover;
        this.queryStart = frame > 0 ? queryStart + 1 : queryStart - queryLength + 1;
        this.queryLength = queryLength;
        this.bitScore = bitScore;
        this.rawScore = rawScore;
        this.subjectID = subjectID;
        this.editOperations = editOperations;
        this.frameshifts = frameshifts;
        this.readName = readName;
        this.filePointer = filePointer;
        this.accessPoint = accessPoint;

        cmpAllAliStatistics(editOperations);
    }

    public void cmpAliStatistics(ArrayList<Byte> editOperations) {
        double matches = 0, len = 0, positives = 0;
        for (Byte opByte : editOperations) {
            int op = opByte & 0xFF;
            switch (op >>> 6) {
                case (0): // handling match
                    positives += (op & 63);
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
                    int bit = (opByte >> 5) & 1;
                    if (bit == 1)
                        positives++;
                    len++;
                    break;
            }
        }
        positiveRate = (positives / len) * 100;
        identity = (matches / len) * 100;
    }

    public int getReadId() {
        return readId;
    }

    public void setReadId(int readId) {
        this.readId = readId;
    }

    public int getTotalRefLength() {
        return getRefLength();
    }

    public double getBitScore() {
        return bitScore;
    }

    public double getIdentity() {
        return identity;
    }

    public double getPositiveRate() {
        return positiveRate;
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

    public int getQueryEnd() {
        return queryStart + queryLength;
    }

    public FrameDirection getFrame() {
        return frame;
    }

    public int getFrameshifts() {
        return frameshifts;
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

    public byte[] getPackedQuerySequence() {
        return packedQuerySequence;
    }

    public int getQueryLength() {
        return queryLength;
    }

    public int getRefLength() {
        return refLength;
    }

    public int getRefCoveredBases() {
        return refCoveredBases;
    }

    public double getRefCoveredPercentage() {
        return ((double) refCoveredBases / refLength) * 100;
    }

    public String getRefName() {
        return refName;
    }

    public int getAccessPoint() {
        return accessPoint;
    }

    public long getFilePointer() {
        return filePointer;
    }

    public void setRefName(String refName) {
        this.refName = refName;
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof AlignmentHit) {
            AlignmentHit h = (AlignmentHit) o;
            if (h.getQueryStart() == queryStart && h.getSubjectID() == subjectID && h.getReadName().equals(readName) && h.getRawScore() == rawScore
                    && h.getRefStart() == refStart && h.getEditOperations() != null && editOperations != null
                    && h.getEditOperations().size() == editOperations.size()) {
                for (int i = 0; i < editOperations.size(); i++) {
                    if (!Objects.equals(editOperations.get(i), h.getEditOperations().get(i)))
                        return false;
                }
                return true;
            }
        }
        return false;
    }

    // FOR DEBUGGING ***********************************

    public void print(String prefix) {
        System.err.println(
                prefix + " " + "\tQB:[" + queryStart + ", " + getQueryEnd() + " ]\tRB:[" + refStart + ", ? ]\tRS: " + rawScore + "\tFR: " + frame);
    }

    public String toString() {
        return ("\tQB:[" + queryStart + ", " + getQueryEnd() + " | " + totalQueryLength + " ]\tRB:[" + refStart + ", " + (refStart + refCoveredBases)
                + " | " + refLength + "]\tRS: " + rawScore + "\tFR: " + frame + "\tID:" + identity + "\tPS:" + positiveRate + "\t" + refName + "\t" + readName);
    }

}
