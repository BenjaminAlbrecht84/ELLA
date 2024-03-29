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

package ella.model.io;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import ella.model.aligner.aligning.extending.FrameRun;
import ella.model.aligner.aligning.extending.FrameShiftHit;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.aligner.utils.MyStrand;
import ella.model.io.SubjectManager.MySubject;

public class DaaWriter {

    private static final HashMap<Character, Integer> nucToIndex;

    static {
        nucToIndex = new HashMap<Character, Integer>();
        nucToIndex.put('A', 0);
        nucToIndex.put('C', 1);
        nucToIndex.put('G', 2);
        nucToIndex.put('T', 3);
    }

    private AtomicLong queryRecords = new AtomicLong(0), aliBlockSize = new AtomicLong(0), refNamesBlockSize = new AtomicLong(0),
            refLengthsBlockSize = new AtomicLong(0);
    private SubjectManager subjectManager = new SubjectManager();
    private File out;

    public DaaWriter(File out, IndexText text) {
        this.out = out;
    }

    public void writeHeader() {

        // filling long-section
        ArrayList<Byte> byteBuffer = new ArrayList<Byte>();
        long magicNumber = Long.valueOf("4327487858190246763");
        write(byteBuffer, readLittleEndian(magicNumber));
        long version = 0;
        write(byteBuffer, readLittleEndian(version));
        long diamondBuild = 0;
        write(byteBuffer, readLittleEndian(diamondBuild));
        long dbSeqs = MyParameters.REFERENCE_SEQUENCES;
        write(byteBuffer, readLittleEndian(dbSeqs));
        long dbSeqsUsed = 0;
        write(byteBuffer, readLittleEndian(dbSeqsUsed));
        long dbLetters = BlastStatisticsHelper.n;
        write(byteBuffer, readLittleEndian(dbLetters));
        long flags = 0;
        write(byteBuffer, readLittleEndian(flags));
        long queryRecords = 0;
        write(byteBuffer, readLittleEndian(queryRecords));

        // filling integer-section
        int modeRank = 3;
        write(byteBuffer, readLittleEndian(modeRank));
        write(byteBuffer, readLittleEndian(MyParameters.SCORING_MATRIX.getGapOpen()));
        write(byteBuffer, readLittleEndian(MyParameters.SCORING_MATRIX.getGapExtend()));
        int reward = 0;
        write(byteBuffer, readLittleEndian(reward));
        int penalty = MyParameters.FRAMESHIFT_PENALTY == null ? 0 : MyParameters.FRAMESHIFT_PENALTY;
        write(byteBuffer, readLittleEndian(penalty));
        int reserved1 = 0;
        write(byteBuffer, readLittleEndian(reserved1));
        int reserved2 = 0;
        write(byteBuffer, readLittleEndian(reserved2));
        int reserved3 = 0;
        write(byteBuffer, readLittleEndian(reserved3));

        // filling double-section
        write(byteBuffer, readLittleEndian(BlastStatisticsHelper.K));
        write(byteBuffer, readLittleEndian(BlastStatisticsHelper.LAMBDA));
        double reserved4 = 0;
        write(byteBuffer, readLittleEndian(reserved4));
        double reserved5 = 0;
        write(byteBuffer, readLittleEndian(reserved5));

        // filling block-section
        String matrixType = MyParameters.SCORING_MATRIX.getType();
        for (int i = 128; i < 144; i++) {
            byte b = i - 128 < matrixType.length() ? (byte) matrixType.charAt(i - 128) : 0;
            write(byteBuffer, readLittleEndian(b));
        }
        for (int i = 144; i < 2192; i++)
            write(byteBuffer, readLittleEndian((byte) 0));
        for (int i = 2192; i < 2448; i++) {
            int rank = i - 2192 + 1;
            byte blockTypeRank = rank < 4 ? (byte) rank : 0;
            write(byteBuffer, readLittleEndian(blockTypeRank));
        }

        // writing-out buffer
        byte[] stream = new byte[byteBuffer.size()];
        for (int i = 0; i < byteBuffer.size(); i++)
            stream[i] = byteBuffer.get(i);
        writeInFile(stream, stream.length, false);

    }

    private final ConcurrentHashMap<Integer, QueryContainer> rankToRuns = new ConcurrentHashMap<Integer, QueryContainer>();
    private final ArrayList<byte[]> streams = new ArrayList<>();
    private int currentRank = 0;
    private long aliCounter = 0;

    public synchronized void writeContainers(ArrayList<QueryContainer> containerSet) {
        for (QueryContainer qC : containerSet) {
            aliCounter += qC.getFrameRuns().size();
            rankToRuns.put(qC.getRank(), qC);
            while (rankToRuns.containsKey(currentRank)) {
                QueryContainer qCWrite = rankToRuns.get(currentRank);
                streams.add(writeHits(qCWrite.getFrameRuns()));
                rankToRuns.remove(currentRank);
                qCWrite.freeMemory();
                currentRank++;
            }
        }
        writeInFile(streams, true);
        streams.clear();
    }

    private final ArrayList<Byte> byteBuffer = new ArrayList<>();
    private final HashMap<String, byte[]> readIDToPackedSeq = new HashMap<>();
    private final ArrayList<int[]> allocPairs = new ArrayList<>();

    private byte[] writeHits(ArrayList<FrameRun> runs) {

        byteBuffer.clear();
        readIDToPackedSeq.clear();
        allocPairs.clear();

        String lastReadName = null;
        int begin = -1, alloc = -1;

        for (FrameRun run : runs) {

            String readName = run.getReadName();
            FrameShiftHit hit = run.getFrameShiftHit();

            if (lastReadName == null || !readName.equals(lastReadName)) {

                if (alloc != -1) {
                    alloc = byteBuffer.size() - 4 - begin;
                    int[] allocPair = {begin, alloc};
                    allocPairs.add(allocPair);
                }

                begin = byteBuffer.size();

                alloc = 0;
                write(byteBuffer, readLittleEndian(alloc));

                int totalQueryLength = run.getTotalQueryLenth();
                write(byteBuffer, readLittleEndian(totalQueryLength));

                String queryName = readName;
                write(byteBuffer, readLittleEndian(queryName));
                write(byteBuffer, (byte) 0);

                byte nFlag = 0;
                write(byteBuffer, nFlag);

                if (!readIDToPackedSeq.containsKey(queryName)) {
                    byte[] packedSequence = run.getPackedSequence();
                    readIDToPackedSeq.put(queryName, packedSequence);
                }
                byte[] packedSequence = readIDToPackedSeq.get(readName);
                write(byteBuffer, packedSequence);

            }

            // int subjectID = run.getProteinLocationIndex();
            int subjectID = subjectManager.addSubject(run.getProteinLocationIndex());
            write(byteBuffer, readLittleEndian(subjectID));

            byte typeFlags = 1 << 1;
            typeFlags |= 1 << 3;
            typeFlags |= 1 << 5;
            typeFlags |= run.getStrand() == MyStrand.FWD ? 0 : 1 << 6;
            write(byteBuffer, typeFlags);

            int rawScore = hit.getRawScore();
            write(byteBuffer, readLittleEndian(rawScore));

            int queryStart = hit.getQueryStart();
            queryStart = run.getStrand() == MyStrand.FWD ? queryStart : run.getTotalQueryLenth() - queryStart - 1;
            write(byteBuffer, readLittleEndian(queryStart));

            int refStart = hit.getRefStart();
            write(byteBuffer, readLittleEndian(refStart));

            ArrayList<Byte> editOperations = hit.getEditOperations();
            for (byte op : editOperations)
                write(byteBuffer, op);
            write(byteBuffer, (byte) 0);

            lastReadName = readName;

        }

        if (alloc != -1) {

            alloc = byteBuffer.size() - 4 - begin;
            int[] allocPair = {begin, alloc};
            allocPairs.add(allocPair);

            byte[] stream = new byte[byteBuffer.size()];
            for (int i = 0; i < byteBuffer.size(); i++)
                stream[i] = byteBuffer.get(i);

            for (int[] p : allocPairs) {
                int counter = 0;
                for (byte b : readLittleEndian(p[1])) {
                    int pos = p[0] + (counter++);
                    stream[pos] = b;
                }
            }

//			writeInFile(stream, stream.length, true);

            aliBlockSize.getAndAdd(stream.length);
            queryRecords.getAndAdd(allocPairs.size());

            return stream;

        }

        return null;

    }

    private void write(ArrayList<Byte> byteBuffer, byte b) {
        byteBuffer.add(b);
    }

    private void write(ArrayList<Byte> byteBuffer, byte[] bytes) {
        for (byte b : bytes)
            byteBuffer.add(b);
    }

    private byte[] readLittleEndian(String s) {
        ByteBuffer buffer = ByteBuffer.allocate(s.length());
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        for (int i = 0; i < s.length(); i++)
            buffer.put((byte) s.charAt(i));
        return buffer.array();
    }

    private byte[] readLittleEndian(byte o) {
        ByteBuffer buffer = ByteBuffer.allocate(1);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.put(o);
        return buffer.array();
    }

    private byte[] readLittleEndian(int o) {
        ByteBuffer buffer = ByteBuffer.allocate(4);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(o);
        return buffer.array();
    }

    private byte[] readLittleEndian(double o) {
        ByteBuffer buffer = ByteBuffer.allocate(8);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putDouble(o);
        return buffer.array();
    }

    private byte[] readLittleEndian(long o) {
        ByteBuffer buffer = ByteBuffer.allocate(8);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        buffer.putLong(o);
        return buffer.array();
    }

    public void writeEnd(IndexText text) {
        try {

            // finishing alignment block
            writeInFile(readLittleEndian((int) 0), 4, true);
            aliBlockSize.getAndAdd(4);

            ArrayList<MySubject> subjects = subjectManager.getOrderedSubjectList();

            // inserting reference names
            ArrayList<Byte> byteBuffer = new ArrayList<Byte>();
            for (MySubject s : subjects) {
                String refName = text.getProteinAccessionFromIndex(s.getLocationIndex());
                write(byteBuffer, readLittleEndian(refName));
                write(byteBuffer, (byte) 0);
            }
            refNamesBlockSize.getAndSet(byteBuffer.size());

            // inserting reference length
            for (MySubject s : subjects) {
                int refLength = text.getProteinLengthFromIndex(s.getLocationIndex());
                write(byteBuffer, readLittleEndian(refLength));
            }
            refLengthsBlockSize.getAndSet(byteBuffer.size() - refNamesBlockSize.get());

            // writing-out buffer
            byte[] stream = new byte[byteBuffer.size()];
            for (int i = 0; i < byteBuffer.size(); i++)
                stream[i] = byteBuffer.get(i);
            writeInFile(stream, stream.length, true);

            // updating #dbSeqsUsed
            writeByteInFile(readLittleEndian((long) subjects.size()), 32);

            // updating #queryRecords
            writeByteInFile(readLittleEndian((long) queryRecords.longValue()), 56);

            // updating alignment block size
            writeByteInFile(readLittleEndian(aliBlockSize.get()), 144);

            // updating refNames block size
            writeByteInFile(readLittleEndian(refNamesBlockSize.get()), 152);

            // updating refLengths block size
            writeByteInFile(readLittleEndian(refLengthsBlockSize.get()), 160);

            System.out.println(aliCounter + " alignments reported");

        } catch (

                Exception e) {
            e.printStackTrace();
        }
    }

    private synchronized void writeByteInFile(byte[] b, long pos) {
        try {
            RandomAccessFile raf = new RandomAccessFile(out, "rw");
            try {
                raf.seek(pos);
                raf.write(b);
            } finally {
                raf.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private synchronized void writeInFile(ArrayList<byte[]> bytes, boolean append) {
        try {
            OutputStream output = null;
            try {
                output = new BufferedOutputStream(new FileOutputStream(out, append));
                for (byte[] b : bytes) {
                    if (b != null)
                        output.write(b, 0, b.length);
                }
            } finally {
                if (output != null)
                    output.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private synchronized void writeInFile(byte[] b, int len, boolean append) {
        try {
            OutputStream output = null;
            try {
                output = new FileOutputStream(out, append);
                output.write(b, 0, len);
            } finally {
                if (output != null)
                    output.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String getTotalQueryDNA(byte[] unpackedSequence) {
        char[] sigma = {'A', 'C', 'G', 'T'};
        StringBuilder buf = new StringBuilder();
        for (byte a : unpackedSequence) {
            buf.append(sigma[a]);
        }
        return buf.toString();
    }

    public String toStringUnpacked(byte[] unpacked) {
        StringBuilder buf = new StringBuilder();
        for (byte a : unpacked)
            buf.append(String.format("%d", a));
        return buf.toString();
    }

    public static byte[] getUnpackedSequence(byte[] packed, int query_len, int bits) {
        byte[] result = new byte[query_len];
        long x = 0;
        int n = 0, l = 0;
        int mask = (1 << bits) - 1;

        for (int i = 0; i < packed.length; i++) {
            x |= (packed[i] & 0xFF) << n;
            n += 8;

            while (n >= bits && l < query_len) {
                result[l] = (byte) (x & mask);
                n -= bits;
                x >>>= bits;
                l++;
            }
        }
        return result;
    }

}
