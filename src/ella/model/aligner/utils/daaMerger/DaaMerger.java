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

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import ella.model.aligner.aligning.Aligner;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.ObservableList;
import ella.model.aligner.utils.SparseString;
import ella.model.io.FastaReader;
import ella.model.io.MyParameters;

public class DaaMerger {

    private int maxProgress, lastProgress = 0;
    private AtomicInteger progress = new AtomicInteger();
    private ObservableList<Object[]> progressList;
    private int progressCounter;
    private Aligner aligner;

    private CountDownLatch latch;
    private ExecutorService executor;

    private long queryCounter, referenceCounter, hitCounter;

    public void run(File daaFile, ArrayList<File> daaFiles, File queryFile, int cores, boolean verbose, File headerFile, ObservableList<Object[]> progressList, Aligner aligner) {

        this.aligner = aligner;
        if (progressList != null) {
            this.progressList = progressList;
            this.progressCounter = progressList.size();
        }
        long time = System.currentTimeMillis();
        aligner.setStatus("Merging batch files...");
        System.out.println("Merging batch files into " + daaFile.getAbsolutePath() + "...");

        this.executor = Executors.newFixedThreadPool(cores);
        Header headerInfo = new Header();
        headerInfo.loadFromDAA(headerFile);

        ArrayList<DaaReader> daaReaders = new ArrayList<DaaReader>();
        for (File f : daaFiles)
            daaReaders.add(new DaaReader(f, false));

        // pre-processing daa files
        System.out.println(">Processing DAA files ");
        maxProgress = getTotalSeqUsed(daaReaders);
        ConcurrentSkipListSet<SubjectEntry> subjectInfoSet = new ConcurrentSkipListSet<SubjectEntry>();
        ArrayList<Thread> processThreads = generateProcessThreads(daaReaders, subjectInfoSet);
        runInParallel(processThreads);
        ArrayList<Object[]> subjectInfos = new ArrayList<Object[]>();
        Iterator<SubjectEntry> it = subjectInfoSet.iterator();
        while (it.hasNext()) {
            SubjectEntry e = it.next();
            Object[] subject = {e.getName(), e.getLength()};
            subjectInfos.add(subject);
        }
        if (verbose)
            System.out.println(subjectInfos.size() + " references processed!");
        referenceCounter = subjectInfos.size();
        reportFinish();
        progressCounter++;

        // parsing read information
        ArrayList<Object[]> readInfos = FastaReader.read(queryFile);

        // writing header of daa file
        DaaMergeWriter daaWriter = new DaaMergeWriter(daaFile);
        daaWriter.writeHeader(headerInfo.getDbSeqs(), headerInfo.getDbLetters(), headerInfo.getGapOpen(), headerInfo.getGapExtend(),
                headerInfo.getK(), headerInfo.getLambda());

        // writing hits into daa file
        aligner.setStatus("Writing final DAA file...");
        System.out.println(">Writing into DAA file: " + daaFile.getAbsolutePath());
        maxProgress = getTotalQueryRecords(daaReaders);
        ArrayList<Thread> batchReaders = new ArrayList<Thread>();
        for (DaaReader reader : daaReaders)
            batchReaders.add(new BatchReader(reader, subjectInfos));
        ArrayList<Hit> hits = new ArrayList<Hit>();
        queryCounter = readInfos.size();
        for (int i = 0; i < readInfos.size(); i++) {

            Object[] readInfo = readInfos.get(i);

            // reading-out hits in parallel
            for (Thread reader : batchReaders)
                ((BatchReader) reader).setReadName(readInfo);
            runInParallel(batchReaders);

            // storing hits
            ArrayList<Hit> batchHits = new ArrayList<Hit>();
            for (Thread reader : batchReaders)
                batchHits.addAll(((BatchReader) reader).getHits());
            if (MyParameters.DO_FILTERING) {
                ArrayList<Hit> dominantHits;
                if (batchHits.size() < HitFilterParallel.TRESHOLD)
                    dominantHits = HitFilter.run(batchHits, headerInfo.getLambda(), headerInfo.getK());
                else {
                    dominantHits = HitFilterParallel.run(batchHits, headerInfo.getLambda(), headerInfo.getK(), cores);
                    System.out.println("STEP 3 - Continuing writing into daa-file: " + daaFile.getAbsolutePath());
                }
                hits.addAll(dominantHits);
            } else {
                hits.addAll(batchHits);
            }

            // writing hits into daa file
            if (hits.size() > 10000 || i == readInfos.size() - 1) {
                hitCounter += hits.size();
                daaWriter.writeHits(hits);
                hits.clear();
            }

        }

        // writing subject info into daa file
        daaWriter.writeEnd(subjectInfos);

        // deleting daa files
        for (File f : daaFiles)
            f.deleteOnExit();

        reportFinish();
        if (verbose)
            System.out.println(hitCounter + " alignments written into DAA-File!");

        executor.shutdown();

        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.out.println("Runtime: " + (runtime / 60) + "min " + (runtime % 60) + "s\n");

    }

    public long getReferenceCounter() {
        return referenceCounter;
    }

    public long getHitCounter() {
        return hitCounter;
    }

    public long getQueryCounter() {
        return queryCounter;
    }

    private ArrayList<Hit> filterForUniqueHits(ArrayList<Hit> batchHits) {
        ArrayList<Hit> uniqueHits = new ArrayList<Hit>();
        for (int i = 0; i < batchHits.size(); i++) {
            boolean isUnique = true;
            for (int j = i + 1; j < batchHits.size(); j++) {
                if (batchHits.get(i).equals(batchHits.get(j))) {
                    isUnique = false;
                    break;
                }
            }
            if (isUnique)
                uniqueHits.add(batchHits.get(i));
        }
        return uniqueHits;
    }

    private int getTotalSeqUsed(ArrayList<DaaReader> daaReader) {
        int sum = 0;
        for (DaaReader reader : daaReader)
            sum += reader.getDaaHeader().getDbSeqsUsed();
        return sum;
    }

    private int getTotalQueryRecords(ArrayList<DaaReader> daaReader) {
        int sum = 0;
        for (DaaReader reader : daaReader)
            sum += reader.getDaaHeader().getNumberOfQueryRecords();
        return sum;
    }

    private void reportProgress(int delta) {
        progress.getAndAdd(delta);
        if (progressList != null) {
            ((SimpleDoubleProperty) progressList.get(progressCounter / 2)[2 + ((progressCounter % 2) * 2)]).setValue((double) progress.get() / (double) maxProgress);
        }
        int p = ((int) ((((double) progress.get() / (double) maxProgress)) * 100) / 10) * 10;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% ");
        }
    }

    private void reportFinish() {
        progress.set(0);
        lastProgress = 0;
        System.out.print(100 + "%\n");
        if (progressList != null)
            ((SimpleDoubleProperty) progressList.get(progressCounter / 2)[2 + ((progressCounter % 2) * 2)]).setValue(1.0);
    }

    public class BatchReader extends Thread {

        private ArrayList<Object[]> subjectInfo;
        private DaaReader daaReader;
        private Object[] readInfo;
        private ArrayList<Hit> hits;
        private int lastIndex = 0;
        private Long filePointer = null;

        private int index = 0;

        public BatchReader(DaaReader daaReader, ArrayList<Object[]> subjectInfo) {
            this.daaReader = daaReader;
            this.subjectInfo = subjectInfo;
        }

        public void run() {

            hits = new ArrayList<Hit>();
            Object[] result = daaReader.parseDAAHitByIndex(index, lastIndex, filePointer);
            ArrayList<DaaHit> daaHits = (ArrayList<DaaHit>) result[0];
            lastIndex = (int) result[1];
            filePointer = (long) result[2];

            if (index < daaReader.getDaaHeader().getNumberOfQueryRecords() && daaHits.get(0).getQueryName().equals(readInfo[0].toString())) {
                for (DaaHit daaHit : daaHits) {
                    int rawScore = daaHit.getRawScore();
                    String subjectName = daaHit.getReferenceName();
                    int refStart = daaHit.getRefStart();
                    String queryName = daaHit.getQueryName();
                    int queryStart = daaHit.getQueryStart();
                    int queryLength = daaHit.getQueryLength();
                    int totalQueryLength = daaHit.getTotalQueryLength();
                    int frame = daaHit.getFrame();
                    ArrayList<Byte> editOperations = daaHit.getEditByteOperations();
                    byte[] dnaSequence = daaHit.getPackedQuerySequence();
                    Hit h = new Hit(frame, rawScore, refStart, queryStart, subjectName, queryName, editOperations, subjectInfo, queryLength,
                            totalQueryLength, dnaSequence);
                    hits.add(h);
                }
                index++;

                if (index % 100 == 0)
                    reportProgress(100);

            }

            latch.countDown();

        }

        public void setReadName(Object[] readInfo) {
            this.readInfo = readInfo;
        }

        public ArrayList<Hit> getHits() {
            return hits;
        }

    }

    public void runInParallel(ArrayList<Thread> threads) {
        latch = new CountDownLatch(threads.size());
        for (Thread t : threads)
            executor.execute(t);
        try {
            latch.await();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public class ProcessThread extends Thread {

        private DaaReader reader;
        private ConcurrentSkipListSet<SubjectEntry> subjectInfo_Set;

        public ProcessThread(DaaReader reader, ConcurrentSkipListSet<SubjectEntry> subjectInfo_Set) {
            this.reader = reader;
            this.subjectInfo_Set = subjectInfo_Set;
        }

        public void run() {

            DaaHeader daaHeader = reader.getDaaHeader();
            try {
                for (int i = 0; i < (int) daaHeader.getDbSeqsUsed(); i++) {
                    SparseString refName = new SparseString(new String(daaHeader.getReferenceName(i)));
                    int refLength = daaHeader.getRefLength(i);
                    subjectInfo_Set.add(new SubjectEntry(refName, refLength));

                    if (i % 100 == 0)
                        reportProgress(100);

                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            latch.countDown();

        }

    }

    public ArrayList<Thread> generateProcessThreads(ArrayList<DaaReader> daaReader, ConcurrentSkipListSet<SubjectEntry> subjectInfo_Set) {
        ArrayList<Thread> processThreads = new ArrayList<Thread>();
        for (DaaReader reader : daaReader)
            processThreads.add(new ProcessThread(reader, subjectInfo_Set));
        return processThreads;

    }

    public class SubjectEntry implements Comparable<SubjectEntry> {

        private SparseString name;
        private int length;

        public SubjectEntry(SparseString name, int length) {
            this.name = name;
            this.length = length;
        }

        @Override
        public int compareTo(SubjectEntry o) {
            return name.compareTo(o.getName());
        }

        public SparseString getName() {
            return name;
        }

        public int getLength() {
            return length;
        }

    }

}
