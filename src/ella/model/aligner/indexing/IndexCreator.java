package ella.model.aligner.indexing;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ella.model.Taskmanager;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.Minimizer;
import ella.model.aligner.utils.PrefixArrayManager;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.streams.MyByteStream;
import ella.model.aligner.utils.streams.MyIntegerStream;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;
import ella.model.io.ByteWriter;
import ella.model.io.FastaCounter;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;

public class IndexCreator {

    private SAConstructor constructor;
    private String message = "";
    private Taskmanager taskmanager;
    private CountDownLatch latch;
    private long maxProgress, progress, lastProgress;
    private int batch, totalBatches;
    private DoubleProperty totalProgress = new SimpleDoubleProperty(0);
    private boolean isStopped = false;

    private final static long ID = 23062016;
    private final static int BUILD_NR = 110;
    private final static int VERSION_NR = 1;

    private IndexText indexText;

    public File run(File inFile, File edbFile, int sepDepth, CyclicSeedShape[] seedShapes, int p, int q, int cores, Long size) {

        sepDepth = sepDepth > 0 ? sepDepth : 1;

        // checking for inconsistency between separation-depth and seed-shapes
        for (CyclicSeedShape seedShape : seedShapes) {
            for (int i = 0; i < sepDepth; i++) {
                if (!seedShape.usePosition(i))
                    throw new IllegalArgumentException(
                            "SeedShape is incompatible with chosen separation depths: " + seedShape.toString() + " vs " + sepDepth);
            }
        }

        // creating index structures
        edbFile = edbFile == null ? createFile(inFile, "edb") : edbFile;
        createIndexData(inFile, sepDepth, seedShapes, p, q, cores, edbFile, size);

        return edbFile;

    }

    private File createFile(File inFile, String suffix) {
        String name = inFile.getName().split("\\.")[0] + "." + suffix;
        String path = inFile.getParent();
        String absolutePath = path + File.separatorChar + name;
        return new File(absolutePath);
    }

    public synchronized void reportProgress(int delta, int step) {
        progress += delta;
        int p = ((int) ((((double) progress / (double) maxProgress)) * 100) / 10) * 10;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% ");
            totalProgress.set((int) Math.round(((p * 0.33) + (step * 33.)) * (1. / totalBatches) + ((batch - 1) * totalBatches)));
            setStatus(message + " (" + p + "%)");
        }
    }

    public void reportFinish(String info, int step) {
        progress = 1;
        lastProgress = 0;
        System.out.print(100 + "% " + info + "\n");
        setStatus(message + " (" + 100 + "%)");
        int p = step == 0 ? 33 : step == 1 ? 66 : 100;
        totalProgress.set(p * (1. / totalBatches) + ((batch - 1) * totalBatches));
    }

    private void createIndexData(File inFile, int sepDepth, CyclicSeedShape[] seedShapes, int p, int q, int cores, File edbFile, Long size) {

        System.out.println("Parsing input file " + inFile.getAbsolutePath());

        long startTime = System.currentTimeMillis();

        System.out.println("> Scanning input file...");
        setStatus(">Scanning input file...");
        Object[] fileInfo = FastaCounter.run(inFile, sepDepth, size, p, q);
        int avgSequenceLength = (int) ((long) fileInfo[1] / (long) fileInfo[0]);
        long numOfSequences = (long) fileInfo[0];
        long[] batchBounderies = (long[]) fileInfo[2];
        long runtime = (System.currentTimeMillis() - startTime) / 1000;
        int batchCounter = 0;
        totalBatches = batchBounderies.length - 1;
        int approxNumOfBatchSequences = (int) (numOfSequences / totalBatches);
        System.out.println("#Letters: " + fileInfo[1] + " | #Sequences: " + fileInfo[0] + " | #Batches: " + totalBatches + " | (" + runtime + "s)");
        System.out.println("Uptime: " + (runtime / 60) + "min " + (runtime % 60) + "s");

        for (batch = 1; batch < batchBounderies.length; batch++) {

            batchCounter++;

            long numOfBatchLetters = 0, numOfBatchSequences = 0;
            long batchSize = batchBounderies[batch] - batchBounderies[batch - 1];
            long batchChunk = batchSize / cores;

            try {

                // parsing input file in parallel
                System.out.println("> Reading-in data of batch " + batchCounter + "/" + totalBatches + "...");
                message = "Batch " + batchCounter + "/" + totalBatches + " - Reading in data...";
                setStatus(">" + message);
                long time = System.currentTimeMillis();
                maxProgress = batchSize;
                ArrayList<InputParser> inputParsers = new ArrayList<InputParser>();
                for (long i = 0; i < cores; i++) {
                    long batchStart = batchBounderies[batch - 1] + i * batchChunk;
                    int approxNumOfChunkSequences = approxNumOfBatchSequences / cores;
                    InputParser parser = new InputParser(inFile, sepDepth, p, q, batchStart, batchChunk, approxNumOfChunkSequences,
                            avgSequenceLength);
                    inputParsers.add(parser);
                }
                runInParallel(inputParsers, cores);
                runtime = (System.currentTimeMillis() - time) / 1000;
                reportFinish("(" + runtime + "s)", 0);

                // adding text offset
                HashMap<ByteArrayWrapper, Long> prefixSizes = new HashMap<ByteArrayWrapper, Long>();
                for (InputParser parser : inputParsers) {
                    numOfBatchLetters += parser.getLetterCounter();
                    numOfBatchSequences += parser.getSequenceCounter();
                    for (ByteArrayWrapper prefix : parser.getPrefixManager().getPrefixes()) {
                        prefixSizes.putIfAbsent(prefix, 0L);
                        prefixSizes.put(prefix, prefixSizes.get(prefix) + parser.getPrefixManager().get(prefix).size());
                    }
                }

                // collecting results of each parser
                System.out.println("> Collecting data of batch " + batchCounter + "/" + totalBatches + "...");
                message = "Batch " + batchCounter + "/" + totalBatches + " - Collecting data...";
                setStatus(">" + message);
                time = System.currentTimeMillis();
                maxProgress = inputParsers.size();
                MyIntegerStream locationBuffer = new MyIntegerStream();
                MyByteStream textBuffer = new MyByteStream();
                for (int i = 0; i < inputParsers.size(); i++) {

                    if (!isStopped) {

                        InputParser parser = inputParsers.get(i);

                        // adding offset to locations
                        for (long j = 0; j < parser.getLocationBuffer().size(); j += 2)
                            parser.getLocationBuffer().addOffset(j, textBuffer.size());
                        locationBuffer.add(parser.getLocationBuffer());
                        parser.freeLocationBuffer();

                        // adding offset to text positions
                        for (ByteArrayWrapper prefix : parser.getPrefixManager().getPrefixes()) {
                            for (long j = 0; j < parser.getPrefixManager().get(prefix).size(); j += 1)
                                parser.getPrefixManager().get(prefix).addOffset(j, textBuffer.size());
                        }

                        // concatenating text
                        textBuffer.add(parser.getTextBuffer());
                        parser.freeTextBuffer();

                    }

                    // reporting progress
                    reportProgress(1, 1);

                }
                runtime = (System.currentTimeMillis() - time) / 1000;
                reportFinish("(" + runtime + "s)", 1);

                System.out.println("> Generating index structures for batch " + batchCounter + "/" + totalBatches + "...");
                message = "Batch " + batchCounter + "/" + totalBatches + " - Generating index structures...";
                setStatus(">" + message);
                time = System.currentTimeMillis();

                // generating index text
                indexText = new IndexText(textBuffer.toArray(), locationBuffer.toArray(), numOfBatchSequences, numOfBatchLetters, avgSequenceLength);
                locationBuffer = null;
                textBuffer = null;
                runtime = System.currentTimeMillis() - time;

                // writing text to file
                IndexWriter.writeIndexFile2File(indexText, edbFile, ID, BUILD_NR, VERSION_NR, batch == 1, totalBatches);

                // constructing suffix arrays
                IndexWriter.writeSuffixTableHeader2File(p, q, sepDepth, edbFile);
                constructor = new SAConstructor(indexText, sepDepth, p, q, cores, edbFile, this);
                Set<ByteArrayWrapper> prefixes = new PrefixArrayManager(sepDepth).getPrefixes();
                for (ByteArrayWrapper prefix : prefixes) {
                    if (prefixSizes.get(prefix) > 0)
                        constructor.addSAConstructor(prefix, seedShapes, indexText, inputParsers, prefixSizes.get(prefix));
                }
                constructor.start();

                // freeing memory
                for (InputParser parser : inputParsers)
                    parser.freeMemory();
                constructor = null;
                indexText = null;

                runtime = (System.currentTimeMillis() - startTime) / 1000;
                System.out.println("Uptime: " + (runtime / 60) + "min " + (runtime % 60) + "s");

            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

        }

        runtime = (System.currentTimeMillis() - startTime) / 1000;
        System.out.println("Runtime: " + (runtime / 60) + "min " + (runtime % 60) + "s\n");

    }

    public void runInParallel(ArrayList<InputParser> threads, int cores) {
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        latch = new CountDownLatch(threads.size());
        for (Runnable r : threads)
            executor.execute(r);
        try {
            latch.await();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        executor.shutdown();
    }

    public void stop() {
        isStopped = true;
        if (constructor != null)
            constructor.stop();
    }

    public class InputParser implements Runnable {

        private File dbFile;
        private int sepDepth, p, q;
        private long start, chunk;
        private int locationBufferSize, textBufferSize;

        private long letterCounter = 0;
        private int sequenceCounter = 0;
        private PrefixArrayManager prefixManager;
        private MyIntegerStream locationBuffer;
        private MyByteStream textBuffer;

        public InputParser(File dbFile, int sepDepth, int p, int q, long start, long chunk, int numOfSequences, int avgSequenceLength) {
            this.dbFile = dbFile;
            this.sepDepth = sepDepth;
            this.p = p;
            this.q = q;
            this.start = start;
            this.chunk = chunk;
            this.locationBufferSize = numOfSequences * 2 + 2;
            locationBufferSize = (int) Math.round(new Double(locationBufferSize) * 1.2);
            this.textBufferSize = numOfSequences * avgSequenceLength;
            textBufferSize = (int) Math.round(new Double(textBufferSize) * 1.2);
        }

        public void freeMemory() {
            prefixManager.clear();
            prefixManager = null;
            locationBuffer = null;
            textBuffer = null;
        }

        public void freeTextBuffer() {
            textBuffer = null;
        }

        public void freeLocationBuffer() {
            locationBuffer = null;
        }

        @Override
        public void run() {

            // calculating all possible prefixes
            prefixManager = new PrefixArrayManager(sepDepth);

            // parsing index data
            locationBuffer = new MyIntegerStream(locationBufferSize);
            textBuffer = new MyByteStream(textBufferSize);

            try {
                RandomAccessFile raf = new RandomAccessFile(dbFile, "r");
                raf.seek(start);
                byte[] readBuffer = new byte[10000];
                StringBuffer sequenceBuffer = new StringBuffer(), idBuffer = new StringBuffer();
                int charCounter = 0, readChars = 0;
                boolean parseSeq = false, parseID = false;
                while ((readChars = raf.read(readBuffer)) != -1) {

                    if (isStopped)
                        break;

                    for (int i = 0; i < readChars; i++) {
                        char c = (char) readBuffer[i];
                        charCounter++;
                        if (!parseSeq && c == '\n')
                            parseSeq = true;
                        else if (c == '>') {
                            if (idBuffer.length() > 0 && sequenceBuffer.length() > 0) {
                                storeProteinData(idBuffer, sequenceBuffer);
                                sequenceCounter++;
                            }
                            if (charCounter > chunk) {
                                break;
                            }
                            parseSeq = false;
                            parseID = true;
                            sequenceBuffer.setLength(0);
                            idBuffer.setLength(0);
                        } else if (parseSeq && c != '\n') {
                            char aa = c;
                            sequenceBuffer.append(aa);
                        } else if (parseID && (c == '\n' || c == ' '))
                            parseID = false;
                        else if (parseID) {
                            idBuffer.append(c);
                        }

                    }
                    reportProgress(readChars, 0);
                    if (readChars < readBuffer.length) {
                        storeProteinData(idBuffer, sequenceBuffer);
                        sequenceCounter++;
                    }
                    if (charCounter > chunk)
                        break;
                }
                readBuffer = null;

                // finishing reading file
                raf.close();

            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

            latch.countDown();

        }

        private void storeProteinData(StringBuffer id, StringBuffer sequenceBuffer) {

            // storing position and length of each minimizer within the sequence
            Minimizer minimizer = new Minimizer(p, q);
            int[] lastMin = null;
            int start = 0;
            for (int l = sequenceBuffer.length(); l >= Math.max(sepDepth, p); l--) {
                int[] coordinates = {start, l};
                minimizer.addPMer(coordinates, getPMer(sequenceBuffer, start, p));
                start++;
                if (start + p - 2 >= q - 1) {
                    int[] min = (int[]) minimizer.getMinimium();
                    if (lastMin == null || min != lastMin) {
                        lastMin = (int[]) min;
                        ByteArrayWrapper kmer = getKMer(sequenceBuffer, min[0], sepDepth);
                        long pos = textBuffer.size() + min[0];
                        prefixManager.addValue(kmer, pos);
                    }
                }
            }

            // writing sequence location and its length
            locationBuffer.add(textBuffer.size());
            locationBuffer.add(sequenceBuffer.length());

            // writing protein block existing of the sequence, the id, and the sequence length
            ByteWriter.write(textBuffer, packSequence(sequenceBuffer.toString()));
            ByteWriter.write(textBuffer, (byte) 127);
            ByteWriter.write(textBuffer, id.toString());
            ByteWriter.write(textBuffer, (byte) 127);

            letterCounter += sequenceBuffer.length();

        }

        public long getTextSize() {
            return textBuffer.size();
        }

        public long getLocationSize() {
            return locationBuffer.size();
        }

        public int getSequenceCounter() {
            return sequenceCounter;
        }

        public long getLetterCounter() {
            return letterCounter;
        }

        public PrefixArrayManager getPrefixManager() {
            return prefixManager;
        }

        public MyIntegerStream getLocationBuffer() {
            return locationBuffer;
        }

        public MyByteStream getTextBuffer() {
            return textBuffer;
        }

    }

    private String getPMer(StringBuffer sequence, int pos, int p) {
        StringBuilder builder = new StringBuilder(p);
        for (int i = 0; i < p; i++) {
            char c = sequence.charAt(pos + i);
            builder.append(c);
        }
        return builder.toString();
    }

    private ByteArrayWrapper getKMer(StringBuffer sequence, int pos, int k) {
        byte[] kmer = new byte[k];
        for (int i = 0; i < k; i++) {
            byte b = 0;
            b |= Alphabet.getIndex(Alphabet.reduceCharacter(sequence.charAt(pos + i)));
            kmer[i] = b;
        }
        return new ByteArrayWrapper(kmer);
    }

    private byte[] packSequence(String protein) {
        byte[] a = new byte[protein.length()];
        for (int i = 0; i < protein.length(); i++) {
            char c = protein.charAt(i);
            byte b = 0;
            b |= Alphabet.getIndex(c);
            a[i] = b;
        }
        return a;
    }

    public IndexText getIndexText() {
        return indexText;
    }

    public double getTotalProgress() {
        return totalProgress.get();
    }

    public DoubleProperty totalProgressProperty() {
        return totalProgress;
    }

    public void setMaxProgress(long maxProgress) {
        this.maxProgress = maxProgress;
    }

    public void setStatus(String status) {
        if (taskmanager != null && !isStopped)
            taskmanager.setStatus(status);
    }

    public void setTaskmanager(Taskmanager taskmanager) {
        this.taskmanager = taskmanager;
    }

}
