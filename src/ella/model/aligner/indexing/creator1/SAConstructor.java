package ella.model.aligner.indexing.creator1;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.indexing.IndexWriter;
import ella.model.aligner.indexing.SuffixArrayTable;
import ella.model.aligner.indexing.creator1.IndexCreator.InputParser;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.aligner.utils.suffixArray.EnhancedSuffixArray;
import ella.model.aligner.utils.suffixArray.MSDRadixSortForSuffixes;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class SAConstructor {

    private boolean isStopped = false;
    private IndexText indexText;
    private CountDownLatch latch;
    private int cores;
    private ExecutorService executor;
    private ArrayList<SuffixArrayConstructor> saConstructors;
    private File outFile;
    private IndexCreator indexCreator;

    private SuffixArrayTable suffixArrayTable;

    public SAConstructor(IndexText text, int sepDepth, int p, int q, int cores, File outFile, IndexCreator indexCreator) {
        this.cores = cores;
        this.saConstructors = new ArrayList<>();
        this.suffixArrayTable = new SuffixArrayTable(text, sepDepth, p, q);
        this.outFile = outFile;
        this.indexText = text;
        this.indexCreator = indexCreator;
    }

    public void start() {

        long time = System.currentTimeMillis();
        indexCreator.setMaxProgress(saConstructors.size());
        latch = new CountDownLatch(saConstructors.size());
        executor = Executors.newFixedThreadPool(cores);
        for (SuffixArrayConstructor cT : saConstructors)
            executor.execute(cT);
        try {
            latch.await();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        executor.shutdown();

        writeSuffixArrayTable();
        long runtime = (System.currentTimeMillis() - time) / 1000;
        indexCreator.reportFinish("("+runtime + "s)", 2);

        // writing pointers
        try {
            IndexWriter.writePointers2File(outFile);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        suffixArrayTable.clearTable();
        suffixArrayTable = null;
        for (SuffixArrayConstructor cT : saConstructors)
            cT.freeMemory();

    }

    public void addSAConstructor(ByteArrayWrapper prefix, CyclicSeedShape[] seedShapes, IndexText indexText, ArrayList<InputParser> inputParsers,
                                 long size) {
        saConstructors.add(new SuffixArrayConstructor(prefix, seedShapes, indexText, inputParsers, size));

    }

    private synchronized void reportSuffixArray(EnhancedSuffixArray esa) {
        suffixArrayTable.addSuffixArray(esa);
        if (Runtime.getRuntime().freeMemory() < (1. / 4.) * Runtime.getRuntime().totalMemory())
            writeSuffixArrayTable();
    }

    private void writeSuffixArrayTable() {
        try {
            IndexWriter.writeSuffixArrayTable2File(outFile, suffixArrayTable, indexText);
            suffixArrayTable.clearTable();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void stop() {
        isStopped = true;
    }

    public class SuffixArrayConstructor implements Runnable {

        private ByteArrayWrapper prefix;
        private CyclicSeedShape[] seedShapes;
        private IndexText text;
        private UnsignedIntArray suffixes;
        private ArrayList<InputParser> parsers;
        private long size;

        public SuffixArrayConstructor(ByteArrayWrapper prefix, CyclicSeedShape[] seedShapes, IndexText text, ArrayList<InputParser> parsers,
                                      long size) {
            this.prefix = prefix;
            this.seedShapes = seedShapes;
            this.text = text;
            this.parsers = parsers;
            this.size = size;
        }

        public void freeMemory() {
            text = null;
            seedShapes = null;
            suffixes = null;
            prefix = null;
        }

        public void run() {

            if (!isStopped) {

                // collecting text positions
                suffixes = new UnsignedIntArray(size);
                int i = 0;
                for (InputParser p : parsers) {
                    for (int j = 0; j < p.getPrefixManager().get(prefix).size(); j++)
                        suffixes.set(i++, p.getPrefixManager().get(prefix).getStream().get(j));
                    p.getPrefixManager().freeMemory(prefix);
                }

                // computing suffix-array
                MSDRadixSortForSuffixes radixSort = new MSDRadixSortForSuffixes();
                for (CyclicSeedShape seedShape : seedShapes) {
                    EnhancedSuffixArray esa = new EnhancedSuffixArray(prefix, seedShape);
                    esa.initSuffixArray(text, suffixes, radixSort);
                    reportSuffixArray(esa);
                }
                radixSort = null;
                suffixes = null;

            }
            latch.countDown();
            indexCreator.reportProgress(1, 2);

        }

        public UnsignedIntArray getSuffixes() {
            return suffixes;
        }

    }

}
