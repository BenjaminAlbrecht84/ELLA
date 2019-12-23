package ella.model.aligner.aligning;

import java.io.File;
import java.util.ArrayList;
import java.util.Stack;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ella.model.Taskmanager;
import ella.model.aligner.aligning.extending.AlignmentExtender;
import ella.model.aligner.aligning.seeding.BasicSeedFinder;
import ella.model.aligner.aligning.seeding.chaining.ChainBuilder;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import ella.model.aligner.aligning.IndexReader.SuffixArrayReader;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.indexing.SuffixArrayTable;
import ella.model.aligner.utils.SparseString;
import ella.model.aligner.utils.daaMerger.DaaMerger;
import ella.model.aligner.utils.frameshiftAligner.DP_Matrix_last;
import ella.model.io.DaaWriter;
import ella.model.io.FastaReader;
import ella.model.io.MyParameters;

public class Aligner {

    public int counter = 0, seed_counter = 0;
    private ArrayList<Long> runtimes = new ArrayList<Long>();

    private ExecutorService executor;
    private CountDownLatch latch;
    private int progressCounter = 0;
    private ObservableList<Object[]> progressList = FXCollections.observableArrayList();
    private long maxProgress, progress, lastProgress;
    private DoubleProperty totalProgress = new SimpleDoubleProperty(0.);
    private String message = "";

    private File indexFile, queryFile, outputFile;
    private int m, cores;
    private boolean isStopped = false;
    private Taskmanager taskmanager;
    private ArrayList<QueryContainer> queryContainerSet;
    private int containerIndex = 0;

    public Aligner(File indexFile, File queryFile, File outputFile, int m, int cores) {
        this.indexFile = indexFile;
        this.queryFile = queryFile;
        this.outputFile = outputFile;
        this.m = m;
        this.cores = cores;
    }

    public void run() {

        long time = System.currentTimeMillis();

        // starting to read index file
        ArrayList<File> daaFiles = new ArrayList<File>();
        IndexReader indexReader = new IndexReader(indexFile);
        int batchCounter = 1;
        while (indexReader.run() && !isStopped) {

            long batchTime = System.currentTimeMillis();

            MyParameters.init(indexReader.getNumOfLetters(), indexReader.getNumOfSequences());
            IndexText text = indexReader.getIndexText();

            System.out.println("# Batch " + batchCounter + "/" + MyParameters.TOTAL_BATCHES);
            setUpProgressList();

            // loading suffix arrays into memory
            setStatus(">Batch " + batchCounter + "/" + MyParameters.TOTAL_BATCHES + " - Loading index structure...");
            message = "Batch " + batchCounter + "/" + MyParameters.TOTAL_BATCHES + " - Loading index structure...";
            System.out.println("> Loading index structure...");
            ArrayList<Runnable> loadingThreads = new ArrayList<Runnable>();
            for (SuffixArrayReader saReader : indexReader.getSaReaders()) {
                loadingThreads.add(new Thread(() -> {
                    if (!isStopped) {
                        saReader.run();
                        reportProgress(1, 1);
                        latch.countDown();
                    }
                }));
            }
            maxProgress = loadingThreads.size();

            runInParallel(loadingThreads, cores);
            reportFinish(1);
            indexReader.freeMemory();
            progressCounter++;

            // processing query file
            queryContainerSet = new ArrayList<>();
            containerIndex = 0;
            for (Object[] o : FastaReader.read(queryFile)) {
                QueryContainer q = new QueryContainer((SparseString) o[0], (byte[]) o[1], (int) o[2], 2, queryContainerSet.size(), this);
                // q.initFrameMinimizers(text, indexReader.getSaTable());
                queryContainerSet.add(q);
            }

            // initializing DAA_Writer
            File daaFile = new File(outputFile.getAbsolutePath().split(".daa")[0] + "_" + daaFiles.size() + ".daa");
            daaFile.deleteOnExit();
            DaaWriter daaWriter = new DaaWriter(daaFile, text);
            daaWriter.writeHeader();

            // setting up alignment threads
            setStatus("Batch " + batchCounter + "/" + MyParameters.TOTAL_BATCHES + " - Aligning to index structure...");
            message = "Batch " + batchCounter + "/" + MyParameters.TOTAL_BATCHES + " - Aligning to index structure...";
            System.out.println("> Aligning to index structure...");
            ArrayList<Runnable> seedingThreads = new ArrayList<Runnable>();
            for (int i = 0; i < cores; i++)
                seedingThreads.add(new AlignmentThread(text, indexReader.getSaTable(), m, daaWriter));

            // computing alignments in parallel
            maxProgress = queryContainerSet.size();
            runInParallel(seedingThreads, cores);
            reportFinish(2);
            progressCounter++;

            // writing end of DAA_File
            daaWriter.writeEnd(text);
            daaFiles.add(daaFile);

//			DAAEvaluator.run(daaFile);

            long runtime = (System.currentTimeMillis() - batchTime) / 1000;
            System.out.println("Runtime: " + (runtime / 60) + "min " + (runtime % 60) + "s\n");
            System.out.println();

            batchCounter++;

        }

        if (!outputFile.getAbsolutePath().endsWith(".daa"))
            outputFile = new File(outputFile.getAbsolutePath() + ".daa");
        DaaMerger daaMerger = new DaaMerger();
        daaMerger.run(outputFile, daaFiles, queryFile, cores, false, daaFiles.get(0), progressList, this);
        String result = daaMerger.getQueryCounter()+" reads aligned against "+daaMerger.getReferenceCounter()+" proteins resulting in "+daaMerger.getHitCounter()+" alignments.";
        System.out.println(result);
        setStatus(result);
        
        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.out.println("Total Runtime: " + (runtime / 60) + "min " + (runtime % 60) + "s\n");
        if (taskmanager != null)
            taskmanager.reportFinish();

    }

    private void setUpProgressList() {
        if (progressList.size() < MyParameters.TOTAL_BATCHES) {
            for (int i = 1; i <= MyParameters.TOTAL_BATCHES; i++) {
                Object[] data = {"Aligning to Batch " + i, "Loading Index", new SimpleDoubleProperty(), "Aligning to Index",
                        new SimpleDoubleProperty()};
                progressList.add(data);
            }
            Object[] data = {"Writing Output", "Processing DAA Files", new SimpleDoubleProperty(), "Merging DAA Files", new SimpleDoubleProperty()};
            progressList.add(data);
        }
    }

    public void stop() {
        isStopped = true;
        while (latch.getCount() != 0)
            latch.countDown();
    }

    private synchronized QueryContainer nextQueryContainer() {
        QueryContainer q = null;
        if (containerIndex < queryContainerSet.size())
            q = queryContainerSet.get(containerIndex++);
        return q;
    }

    public class AlignmentThread implements Runnable {

        private final DP_Matrix dpMatrix = new DP_Matrix();
        private final BasicSeedFinder basicSeedFinder = new BasicSeedFinder();
        private final AlignmentExtender alignmentExtender = new AlignmentExtender();
        private final ChainBuilder chainBuilder = new ChainBuilder();

        private final IndexText text;
        private final SuffixArrayTable saTable;
        private final DaaWriter daaWriter;
        private final int m, weightThreshold;
        private final ArrayList<QueryContainer> containerSet = new ArrayList<>();

        public AlignmentThread(IndexText text, SuffixArrayTable saTable, int m, DaaWriter daaWriter) {
            this.text = text;
            this.saTable = saTable;
            this.m = m;
            this.daaWriter = daaWriter;
            this.weightThreshold = Math.max(1, (int) (Math.round(100000. / (double) cores)));
        }

        @Override
        public void run() {
            QueryContainer q;
            int totalQueryWeight = 0;
            while ((q = nextQueryContainer()) != null) {
                if (!isStopped) {
                    q.computeSeedings(text, saTable, m, basicSeedFinder);
                    q.chainSeedings(text, chainBuilder);
                    q.computeGappedAlignments(text, dpMatrix, alignmentExtender);
                    containerSet.add(q);
                    totalQueryWeight += q.getWeight();
                }
                reportProgress(1, 2);
                if (totalQueryWeight > weightThreshold) {
                    daaWriter.writeContainers(containerSet);
                    containerSet.clear();
                    totalQueryWeight = 0;
                }
            }
            daaWriter.writeContainers(containerSet);
            latch.countDown();
        }

    }

    public void runInParallel(ArrayList<Runnable> threads, int cores) {
        executor = Executors.newFixedThreadPool(cores);
        latch = new CountDownLatch(threads.size());
        for (Runnable r : threads)
            executor.execute(r);
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        executor.shutdown();
    }

    private void reportProgress(int delta, int step) {
        progress += delta;
        int p = ((int) ((((double) progress / (double) maxProgress)) * 100) / 10) * 10;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% ");
            ((SimpleDoubleProperty) progressList.get(progressCounter / 2)[2 + ((progressCounter % 2) * 2)])
                    .setValue((double) progress / (double) maxProgress);
            totalProgress.setValue(((double) p / 2.) + 50. * (step - 1));
            setStatus(message + " (" + p + "%)");
        }
    }

    private void reportFinish(int step) {
        progress = 0;
        lastProgress = 0;
        System.out.print(100 + "%\n");
        totalProgress.setValue(50. * step);
        ((SimpleDoubleProperty) progressList.get(progressCounter / 2)[2 + ((progressCounter % 2) * 2)]).setValue(1);
    }

    public void addRuntime(int i, long time) {
        for (int j = runtimes.size(); j <= i; j++)
            runtimes.add(0L);
        long newTime = runtimes.get(i) + time;
        runtimes.set(i, newTime);
    }

    public void setStatus(String status) {
        if (taskmanager != null && !isStopped)
            taskmanager.setStatus(status);
    }

    public DoubleProperty totalProgressProperty() {
        return totalProgress;
    }

    public ObservableList<Object[]> getProgressList() {
        return progressList;
    }

    public int getCores() {
        return cores;
    }

    public void setTaskmanager(Taskmanager taskmanager) {
        this.taskmanager = taskmanager;
    }
}
