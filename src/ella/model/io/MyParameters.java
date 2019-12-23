package ella.model.io;

import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.aligner.utils.ScoringMatrix;
import ella.model.aligner.utils.CyclicSeedShape;

import java.io.File;

public class MyParameters {

    //input parameters
    public static int CPUS = Runtime.getRuntime().availableProcessors();

    // technical parameters
    public static long MAX_MEMORY = Runtime.getRuntime().maxMemory();
    public final static int INITIAL_STREAM_CAPACITIY = 16;

    // database properties
    public static int TOTAL_BATCHES = 1;

    // search parameters
    public static int MULTIPLICITY = 10;
    public static int E = 10;
    public static final int MAX_LCP = 127;

    // output parameters
    public static int MIN_COVERAGE = 0;
    public static double MAX_EVALUE = 0.001;
    public static int MIN_BITSCORE = 50;
    public static double MIN_PROPORTION_COVERAGE = 0.9;
    public static double MIN_PROPORTION_SCORE = 0.9;
    public static boolean DO_FILTERING = false;

    // different seed shapes adopted from DIAMOND and LAST
    private static final CyclicSeedShape SEED_SHAPE_0 = new CyclicSeedShape("111111111111111");
    private static final CyclicSeedShape SEED_SHAPE_1 = new CyclicSeedShape("111101110111");
    private static final CyclicSeedShape SEED_SHAPE_2 = new CyclicSeedShape("111011010010111");

    // seed shapes combinations
    public static final CyclicSeedShape[] NO_SEED_SHAPE = {SEED_SHAPE_0};
    public static final CyclicSeedShape[] DEFAULT_SEED_SHAPES = {SEED_SHAPE_1, SEED_SHAPE_2};

    // aligning parameters
    public static ScoringMatrix SCORING_MATRIX = new ScoringMatrix("BLOSUM80", 11, 2);
    public static Integer FRAMESHIFT_PENALTY = 15;
    public static int X_BESTSCORE_DROP = 72;
    public static long REFERENCE_LETTERS, REFERENCE_SEQUENCES;

    public static void setScoringMatrix(String type, int gop, int gep) {
        SCORING_MATRIX = new ScoringMatrix(type, gop, gep);
    }

    public static void init(long referenceLetters, long referenceSequences) {
        REFERENCE_LETTERS = referenceLetters;
        REFERENCE_SEQUENCES = referenceSequences;
        BlastStatisticsHelper.init(SCORING_MATRIX, referenceLetters);
        double maxE = Math.pow(10, 18) / (2 * BlastStatisticsHelper.n * Math.pow(10, 6));
        X_BESTSCORE_DROP = (int) Math.round((Math.log(BlastStatisticsHelper.K * Math.pow(10, 18)) - Math.log(maxE)) / BlastStatisticsHelper.LAMBDA);
    }

}
