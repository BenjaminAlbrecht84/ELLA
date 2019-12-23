package ella.model.aligner.aligning.seeding;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.aligner.utils.MyFrame;
import ella.model.aligner.utils.ScoringMatrix;
import ella.model.aligner.utils.wrapper.ObjectArrayWrapper;
import ella.model.io.MyParameters;

public class MatchExtender {

    private static String query = "Escherichia-coli-str.-K-12-substr.-MG1655,-synthetic-genome_4174425_aligned_840_R_0_12711_0";
    private static String reference = "NP_418414.1";

    private final ScoringMatrix matrix = MyParameters.SCORING_MATRIX;
    private final int[] leftExtension = new int[2], rightExtension = new int[2];
    private final Object[] extension = new Object[2];

    public Object[] run(Object[] match, IndexText text, int qStart, String qSeq, MyFrame frame, int m, Aligner aligner, DiagonalTable diagonalTable) {

        int q = qStart;
        long r = (long) match[0];

        int MIN_SEED_SCORE = BlastStatisticsHelper.minSeedScore(qSeq.length());
//        int MIN_SEED_SCORE = BlastStatisticsHelper.minSeedScoreLast(m);
        int SEED_DROPOFF_SCORE = (int) Math.round((1. / BlastStatisticsHelper.LAMBDA) * 10.); // (int) Math.round(MIN_SEED_SCORE * 0.5);

        // // extending match
        doExtension(q, r, qSeq, text, SEED_DROPOFF_SCORE, aligner, diagonalTable);
        int maxScore = ((int[]) extension[0])[0] + ((int[]) extension[1])[0];

        if (maxScore < MIN_SEED_SCORE
                || !isOptimalExtension(q, r, qSeq, text, ((int[]) extension[0])[1], ((int[]) extension[1])[1], SEED_DROPOFF_SCORE))
            return null;

        aligner.seed_counter++;

        // extending match a second time now with a lower drop-off score, this provokes shorter but better matching seeds
        doExtension(q, r, qSeq, text, (int) Math.round(MIN_SEED_SCORE * 0.1), aligner, diagonalTable);
        maxScore = ((int[]) extension[0])[0] + ((int[]) extension[1])[0];

        // returning extended match as seed
        int queryStart = Math.abs(frame.getNumericalID()) - 1 + 3 * (q - ((int[]) extension[0])[1]);
        long indexStart = (long) match[0] - ((int[]) extension[0])[1];
        ObjectArrayWrapper info = new ObjectArrayWrapper(text.getProteinInfo((long) match[0]));
        int refStart = (int) (indexStart - (long) info.getData()[0]);
        Seed extSeed = new Seed(frame, queryStart, refStart, indexStart, ((int[]) extension[1])[1] + ((int[]) extension[0])[1], maxScore);

        Object[] result = {info, extSeed};
        return result;

    }

    private void doExtension(int q, long r, String qSeq, IndexText text, int SEED_DROPOFF_SCORE, Aligner aligner, DiagonalTable diagonalTable) {

        // extending match to the left
        extendToTheLeft(q, r, qSeq, text, SEED_DROPOFF_SCORE, aligner, diagonalTable);

        // extending match to the right
        extendToTheRight(q, r, qSeq, text, SEED_DROPOFF_SCORE, aligner, diagonalTable);

        extension[0] = leftExtension;
        extension[1] = rightExtension;
//        Object[] extension = {leftExtension, rightExtension};
//        return extension;

    }

    private void extendToTheLeft(int q, long r, String qSeq, IndexText text, int SEED_DROPOFF_SCORE, Aligner aligner, DiagonalTable diagonalTable) {

        int leftMaxScore = 0, leftScore = 0, left = 0, maxLeft = 0;
        while (q - left >= 0 && r - left >= 0 && leftMaxScore - leftScore < SEED_DROPOFF_SCORE) {

            // assessing alignment score
            char c1 = qSeq.charAt(q - left);
            char c2 = text.readAA(r - left);
            if (c2 == 127)
                break;
            int score = matrix.getScore(c1, c2);

            // processing score
            leftScore += score;
            if (leftScore < leftMaxScore - SEED_DROPOFF_SCORE)
                break;
            if (leftScore > leftMaxScore) {
                leftMaxScore = leftScore;
                maxLeft = left;
            }
            left++;

        }

        leftExtension[0] = leftMaxScore;
        leftExtension[1] = maxLeft;

//        int[] result = {leftMaxScore, maxLeft};
//        return result;

    }

    private void extendToTheRight(int q, long r, String qSeq, IndexText text, int SEED_DROPOFF_SCORE, Aligner aligner, DiagonalTable diagonalTable) {

        int rightMaxScore = 0, rightScore = 0, right = 1, maxRight = 1;
        while (r + right < text.length() && q + right < qSeq.length() && rightMaxScore - rightScore < SEED_DROPOFF_SCORE) {

            // assessing alignment score
            char c1 = qSeq.charAt(q + right);
            char c2 = text.readAA(r + right);
            if (c2 == 127)
                break;
            int score = matrix.getScore(c1, c2);

            // processing score
            rightScore += score;
            if (rightScore < rightMaxScore - SEED_DROPOFF_SCORE)
                break;
            if (rightScore > rightMaxScore) {
                rightMaxScore = rightScore;
                maxRight = right;
            }
            right++;

        }

        rightExtension[0] = rightMaxScore;
        rightExtension[1] = maxRight;

//        int[] result = {rightMaxScore, maxRight};
//        return result;

    }

    private boolean isOptimalExtension(int q, long r, String qSeq, IndexText text, int maxLeft, int maxRight, int SEED_DROPOFF_SCORE) {

        // checking for negative prefixes/suffixes and for segments with score smaller than -SEED_DROPOFF_SCORE
        int prefixScore = 0, maxScore = 0;
        for (int i = maxLeft; i >= 0; i--) {
            char c1 = qSeq.charAt(q - i);
            char c2 = text.readAA(r - i);
            prefixScore += matrix.getScore(c1, c2);
            if (prefixScore > maxScore)
                maxScore = prefixScore;
            else if (prefixScore < maxScore - SEED_DROPOFF_SCORE || prefixScore <= 0) {
                return false;
            }
        }
        for (int i = 1; i <= maxRight; i++) {
            char c1 = qSeq.charAt(q + i);
            char c2 = text.readAA(r + i);
            prefixScore += matrix.getScore(c1, c2);
            if (prefixScore > maxScore)
                maxScore = prefixScore;
            else if (prefixScore < maxScore - SEED_DROPOFF_SCORE || prefixScore <= 0 || i == maxRight) {
                return false;
            }
        }

        return true;
    }

}
