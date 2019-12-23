package ella.model.aligner.utils.gotoAligner;

import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.ScoringMatrix;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift;
import ella.model.io.MyParameters;

public class GotohAligner_Mod1 {

    private final static int MIN_VALUE = Integer.MIN_VALUE / 2;
    private final static int MAX_VALUE = Integer.MAX_VALUE / 2;
    private int gapOpen = MyParameters.SCORING_MATRIX.getGapOpen(), gapExtend = MyParameters.SCORING_MATRIX.getGapExtend();
    private ScoringMatrix scoringMatrix = MyParameters.SCORING_MATRIX;
    private int[] minArray;

    public Object[] run(String s1, String s2, int frame, Zheng_XDrop_Frameshift.AliMode aliMode) {

        String dna = s1;

        s1 = s1.replace('T', 'U');
        s1 = CodonTranslator.translate(s1, aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT ? frame : 1);
        if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT) {
            s1 = reverseSequence(s1);
            s2 = reverseSequence(s2);
        }

        if (s2.length() == 0 || s1.length() == 0)
            return cmpTrivialAlignment(s1, s2, frame, aliMode);

        int n1 = s1.length() + 1;
        int n2 = s2.length() + 1;

        int[][] D1 = new int[n1][n2];
        int[][] P1 = new int[n1][n2];
        int[][] Q1 = new int[n1][n2];

        minArray = new int[n2];
        for (int i = 0; i < minArray.length; i++)
            minArray[i] = MIN_VALUE;

        // initialization
        D1[0][0] = 0;
        P1[0][0] = MIN_VALUE;
        Q1[0][0] = MIN_VALUE;
        for (int i = 1; i < n1; i++) {
            D1[i][0] = -w(i);
            P1[i][0] = MIN_VALUE;
            Q1[i][0] = MIN_VALUE;
        }
        for (int j = 1; j < n2; j++) {
            D1[0][j] = -w(j);
            P1[0][j] = MIN_VALUE;
            Q1[0][j] = MIN_VALUE;
        }
        for (int i = 1; i < n1; i++) {
            for (int j = 1; j < n2; j++) {
                D1[i][j] = MIN_VALUE;
                P1[i][j] = MIN_VALUE;
                Q1[i][j] = MIN_VALUE;
            }
        }

        // xDrop alignment step
        int[] cell = {-1, -1, 0};
        int dropOutCounter = 0;
        int maxScore = MIN_VALUE;
        int k = 1, L = 1, U = 1;
        while (k < n1 + n2 - 1) {

            k = k + 1;
            dropOutCounter = L > U ? dropOutCounter + 1 : 0;
            if (dropOutCounter == 4)
                break;

            int L_Prime = -1, U_Prime = -2;
            for (int j = L; j < U + 1; j++) {

                int i = k - j;

                int match = add(D1[i - 1][j - 1], scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)));
                int pMax = Math.max(sub(D1[i - 1][j], w(1)), sub(P1[i - 1][j], gapExtend));
                int qMax = Math.max(sub(D1[i][j - 1], w(1)), sub(Q1[i][j - 1], gapExtend));
                int dMax = Math.max(match, Math.max(pMax, qMax));

                if (dMax > maxScore - MyParameters.X_BESTSCORE_DROP) {

                    P1[i][j] = pMax;
                    Q1[i][j] = qMax;
                    D1[i][j] = dMax;

                    if (dMax >= maxScore) {
                        maxScore = dMax;
                        if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT || aliMode == Zheng_XDrop_Frameshift.AliMode.RIGHT) {
                            cell[0] = i;
                            cell[1] = j;
                            cell[2] = dMax;
                        }
                    }
                    if (aliMode == Zheng_XDrop_Frameshift.AliMode.MIDDLE && (i == n1 - 1 && j == n2 - 1)) {
                        cell[0] = i;
                        cell[1] = j;
                        cell[2] = dMax;
                    }

                    L_Prime = L_Prime == -1 ? j : L_Prime;
                    U_Prime = j;

                }

            }

            if (U_Prime < L_Prime)
                break;
            L = Math.max(L_Prime, k + 2 - n1);
            U = Math.min(U_Prime + 1, n2 - 1);

        }

        if (cell == null || (cell[0] == -1 && cell[1] == -1)) {
            Object[] emptyResult = {"", "", "", 0, s1.length(), s2.length()};
            return emptyResult;
        }

        int aliScore = cell[2];
        int qRightNotAligned = n1 - (cell[0] + 1);
        int rRightNotAligned = n2 - (cell[1] + 1);

        StringBuilder a1 = new StringBuilder("");
        StringBuilder a2 = new StringBuilder("");
        StringBuilder[] ali = {a1, a2};
        while (cell != null)
            cell = backtrace(D1, P1, Q1, cell[0], cell[1], ali, s1, s2, aliMode);
        if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT) {
            ali[0].reverse();
            ali[1].reverse();
        }
        StringBuilder a3 = new StringBuilder(a1.length());
        for (int i = 0; i < a1.length(); i++) {
            if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT)
                a3.append(Math.abs(frame));
            else
                a3.append(1);
        }

        Object[] res = {ali[0].toString(), ali[1].toString(), a3.toString(), aliScore, qRightNotAligned, rRightNotAligned};

        return res;

    }

    private Object[] cmpTrivialAlignment(String s1, String s2, int frame, Zheng_XDrop_Frameshift.AliMode aliMode) {

        if (aliMode != Zheng_XDrop_Frameshift.AliMode.MIDDLE) {
            Object[] emptyResult = {"", "", "", 0, s1.length(), s2.length()};
            return emptyResult;
        }

        int length = Math.max(s1.length(), s2.length());
        int score = length == 0 ? 0 : -gapOpen - length * gapExtend;

        if (score < 0 - MyParameters.X_BESTSCORE_DROP) {
            Object[] emptyResult = {"", "", "", 0, s1.length(), s2.length()};
            return emptyResult;
        }

        String frameString = aliMode != Zheng_XDrop_Frameshift.AliMode.LEFT ? frameString(length, 1) : frameString(length, Math.abs(frame));
        Object[] res = {s1.length() == 0 ? gapString(s2.length()) : s1, s2.length() == 0 ? gapString(s1.length()) : s2, frameString, score, 0, 0};
        return res;

    }

    private String frameString(int l, int frame) {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < l; i++)
            buf.append(frame);
        return buf.toString();
    }

    private void fillArray(int[] a, int fillFrom, int fillTo) {
        System.arraycopy(minArray, 0, a, fillFrom, fillTo - fillFrom);
    }

    private String reverseSequence(String s) {
        StringBuilder rev = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; i--)
            rev.append(s.charAt(i));
        return rev.toString();
    }

    private int[] backtrace(int[][] D1, int[][] P1, int[][] Q1, int i, int j, StringBuilder[] alignment, String s1, String s2, Zheng_XDrop_Frameshift.AliMode aliMode) {

        // top-left corner reached
        if (i != 0 && j != 0) {

            // top border reached
            if (i == 0) {
                StringBuffer subSeq = new StringBuffer("");
                for (int k = j - 1; k >= 0; k--)
                    subSeq.insert(0, s2.charAt(k));
                alignment[0] = alignment[0].insert(0, gapString(j));
                alignment[1] = alignment[1].insert(0, subSeq);
                // backtrace(D1, P1, Q1, i, 0, alignment, s1, s2);
                int[] nextCell = {i, 0};
                return nextCell;
            }

            // left border reached
            else if (j == 0) {
                StringBuffer subSeq = new StringBuffer("");
                for (int k = i - 1; k >= 0; k--)
                    subSeq.insert(0, s1.charAt(k));
                alignment[0] = alignment[0].insert(0, subSeq);
                alignment[1] = alignment[1].insert(0, gapString(i));
                // backtrace(D1, P1, Q1, 0, j, alignment, s1, s2);
                int[] nextCell = {0, j};
                return nextCell;
            } else {

                // checking for match
                if (i != 0 && j != 0 && D1[i - 1][j - 1] + scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)) == D1[i][j]) {
                    alignment[0] = alignment[0].insert(0, s1.charAt(i - 1));
                    alignment[1] = alignment[1].insert(0, s2.charAt(j - 1));
                    // backtrace(D1, P1, Q1, i - 1, j - 1, alignment, s1, s2);
                    int[] nextCell = {i - 1, j - 1};
                    return nextCell;
                } else {

                    // checking left side
                    Integer l = null;
                    int score = D1[i][j];
                    for (int k = j - 1; k >= 0; k--) {
                        if (D1[i][k] - w(j - k) == D1[i][j]) {
                            l = k;
                            break;
                        }
                        score += gapExtend;
                        if (Q1[i][k] != score)
                            break;
                    }
                    if (l != null) {
                        StringBuffer subSeq = new StringBuffer("");
                        for (int k = j - 1; k >= l; k--)
                            subSeq.insert(0, s2.charAt(k));
                        alignment[0] = alignment[0].insert(0, gapString(j - l));
                        alignment[1] = alignment[1].insert(0, subSeq);
                        // backtrace(D1, P1, Q1, i, l, alignment, s1, s2);
                        int[] nextCell = {i, l};
                        return nextCell;
                    } else {

                        // checking top side
                        Integer t = null;
                        score = D1[i][j];
                        for (int k = i - 1; k >= 0; k--) {
                            if (D1[k][j] - w(i - k) == D1[i][j]) {
                                t = k;
                                break;
                            }
                            score += gapExtend;
                            if (P1[k][j] != score)
                                break;
                        }
                        if (t != null) {
                            StringBuffer subSeq = new StringBuffer("");
                            for (int k = i - 1; k >= t; k--)
                                subSeq.insert(0, s1.charAt(k));
                            alignment[0] = alignment[0].insert(0, subSeq);
                            alignment[1] = alignment[1].insert(0, gapString(i - t));
                            // backtrace(D1, P1, Q1, t, j, alignment, s1, s2);
                            int[] nextCell = {t, j};
                            return nextCell;
                        }

                    }
                }
            }

        }

        return null;

    }

    private Object gapString(Integer l) {
        StringBuffer buf = new StringBuffer();
        for (int i = 0; i < l; i++)
            buf.append("-");
        return buf.toString();
    }

    // taking care of overflow
    private int add(int a, int b) {
        if (a == MAX_VALUE || b == MAX_VALUE)
            return MAX_VALUE;
        return a + b;
    }

    // taking care of underflow
    private int sub(int a, int b) {
        if (a == MIN_VALUE || b == MIN_VALUE)
            return MIN_VALUE;
        return a - b;
    }

    private int w(int k) {
        return gapExtend * k + gapOpen;
    }

    // JUST FOR DEBUG
    // ***************************************************************

    private void printMatrix(int[][] m, String s1, String s2) {

        System.out.print("\t\t");
        for (int i = 0; i < s2.length(); i++)
            System.out.print(s2.charAt(i) + "\t");
        System.out.println();

        for (int i = 0; i < m.length; i++) {
            if (i > 0)
                System.out.print(s1.charAt(i - 1) + "\t");
            else
                System.out.print("\t");
            for (int j = 0; j < m[0].length; j++) {
                String num = m[i][j] == MIN_VALUE ? "MIN" : String.valueOf(m[i][j]);
                System.out.print(num + "\t");
            }
            System.out.println();
        }

    }

}
