package ella.model.aligner.utils.gotoAligner;

import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.ScoringMatrix;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift;
import ella.model.io.MyParameters;

public class GotohAligner {

    private int gapOpen = MyParameters.SCORING_MATRIX.getGapOpen(), gapExtend = MyParameters.SCORING_MATRIX.getGapExtend();
    private ScoringMatrix scoringMatrix = MyParameters.SCORING_MATRIX;
    private int[] minArray;

    public Object[] run(String s1, String s2, int frame, Zheng_XDrop_Frameshift.AliMode aliMode) {

        s1 = s1.replace('T', 'U');
        s1 = CodonTranslator.translate(s1, aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT ? frame : 1);
        if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT) {
            s1 = reverseSequence(s1);
            s2 = reverseSequence(s2);
        }

        int n1 = s1.length() + 1;
        int n2 = s2.length() + 1;

        int[][] D1 = new int[n1][n2];
        int[][] P1 = new int[n1][n2];
        int[][] Q1 = new int[n1][n2];

        minArray = new int[n2];
        for (int i = 0; i < minArray.length; i++)
            minArray[i] = Integer.MIN_VALUE;

        // initialization
        D1[0][0] = 0;
        P1[0][0] = Integer.MIN_VALUE;
        Q1[0][0] = Integer.MIN_VALUE;
        for (int i = 1; i < n1; i++) {
            D1[i][0] = -w(i);
            P1[i][0] = Integer.MIN_VALUE;
            Q1[i][0] = Integer.MIN_VALUE;
        }
        for (int j = 1; j < n2; j++) {
            D1[0][j] = -w(j);
            P1[0][j] = Integer.MIN_VALUE;
            Q1[0][j] = Integer.MIN_VALUE;
        }

        // recursion
        for (int i = 1; i < n1; i++) {
            for (int j = 1; j < n2; j++) {
                int match = add(D1[i - 1][j - 1], scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)));
                P1[i][j] = Math.max(sub(D1[i - 1][j], w(1)), sub(P1[i - 1][j], gapExtend));
                Q1[i][j] = Math.max(sub(D1[i][j - 1], w(1)), sub(Q1[i][j - 1], gapExtend));
                D1[i][j] = Math.max(match, Math.max(P1[i][j], Q1[i][j]));
            }
        }

        int[] cell = {n1 - 1, n2 - 1};
        if (aliMode != Zheng_XDrop_Frameshift.AliMode.MIDDLE)
            cell = maxBorderEntry(D1, n1 - 1, n2 - 1);
        int aliScore = D1[cell[0]][cell[1]];
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

    private void fillArray(int[] a, int fillFrom, int fillTo) {
        System.arraycopy(minArray, 0, a, fillFrom, fillTo - fillFrom);
    }

    private String reverseSequence(String s) {
        StringBuilder rev = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; i--)
            rev.append(s.charAt(i));
        return rev.toString();
    }

    private int[] maxBorderEntry(int[][] M, int n, int m) {
        int max = Integer.MIN_VALUE;
        int[] cell = {0, 0};

        // searching last column
        for (int i = 0; i < M.length; i++) {
            if (M[i][m] > max) {
                max = M[i][m];
                cell[0] = i;
                cell[1] = m;
            }
        }

        // searching last row
        for (int j = 0; j < M[0].length; j++) {
            if (M[n][j] > max) {
                max = M[n][j];
                cell[0] = n;
                cell[1] = j;
            }
        }
        return cell;
    }

    private int[] backtrace(int[][] D1, int[][] P1, int[][] Q1, int i, int j, StringBuilder[] alignment, String s1, String s2, Zheng_XDrop_Frameshift.AliMode aliMode) {

        // System.out.println("\n["+i+","+j+"]: "+D1[i][j]);

        // top-left corner reached
        if (i != 0 && j != 0) {

            // top border reached
            if (i == 0) {
                StringBuffer subSeq = new StringBuffer("");
                for (int k = j - 1; k >= 0; k--)
                    subSeq.insert(0, s2.charAt(k));
                alignment[0] = alignment[0].insert(0, gapString(j));
                alignment[1] = alignment[1].insert(0, subSeq);
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
                int[] nextCell = {0, j};
                return nextCell;
            } else {

                // checking for match
                if (i != 0 && j != 0 && D1[i - 1][j - 1] + scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)) == D1[i][j]) {
                    alignment[0] = alignment[0].insert(0, s1.charAt(i - 1));
                    alignment[1] = alignment[1].insert(0, s2.charAt(j - 1));
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
        if (a == Integer.MAX_VALUE || b == Integer.MAX_VALUE)
            return Integer.MAX_VALUE;
        return a + b;
    }

    // taking care of underflow
    private int sub(int a, int b) {
        if (a == Integer.MIN_VALUE || b == Integer.MIN_VALUE)
            return Integer.MIN_VALUE;
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
                System.out.print(m[i][j] + "\t");
            }
            System.out.println();
        }

    }

}
