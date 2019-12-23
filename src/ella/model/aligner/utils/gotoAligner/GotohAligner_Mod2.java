package ella.model.aligner.utils.gotoAligner;

import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.ScoringMatrix;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift;
import ella.model.io.MyParameters;

import java.util.ArrayList;

public class GotohAligner_Mod2 {

    private final static int MIN_VALUE = Integer.MIN_VALUE / 2;
    private final static int MAX_VALUE = Integer.MAX_VALUE / 2;
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

        if (s2.length() == 0 || s1.length() == 0)
            return cmpTrivialAlignment(s1, s2, frame, aliMode);

        int n1 = s1.length() + 1;
        int n2 = s2.length() + 1;

        // xDrop alignment step
        GotohMatrix matrix = new GotohMatrix();
        matrix.init(n1, n2);
        int[] cell = runAlignment(matrix, s1, s2, n1, n2, aliMode);

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
            cell = backtrace(matrix, cell[0], cell[1], ali, s1, s2, aliMode);
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

    private int[] runAlignment(GotohMatrix matrix, String s1, String s2, int n1, int n2, Zheng_XDrop_Frameshift.AliMode aliMode) {

        // init diagonals
        ArrayList<int[][]> prevDiagonals = new ArrayList<>();
        prevDiagonals.add(null);
        prevDiagonals.add(null);

        int[] cell = {-1, -1, 0};
        int dropOutCounter = 0;
        int maxScore = MIN_VALUE;
        int k = 1, L = 1, U = 1;
        while (k < n1 + n2 - 2) {

            k = k + 1;
            dropOutCounter = L > U ? dropOutCounter + 1 : 0;
            if (dropOutCounter == 4)
                break;

            int[][] diagonals = matrix.nextDiagonal(U - L + 3);
            for (int i = 0; i < 3; i++) {
                diagonals[i][0] = U >= L ? L : n2 + 1;
                diagonals[i][1] = U >= L ? U - L + 3 : 2;
            }

            getMaximum(k, L, U, s1, s2, matrix, prevDiagonals, diagonals);
            int L_Prime = -1, U_Prime = -2;
            for (int j = L; j < U + 1; j++) {
                int i = k - j;
                int recMax = max(max(diagonals[0][j - L + 2], diagonals[1][j - L + 2]), diagonals[2][j - L + 2]);
                int[] recResult = storeMaximum(recMax, maxScore, diagonals, cell, aliMode, i, j, L, n1, n2, L_Prime, U_Prime);
                if (recResult != null) {
                    L_Prime = recResult[0];
                    U_Prime = recResult[1];
                    maxScore = recResult[2];
                }
            }

            if (U_Prime < L_Prime)
                break;
            L = max(L_Prime, k + 2 - n1);
            U = min(U_Prime + 1, n2 - 1);
            prevDiagonals.remove(0);
            prevDiagonals.add(diagonals);

        }

//        matrix.printMatrix(s1, s2, n1, n2);

        return cell;

    }

    private void getMaximum(int k, int L, int U, String s1, String s2, GotohMatrix matrix, ArrayList<int[][]> prevDiagonals, int[][] diagonals) {
        int[] D2 = prevDiagonals.get(0) != null ? prevDiagonals.get(0)[0] : null;
        int[] D1 = prevDiagonals.get(1) != null ? prevDiagonals.get(1)[0] : null;
        int[] P1 = prevDiagonals.get(1) != null ? prevDiagonals.get(1)[1] : null;
        int[] Q1 = prevDiagonals.get(1) != null ? prevDiagonals.get(1)[2] : null;
        caseD(L, U, k, diagonals, s1, s2, D2, matrix);
        caseP(L, U, k, diagonals, D1, P1, matrix);
        caseQ(L, U, k, diagonals, D1, Q1, matrix);
    }

    private void caseQ(int L, int U, int k, int[][] diagonals, int[] D1, int[] Q1, GotohMatrix matrix) {
        if (D1 != null && Q1 != null) {
            int[] q = new int[U - L + 3];
            int fromIndex = max(L, D1[0] + 1);
            int toIndex = min(U + 1, D1[0] + D1[1] - 2 + 1);
            for (int j = L; j < fromIndex; j++)
                q[j - L + 2] = max(sub(matrix.D(k - j, j - 1, D1), w(1)), sub(matrix.Q(k - j, j - 1, Q1), gapExtend));
            for (int j = fromIndex; j < toIndex; j++)
                q[j - L + 2] = D1[j - D1[0] + 2 - 1] - (gapExtend + gapOpen);
            for (int j = fromIndex; j < toIndex; j++)
                q[j - L + 2] = max(q[j - L + 2], Q1[j - Q1[0] + 2 - 1] - gapExtend);
            for (int j = toIndex; j < U + 1; j++)
                q[j - L + 2] = max(sub(matrix.D(k - j, j - 1, D1), w(1)), sub(matrix.Q(k - j, j - 1, Q1), gapExtend));
            System.arraycopy(q, 2, diagonals[2], 2, q.length - 2);
        } else {
            for (int j = L; j < U + 1; j++)
                diagonals[2][j - L + 2] = max(sub(matrix.D(k - j, j - 1, D1), w(1)), sub(matrix.Q(k - j, j - 1, Q1), gapExtend));
        }
    }

    private void caseP(int L, int U, int k, int[][] diagonals, int[] D1, int[] P1, GotohMatrix matrix) {
        if (D1 != null && P1 != null) {
            int[] p = new int[U - L + 3];
            int fromIndex = max(L, D1[0]);
            int toIndex = min(U + 1, D1[0] + D1[1] - 2);
            for (int j = L; j < fromIndex; j++)
                p[j - L + 2] = max(sub(matrix.D(k - j - 1, j, D1), w(1)), sub(matrix.P(k - j - 1, j, P1), gapExtend));
            for (int j = fromIndex; j < toIndex; j++)
                p[j - L + 2] = D1[j - D1[0] + 2] - (gapExtend + gapOpen);
            for (int j = fromIndex; j < toIndex; j++)
                p[j - L + 2] = max(p[j - L + 2], P1[j - P1[0] + 2] - gapExtend);
            for (int j = toIndex; j < U + 1; j++)
                p[j - L + 2] = max(sub(matrix.D(k - j - 1, j, D1), w(1)), sub(matrix.P(k - j - 1, j, P1), gapExtend));
            System.arraycopy(p, 2, diagonals[1], 2, p.length - 2);
        } else {
            for (int j = L; j < U + 1; j++)
                diagonals[1][j - L + 2] = max(sub(matrix.D(k - j - 1, j, D1), w(1)), sub(matrix.P(k - j - 1, j, P1), gapExtend));
        }
    }

    private void caseD(int L, int U, int k, int[][] diagonals, String s1, String s2, int[] D2, GotohMatrix matrix) {
        if (D2 != null) {
            int[] d = new int[U - L + 3];
            int fromIndex = max(L, D2[0] + 1);
            int toIndex = min(U + 1, D2[0] + D2[1] - 2 + 1);
            for (int j = L; j < fromIndex; j++)
                d[j - L + 2] = add(matrix.D(k - j - 1, j - 1, D2), scoringMatrix.getScore(s1.charAt(k - j - 1), s2.charAt(j - 1)));
            int[] m = new int[toIndex - fromIndex];
            for (int j = fromIndex; j < toIndex; j++)
                m[j - fromIndex] = scoringMatrix.getScore(s1.charAt(k - j - 1), s2.charAt(j - 1));
            for (int j = fromIndex; j < toIndex; j++)
                d[j - L + 2] = D2[j - D2[0] + 2 - 1] + m[j - fromIndex];
            for (int j = toIndex; j < U + 1; j++)
                d[j - L + 2] = add(matrix.D(k - j - 1, j - 1, D2), scoringMatrix.getScore(s1.charAt(k - j - 1), s2.charAt(j - 1)));
            System.arraycopy(d, 2, diagonals[0], 2, d.length - 2);
        } else {
            for (int j = L; j < U + 1; j++)
                diagonals[0][j - L + 2] = add(matrix.D(k - j - 1, j - 1, D2), scoringMatrix.getScore(s1.charAt(k - j - 1), s2.charAt(j - 1)));
        }
    }

    private int max(int a, int b) {
        return a > b ? a : b;
    }

    private int min(int a, int b) {
        return a < b ? a : b;
    }

    private int[] storeMaximum(int recMax, int maxScore, int[][] diagonals, int[] cell, Zheng_XDrop_Frameshift.
            AliMode aliMode, int i, int j, int L, int n1, int n2, int L_Prime, int U_Prime) {
        if (recMax > maxScore - MyParameters.X_BESTSCORE_DROP) {

            if (recMax >= maxScore) {
                maxScore = recMax;
                if (aliMode == Zheng_XDrop_Frameshift.AliMode.LEFT || aliMode == Zheng_XDrop_Frameshift.AliMode.RIGHT) {
                    cell[0] = i;
                    cell[1] = j;
                    cell[2] = recMax;
                }
            }
            if (aliMode == Zheng_XDrop_Frameshift.AliMode.MIDDLE && (i == n1 - 1 && j == n2 - 1)) {
                cell[0] = i;
                cell[1] = j;
                cell[2] = recMax;
            }

            L_Prime = L_Prime == -1 ? j : L_Prime;
            U_Prime = j;

            int[] result = {L_Prime, U_Prime, maxScore};
            return result;

        }

        return null;

    }

    private Object[] cmpTrivialAlignment(String s1, String s2, int frame, Zheng_XDrop_Frameshift.AliMode aliMode) {

        if (aliMode != Zheng_XDrop_Frameshift.AliMode.MIDDLE) {
            Object[] emptyResult = {"", "", "", 0, s1.length(), s2.length()};
            return emptyResult;
        }

        int length = max(s1.length(), s2.length());
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

    private int[] backtrace(GotohMatrix matrix, int i, int j, StringBuilder[] alignment, String s1, String
            s2, Zheng_XDrop_Frameshift.AliMode aliMode) {

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
                if (i != 0 && j != 0 && matrix.D(i - 1, j - 1) + scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)) == matrix.D(i, j)) {
                    alignment[0] = alignment[0].insert(0, s1.charAt(i - 1));
                    alignment[1] = alignment[1].insert(0, s2.charAt(j - 1));
                    int[] nextCell = {i - 1, j - 1};
                    return nextCell;
                } else {

                    // checking left side
                    Integer l = null;
                    int score = matrix.D(i, j);
                    for (int k = j - 1; k >= 0; k--) {
                        if (matrix.D(i, k) - w(j - k) == matrix.D(i, j)) {
                            l = k;
                            break;
                        }
                        score += gapExtend;
                        if (matrix.Q(i, k) != score)
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
                        score = matrix.D(i, j);
                        for (int k = i - 1; k >= 0; k--) {
                            if (matrix.D(k, j) - w(i - k) == matrix.D(i, j)) {
                                t = k;
                                break;
                            }
                            score += gapExtend;
                            if (matrix.P(k, j) != score)
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
