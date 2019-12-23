package ella.model.aligner.utils.frameshiftAligner;

import java.util.ArrayList;
import java.util.Arrays;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift.AliMode;
import ella.model.io.MyParameters;

public class Aligner_Mod9 {

    private final static int MIN_VALUE = Integer.MIN_VALUE / 2;
    private String protein, dna;
    private String[] frameProteins;
    private int[] minArray;

    private boolean verbose = false;

    private int gop, gep;
    private Integer F;

    private ArrayList<Long> runtimes = new ArrayList<Long>();

    public Aligner_Mod9() {
        this.gop = MyParameters.SCORING_MATRIX.getGapOpen();
        this.gep = MyParameters.SCORING_MATRIX.getGapExtend();
        this.F = MyParameters.FRAMESHIFT_PENALTY;
    }

    public Object[] run(String s1, String s2, AliMode mode, DP_Matrix matrix, Aligner aligner) {

        this.dna = s1.replace('T', 'U');
        this.protein = s2;

        if (protein.length() == 0 || dna.length() <= 2)
            return cmpTrivialAlignment(mode);

        if (mode == AliMode.LEFT) {
            protein = reverseSequence(protein);
            dna = reverseSequence(dna);
        }

        frameProteins = new String[3];
        for (int i = 0; i < 3; i++)
            frameProteins[i] = CodonTranslator.translate(dna, i + 1);

        int n = protein.length() + 1;
        int m = dna.length() + 4;

        // initializing matrices
        minArray = new int[n];
        Arrays.setAll(minArray, i -> MIN_VALUE);
//        DP_Matrix matrix2 = new DP_Matrix();
        matrix.init(n, m);


        // filling Matrix
        int[] startingCell = fillingMatrix(n, m, matrix, mode, aligner);
        if (startingCell == null || (startingCell[0] == -1 && startingCell[1] == -1)) {
            Object[] emptyResult = {"", "", "", 0, dna.length(), protein.length()};
            return emptyResult;
        }

        // starting traceback
        Object[] tracebackResult = traceback(matrix, mode, startingCell, n, m);

        return tracebackResult;

    }

    private String reverseSequence(String s) {
        StringBuilder rev = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; i--)
            rev.append(s.charAt(i));
        return rev.toString();
    }

    private Object[] cmpTrivialAlignment(AliMode mode) {
        if (dna.length() < 3) {
            int score = protein.length() == 0 ? 0 : -gop - protein.length() * gep;
            Object[] res = {gapString(protein.length()), protein, frameString(protein.length(), 1), score, dna.length(), 0};
            return res;
        } else {
            StringBuffer p = new StringBuffer();
            for (int i = 0; i < dna.length() - 2; i += 3)
                p.append(CodonTranslator.translateCodon(dna.charAt(i), dna.charAt(i + 1), dna.charAt(i + 2)));
            Object[] res = {p, gapString(p.length()), frameString(protein.length(), 1), -gop - p.length() * gep, 0, 0};
            return res;
        }
    }

    private final int[] T_Primes = new int[4], T = new int[4], tracebackCell = new int[3];
    private final int[][] B = new int[4][2];

    private void initHelpArrays() {
        for (int i = 0; i < 4; i++) {
            T_Primes[i] = MIN_VALUE;
            T[i] = 0;
        }
        tracebackCell[0] = -1;
        tracebackCell[1] = -1;
        tracebackCell[2] = 0;
        for (int i = 0; i < 4; i++) {
            B[i][0] = 4;
            B[i][1] = 4 + i;
        }
    }

    private int[] fillingMatrix(int n, int m, DP_Matrix matrix, AliMode mode, Aligner aligner) {

        initHelpArrays();
//        int[] T_Primes = {MIN_VALUE, MIN_VALUE, MIN_VALUE, MIN_VALUE};
//        int[] tracebackCell = {-1, -1, 0};
//        int[] T = {0, 0, 0, 0};
//        int[][] B = {{4, 4}, {4, 5}, {4, 6}, {4, 7}};

        int k = 4, L = 4, U = 4;
        int dropOutCounter = 0;
        int maxScore = MIN_VALUE;

        while (k < m + n - 4) {

            // printMatrix(X, protein, dna);

            k = k + 1;
            int L_Prime = -1, U_Prime = -2;
            int[][] diagonals = getNextDiagonal(matrix, L, U);
            for (int i = 0; i < 3; i++) {
                diagonals[i][0] = U >= L ? L : m + 1;
                diagonals[i][1] = U >= L ? U - L + 3 : 2;
            }

            if (L > U) {
                if (dropOutCounter++ == 4)
                    return tracebackCell;
            } else {

                caseY(k, L, U, matrix, diagonals, aligner);
                caseD(k, L, U, matrix, diagonals, aligner);
                caseX(k, L, U, matrix, diagonals);

                int[] a = new int[U - L + 1];
                caseAFD(k, L, U, matrix, a, aligner);
                if (F != null) {
                    caseAF2(k, L, U, matrix, a, aligner);
                    caseAF4(k, L, U, matrix, a, aligner);
                }
                for (int j = L; j < U + 1; j++)
                    diagonals[2][j - L + 2] = max(diagonals[2][j - L + 2], a[j - L] + matchScore(k - j, j, false, mode));

                Object[] result = caseFinal(diagonals, U, L, maxScore, T_Primes, mode, tracebackCell, k, L_Prime, U_Prime, n, m, matrix);
                L_Prime = (int) result[0];
                U_Prime = (int) result[1];
                maxScore = (int) result[2];
//                tracebackCell = (int[]) result[3];
                diagonals = (int[][]) result[4];

            }

            int[] B_Prime = new int[2];
            B_Prime[0] = L_Prime;
            B_Prime[1] = L_Prime >= 0 ? U_Prime + 4 : U_Prime;
            B[(k - 1) & 3] = B_Prime;
            L = max(B[k & 3][0], k + 2 - n);
            U = min(B[k & 3][1], m - 3);

            T[k & 3] = T_Primes[k & 3];
            T_Primes[k & 3] = MIN_VALUE;

        }

        // printMatrix(X, protein, dna);
        // matrix.printMatrix(protein, dna);

        return tracebackCell;

    }

    private int[][] getNextDiagonal(DP_Matrix matrix, int L, int U) {
        return matrix.nextDiagonal(U - L + 3);
    }

    private final Object[] caseFinalResult = new Object[5];

    private Object[] caseFinal(int[][] diagonals, int U, int L, int maxScore, int[] T_Primes, AliMode mode, int[] tracebackCell, int k, int L_Prime,
                               int U_Prime, int n, int m, DP_Matrix matrix) {
        for (int j = L; j < U + 1; j++) {
            int xMax = diagonals[2][j - L + 2];
            maxScore = maxScore < xMax ? xMax : maxScore;

            if (xMax > maxScore - MyParameters.X_BESTSCORE_DROP) {

                if (xMax > T_Primes[k & 3]) {
                    T_Primes[k & 3] = xMax;
                    if (mode == AliMode.LEFT || mode == AliMode.RIGHT) {
                        if (xMax > tracebackCell[2]) {
//                            int[] cell = {k - j, j, xMax};
//                            tracebackCell = cell;
                            tracebackCell[0] = k - j;
                            tracebackCell[1] = j;
                            tracebackCell[2] = xMax;
                        }
                    } else {
                        if (k - j > tracebackCell[0] || (k - j == tracebackCell[0] && (xMax > tracebackCell[2] || j > tracebackCell[1] + 2))) {
//                            int[] cell = {k - j, j, xMax};
//                            tracebackCell = cell;
                            tracebackCell[0] = k - j;
                            tracebackCell[1] = j;
                            tracebackCell[2] = xMax;
                        }
                    }
                }

                L_Prime = L_Prime == -1 ? j : L_Prime;
                U_Prime = j;

            } else {
                diagonals[0][j - L + 2] = MIN_VALUE;
                diagonals[1][j - L + 2] = MIN_VALUE;
                diagonals[2][j - L + 2] = MIN_VALUE;
            }
        }
//        Object[] result = {L_Prime, U_Prime, maxScore, tracebackCell, diagonals};
//        return result;
        caseFinalResult[0] = L_Prime;
        caseFinalResult[1] = U_Prime;
        caseFinalResult[2] = maxScore;
        caseFinalResult[3] = tracebackCell;
        caseFinalResult[4] = diagonals;
        return caseFinalResult;
    }

    private void caseAFD(int k, int L, int U, DP_Matrix matrix, int[] a, Aligner aligner) {
        if (L > 6 && k - U - 1 > 0) {
            int[] X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 3);
            int x0 = X_diagonal[0], x1 = X_diagonal[1];
            fillArray(a, 0, min(max(L, X_diagonal[0] + 3), U + 1) - L);
            int fromIndex = max(L, X_diagonal[0] + 3), toIndex = min(U + 1, x0 + 3 + x1 - 2);
            if (toIndex > fromIndex)
                System.arraycopy(X_diagonal, fromIndex - x0 - 3 + 2, a, fromIndex - L, toIndex - fromIndex);
            fillArray(a, max(min(U + 1, x0 + 3 + x1 - 2), L) - L, U + 1 - L);
        } else {
            for (int j = L; j < U + 1; j++)
                a[j - L] = matrix.X(k - j - 1, j - 3);
        }
    }

    private void caseAF4(int k, int L, int U, DP_Matrix matrix, int[] a, Aligner aligner) {
        if (L > 7 && k - U - 1 > 0) {
            int[] X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 4);
            int x0 = X_diagonal[0], x1 = X_diagonal[1];
            fillArray(a, 0, min(max(L, x0 + 4), U + 1) - L);
            int fromIndex = max(L, x0 + 4), toIndex = min(U + 1, x0 + 4 + x1 - 2);
            for (int j = fromIndex; j < toIndex; j++)
                a[j - L] = max(a[j - L], X_diagonal[j - x0 - 4 + 2] - F);
            fillArray(a, max(min(U + 1, x0 + 4 + x1 - 2), L) - L, U + 1 - L);
        } else {
            for (int j = max(L, 8); j < U + 1; j++)
                a[j - L] = max(a[j - L], matrix.X(k - j - 1, j - 4) - F);
        }
    }

    private void caseAF2(int k, int L, int U, DP_Matrix matrix, int[] a, Aligner aligner) {
        if (L > 5 && k - U - 1 > 0) {
            int[] X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 2);
            int x0 = X_diagonal[0], x1 = X_diagonal[1];
            fillArray(a, 0, min(max(L, x0 + 2), U + 1) - L);
            int fromIndex = max(L, x0 + 2), toIndex = min(U + 1, x0 + 2 + x1 - 2);
            for (int j = fromIndex; j < toIndex; j++)
                a[j - L] = max(a[j - L], X_diagonal[j - x0 - 2 + 2] - F);
            fillArray(a, max(min(U + 1, x0 + 2 + x1 - 2), L) - L, U + 1 - L);
        } else {
            for (int j = max(L, 6); j < U + 1; j++)
                a[j - L] = max(a[j - L], matrix.X(k - j - 1, j - 2) - F);
        }
    }

    private void caseX(int k, int L, int U, DP_Matrix matrix, int[][] diagonals) {
        System.arraycopy(diagonals[1], 2, diagonals[2], 2, U - L + 1);
        for (int j = L; j < U + 1; j++)
            diagonals[2][j - L + 2] = max(diagonals[2][j - L + 2], diagonals[0][j - L + 2]);
    }

    private void caseD(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, Aligner aligner) {
        if (L > 6 && k - U > 0) {
            int[] X_diagonal = matrix.getDiagonals(2, k - L, L - 3);
            int x0 = X_diagonal[0], x1 = X_diagonal[1];
            fillArray(diagonals[0], 2, min(max(L, x0 + 3), U + 1) - L + 2);
            caseD1(k, L, U, matrix, diagonals, x0, x1);
            caseD2(k, L, U, matrix, diagonals, x0, x1);
            fillArray(diagonals[0], max(min(U + 1, x0 + 3 + x1 - 2), L) - L + 2, U + 1 - L + 2);
        } else {
            for (int j = L; j < U + 1; j++)
                diagonals[0][j - L + 2] = max(matrix.X(k - j, j - 3) - gop, matrix.D(k - j, j - 3) - gep);
        }
    }

    private void caseD1(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, int x0, int x1) {
        int[] X_diagonal = matrix.getDiagonals(2, k - L, L - 3);
        int fromIndex = max(L, x0 + 3), toIndex = min(U + 1, x0 + 3 + x1 - 2);
        for (int j = fromIndex; j < toIndex; j++)
            diagonals[0][j - L + 2] = X_diagonal[j - x0 - 3 + 2] - gop;
    }

    private void caseD2(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, int x0, int x1) {
        int[] D_diagonal = matrix.getDiagonals(0, k - L, L - 3);
        int d0 = D_diagonal[0];
        int fromIndex = max(L, x0 + 3), toIndex = min(U + 1, x0 + 3 + x1 - 2);
        for (int j = fromIndex; j < toIndex; j++)
            diagonals[0][j - L + 2] = max(diagonals[0][j - L + 2], D_diagonal[j - d0 - 3 + 2] - gep);
    }

    private void caseY(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, Aligner aligner) {
        if (L > 3 && k - U - 1 > 0) {
            int[] X_diagonal = matrix.getDiagonals(2, k - L - 1, L);
            int x0 = X_diagonal[0], x1 = X_diagonal[1];
            fillArray(diagonals[1], 2, min(max(L, x0), U + 1) - L + 2);
            caseY1(k, L, U, matrix, diagonals, x0, x1);
            caseY2(k, L, U, matrix, diagonals, x0, x1);
            fillArray(diagonals[1], max(min(U + 1, x0 + x1 - 2), L) - L + 2, U + 1 - L + 2);
        } else {
            for (int j = L; j < U + 1; j++)
                diagonals[1][j - L + 2] = max(matrix.X(k - j - 1, j) - gop, matrix.Y(k - j - 1, j) - gep);
        }
    }

    private void caseY1(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, int x0, int x1) {
        int[] X_diagonal = matrix.getDiagonals(2, k - L - 1, L);
        int fromIndex = max(L, X_diagonal[0]), toIndex = min(U + 1, x0 + x1 - 2);
        for (int j = fromIndex; j < toIndex; j++)
            diagonals[1][j - L + 2] = X_diagonal[j - x0 + 2] - gop;
    }

    private void caseY2(int k, int L, int U, DP_Matrix matrix, int[][] diagonals, int x0, int x1) {
        int[] Y_diagonal = matrix.getDiagonals(1, k - L - 1, L);
        int y0 = Y_diagonal[0];
        int fromIndex = max(L, x0), toIndex = min(U + 1, x0 + x1 - 2);
        for (int j = fromIndex; j < toIndex; j++)
            diagonals[1][j - L + 2] = max(diagonals[1][j - L + 2], Y_diagonal[j - y0 + 2] - gep);
    }

    private int min(int a, int b) {
        return a < b ? a : b;
    }

    private int max(int a, int b) {
        return a > b ? a : b;
    }

    private void fillArray(int[] a, int fillFrom, int fillTo) {
        System.arraycopy(minArray, 0, a, fillFrom, fillTo - fillFrom);
    }

    private Object[] traceback(DP_Matrix matrix, AliMode mode, int[] startingCell, int n, int m) {

        StringBuilder ali1 = new StringBuilder(n + m);
        StringBuilder ali2 = new StringBuilder(n + m);
        StringBuilder frames = new StringBuilder(n + m);
        StringBuilder[] alignment = {ali1, ali2, frames};

        int[] cell = startingCell;

        ArrayList<Integer> textIndices = new ArrayList<Integer>(ali1.length() * ali2.length());
        while (cell != null) {
            textIndices.add(cell[1]);
            cell = tracebackRec(matrix, cell[0], cell[1], alignment, mode);
        }

        // System.out.println("-----\n" + alignment[0] + "\n" + alignment[1] + "\n" + alignment[2]);

        if (mode == AliMode.LEFT) {
            alignment[0].reverse();
            alignment[1].reverse();
            alignment[2].reverse();
        }

        int qRightNotAligned = m - (startingCell[1] + 3);

        if (mode == AliMode.LEFT) {
            ArrayList<Integer> frameChanges = new ArrayList<Integer>(textIndices.size());
            int lastFrameChange = -1;
            for (int i = 0; i < textIndices.size(); i++) {
                int index = textIndices.get(i);
                int frameChange = ((m - index) % 3) + 1;
                if (frameChanges.isEmpty() || frameChange != lastFrameChange)
                    frameChanges.add(frameChange);
                lastFrameChange = frameChange;
            }
            StringBuilder shiftedFrames = new StringBuilder(alignment[2].length());
            int pointer = 0;
            char lastFrame = '-';
            for (int i = 0; i < alignment[2].length(); i++) {
                char frame = alignment[2].charAt(i);
                if (lastFrame != '-' && frame != lastFrame)
                    pointer++;
                lastFrame = frame;
                shiftedFrames.append(frameChanges.get(pointer) + "");
            }
            alignment[2] = shiftedFrames;
        }

        int qNotAligned = qRightNotAligned;
        qNotAligned -= Integer.parseInt(alignment[2].charAt(0) + "") - 1;
        qNotAligned /= 3;

        int rNotAligned = n - (startingCell[0] + 1);
        rNotAligned = rNotAligned < 0 ? 0 : rNotAligned;

        int aliScore = startingCell[2];
        Object[] res = {alignment[1].toString(), alignment[0].toString(), alignment[2].toString(), aliScore, qNotAligned, rNotAligned};

        return res;

    }

    private int[] tracebackRec(DP_Matrix matrix, int i, int j, StringBuilder[] alignment, AliMode mode) {

        if (verbose) {
            System.out.println(i + " " + j);
            System.out.println(alignment[0].toString());
            System.out.println(alignment[1].toString());
            System.out.println(alignment[2].toString());
        }

        if (i > 0 || j > 3) {

            int val = matrix.X(i, j);

            // checking for match
            if (i > 0 && j > 3) {

                int score = matchScore(i, j, false, mode);

                // checking for hit in same frame
                if ((j - 3) > 0 && matrix.X(i - 1, j - 3) == val - score) {
                    int[] prevCell = {i - 1, j - 3};
                    alignment[0].insert(0, protein.charAt(i - 1));
                    alignment[1].insert(0, getAA(j - 4, mode));
                    alignment[2].insert(0, ((j - 4) % 3) + 1);
                    return prevCell;
                }

                // checking for hit in different frames
                if (F != null) {
                    if ((j - 2) >= 2 && matrix.X(i - 1, j - 2) == val - score + F) {
                        int[] prevCell = {i - 1, j - 2};
                        alignment[0].insert(0, protein.charAt(i - 1));
                        alignment[1].insert(0, getAA(j - 4, mode));
                        alignment[2].insert(0, ((j - 4) % 3) + 1);
                        return prevCell;
                    }
                    if ((j - 4) >= 4 && matrix.X(i - 1, j - 4) == val - score + F) {
                        int[] prevCell = {i - 1, j - 4};
                        alignment[0].insert(0, protein.charAt(i - 1));
                        alignment[1].insert(0, getAA(j - 4, mode));
                        alignment[2].insert(0, ((j - 4) % 3) + 1);
                        return prevCell;
                    }
                }

            }

            // checking step-wise upper and left direction
            int gap_length = 0;
            while (true) {

                // checking upper direction
                Integer u = null;
                int pos = i - 1 - gap_length;
                if (pos >= 0) {
                    if (matrix.X(pos, j) - (gop + (i - 1 - pos) * gep) == val)
                        u = pos;
                    if (u != null) {
                        int[] prevCell = {u, j};
                        StringBuilder subSeq = new StringBuilder("");
                        for (int k = i; k > u; k--)
                            subSeq.insert(0, protein.charAt(k - 1));
                        alignment[0].insert(0, subSeq);
                        alignment[1].insert(0, gapString(i - u));
                        if (j > 3)
                            alignment[2].insert(0, frameString(i - u, ((j - 4) % 3) + 1));
                        else
                            alignment[2].insert(0, frameString(i - u, j));
                        return prevCell;
                    }
                }

                // checking left direction
                Integer t = null;
                pos = j - 3 - (gap_length * 3);
                if (pos >= 0) {
                    int c = (j - pos) / 3;
                    if (matrix.X(i, pos) - (gop + (c - 1) * gep) == val)
                        t = pos;
                    if (t != null) {
                        int[] prevCell = {i, t};
                        StringBuilder subSeq = new StringBuilder("");
                        for (int k = j; k > t; k -= 3)
                            subSeq.insert(0, getAA(k - 4, mode));
                        alignment[0].insert(0, gapString((j - t) / 3));
                        alignment[1].insert(0, subSeq);
                        alignment[2].insert(0, frameString((j - t) / 3, ((j - 4) % 3) + 1));
                        return prevCell;
                    }
                }

                gap_length++;

            }

        }

        return null;

    }

    private Object gapString(int l) {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < l; i++)
            buf.append("-");
        return buf.toString();
    }

    private Object frameString(int l, int frame) {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < l; i++)
            buf.append(frame);
        return buf.toString();
    }

    private int matchScore(int i, int j, boolean debug, AliMode mode) {
        char a = protein.charAt(i - 1);
        char b = mode != AliMode.LEFT ? frameProteins[(j - 4) % 3].charAt((j - 4) / 3)
                : CodonTranslator.translateCodon(dna.charAt(j - 2), dna.charAt(j - 3), dna.charAt(j - 4));
        return MyParameters.SCORING_MATRIX.getScore(a, b);
    }

    // private int matchScore(int i, int j, boolean debug, AliMode mode) {
    // char a = protein.charAt(i - 1);
    // char b = mode != AliMode.LEFT ? CodonTranslator_Array.translateCodon(dna.charAt(j - 4), dna.charAt(j - 3), dna.charAt(j - 2))
    // : CodonTranslator_Array.translateCodon(dna.charAt(j - 2), dna.charAt(j - 3), dna.charAt(j - 4));
    // return MyParameters.SCORING_MATRIX.getScore(a, b);
    // }

    private char getAA(int i, AliMode mode) {

        StringBuilder codon = new StringBuilder(3);
        for (int pos = i; pos < i + 3; pos++)
            codon.append(dna.charAt(pos));
        char aa = mode != AliMode.LEFT ? CodonTranslator.translateCodon(codon.toString())
                : CodonTranslator.translateCodon(codon.reverse().toString());

        return aa;
    }

    // *********************************************************************

    private void printResult(Object[] res) {

        String ref = (String) res[0];
        String query = (String) res[1];
        String frames = (String) res[2];
        int score = (int) res[3];
        int qAligned = (int) res[4];

        System.out.println(ref + "\n" + query + "\n" + frames + "\n" + score + "\n" + qAligned);

    }

}
