package ella.model.aligner.utils.gotoAligner;

import ella.model.io.MyParameters;

import java.util.ArrayList;

public class GotohMatrix {

    private int MIN_VALUE = Integer.MIN_VALUE / 2;

    private int pointer = 0, offset = 0;
    private ArrayList<int[][]> diagonals = new ArrayList<>();
    private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

    public void init(int n, int m) {
        pointer = 0;
        diagonals.ensureCapacity(n + m);
    }

    public int[][] nextDiagonal(int size) {
        if (pointer == diagonals.size())
            diagonals.add(new int[3][size > 2 ? size + offset : 2]);
        return diagonals.get(pointer++);
    }

    public void addDiagonal(int[][] d) {
        diagonals.add(d);
    }

    public int D(int i, int j) {
        return D(i, j, null);
    }

    public int P(int i, int j) {
        return P(i, j, null);
    }

    public int Q(int i, int j) {
        return Q(i, j, null);
    }

    public int D(int i, int j, int[] D) {

        // handle first cell
        if (i == 0 && j == 0)
            return 0;

        // handling 1st row
        if (i == 0)
            return -gop - j * gep;

        // handling 1st column
        if (j == 0)
            return -gop - i * gep;

        // handling diagonals
        if (D != null)
            return diagonal(i, j, D);
        return diagonal(0, i, j);
    }

    public int P(int i, int j, int[] P) {

        // handling first column/row
        if (i == 0 || j == 0)
            return MIN_VALUE;

        // handling diagonals
        if (P != null)
            return diagonal(i, j, P);
        return diagonal(1, i, j);

    }

    public int Q(int i, int j, int[] Q) {

        // handling first column/row
        if (i == 0 || j == 0)
            return MIN_VALUE;

        // handling diagonals
        if (Q != null)
            return diagonal(i, j, Q);
        return diagonal(2, i, j);

    }

    public int[] getDiagonals(int type, int i, int j) {
        int d = i + j - 5;
        if (d >= 0 && d < diagonals.size())
            return diagonals.get(d)[type];
        return null;
    }

    public int diagonal(int i, int j, int[] a) {
        int k = j - a[0] + 2;
        if (k >= 2 && k < a[1])
            return a[k];
        return MIN_VALUE;
    }

    public int diagonal(int type, int i, int j) {
        int d = i + j - 2;
        if (d < diagonals.size()) {
            int[] diagonal = diagonals.get(d)[type];
            int k = j - diagonal[0] + 2;
            if (k >= 2 && k < diagonal[1])
                return diagonal[k];
        }
        return MIN_VALUE;
    }

    public void freeMemory() {
        diagonals = null;
    }

    // ONLY FOR DEBUGGING ***********************************

    public void printMatrix(String s1, String s2, int n, int m) {

        System.out.print("\t\t");
        for (int i = 0; i < s2.length(); i++)
            System.out.print(i + "\t");
        System.out.println();

        System.out.print("\t\t\t");
        for (int i = 0; i < s2.length(); i++)
            System.out.print(s2.charAt(i) + "\t");
        System.out.println();

        for (int i = 0; i < n; i++) {
            System.out.print(i + "\t");
            if (i > 0)
                System.out.print(s1.charAt(i - 1) + "\t");
            else
                System.out.print("\t");
            for (int j = 0; j < m; j++) {
                String val = D(i, j) == MIN_VALUE ? "MIN" : String.valueOf(D(i, j));
                System.out.print(val + "\t");
            }
            System.out.println();
        }

    }

}
