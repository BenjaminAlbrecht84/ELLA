package ella.model.aligner.utils.frameshiftAligner;

import java.util.ArrayList;

import ella.model.io.MyParameters;

public class DP_Matrix {

    private int MIN_VALUE = Integer.MIN_VALUE / 2;

    private int pointer = 0, offset = 0;
    private ArrayList<int[][]> diagonals = new ArrayList<>();
    private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

    public void init(int n, int m) {
        pointer = 0;
        diagonals.ensureCapacity(n + m);
    }

    public int[][] nextDiagonal(int size) {
        int s = size > 2 ? size + offset : 2;
        if (pointer == diagonals.size())
            diagonals.add(new int[3][s]);
        else if (diagonals.get(pointer)[0].length < s) {
            diagonals.remove(pointer);
            diagonals.add(pointer, new int[3][s]);
        }
        return diagonals.get(pointer++);
    }

    public void addDiagonal(int[][] d) {
        diagonals.add(d);
    }

    public int D(int i, int j) {
        return D(i, j, null);
    }

    public int Y(int i, int j) {
        return Y(i, j, null);
    }

    public int X(int i, int j) {
        return X(i, j, null);
    }

    public int D(int i, int j, int[] D) {

        // handling 1st row
        if (i == 0) {
            if (j <= 3)
                return 0;
            else
                return -gop - (((j - 1) / 3) - 1) * gep;
        }

        // handling 1st column
        if (i > 0 && j == 0)
            return MIN_VALUE;

        // handling column 2-4
        if (j <= 3)
            return -gop - ((i - 1) * gep) - gop;

        // handling diagonals
        if (D != null)
            return diagonal(2, i, j, D);
        return diagonal(0, i, j);
    }

    public int Y(int i, int j, int[] Y) {

        // handling first row
        if (i == 0) {
            if (j <= 3)
                return 0;
            else
                return MIN_VALUE;
        }

        // handling 1st column
        if (i > 0 && j == 0)
            return MIN_VALUE;

        // handling column 2-4
        if (j <= 3)
            return -gop - ((i - 1) * gep);

        // handling diagonals
        if (Y != null)
            return diagonal(2, i, j, Y);
        return diagonal(1, i, j);

    }

    public int X(int i, int j, int[] X) {

        // handling first row
        if (i == 0) {
            if (j <= 3)
                return 0;
            else
                return -gop - (((j - 1) / 3) - 1) * gep;
        }

        // handling 1st column
        if (i > 0 && j == 0)
            return MIN_VALUE;

        // handling column 2-4
        if (j <= 3)
            return -gop - ((i - 1) * gep);

        // handling diagonals
        if (X != null)
            return diagonal(2, i, j, X);
        return diagonal(2, i, j);

    }

    public int[] getDiagonals(int type, int i, int j) {
        int d = i + j - 5;
        if (d >= 0 && d < diagonals.size())
            return diagonals.get(d)[type];
        return null;
    }

    public int diagonal(int type, int i, int j, int[] a) {
        int k = j - a[0] + 2;
        if (k >= 2 && k < a[1])
            return a[k];
        return MIN_VALUE;
    }

    public int diagonal(int type, int i, int j) {
        int d = i + j - 5;
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
        for (int i = 0; i < s2.length() + 4; i++)
            System.out.print(i + "\t");
        System.out.println();

        System.out.print("\t\t\t\t\t\t");
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
                String val = X(i, j) == MIN_VALUE ? "-IN" : String.valueOf(X(i, j));
                System.out.print(val + "\t");
            }
            System.out.println();
        }

    }

}
