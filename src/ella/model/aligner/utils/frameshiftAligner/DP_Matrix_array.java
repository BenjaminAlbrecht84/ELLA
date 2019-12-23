package ella.model.aligner.utils.frameshiftAligner;

import ella.model.io.MyParameters;

public class DP_Matrix_array {

	private int MIN_VALUE = Integer.MIN_VALUE / 2;

	private int pointer = 0, offset = 10;
	private int[][][] diagonals;
	private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

	public void init(int n, int m) {
		pointer = 0;
		diagonals = new int[n * m][][];
	}

	public int[][] nextDiagonal(int size) {
		if (diagonals[pointer] == null)
			diagonals[pointer] = new int[3][size > 2 ? size + offset : 2];
		ensureSize(size);
		return diagonals[pointer++];
	}

	private void ensureSize(int size) {
		if (size > diagonals[pointer][0].length) {
			int[][] diagonal = new int[3][size + offset];
			diagonals[pointer] = diagonal;
		}
	}

	public int D(int i, int j) {

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
		return diagonal(0, i, j);
	}

	public int Y(int i, int j) {

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
		return diagonal(1, i, j);

	}

	public int X(int i, int j) {

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
		return diagonal(2, i, j);

	}

	public int[] getDiagonals(int type, int i, int j) {
		int d = i + j - 5;
		if (d >= 0 && d < diagonals.length)
			return diagonals[d][type];
		return null;
	}

	public int diagonal(int type, int i, int j) {
		int d = i + j - 5;
		if (d < diagonals.length) {
			int[] diagonal = diagonals[d][type];
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
				int val = X(i, j) == MIN_VALUE ? -1000 : X(i, j);
				System.out.print(val + "\t");
			}
			System.out.println();
		}

	}

}
