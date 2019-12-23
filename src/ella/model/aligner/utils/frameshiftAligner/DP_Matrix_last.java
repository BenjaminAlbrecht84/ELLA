package ella.model.aligner.utils.frameshiftAligner;

import ella.model.aligner.utils.streams.MyDiagonalStream;
import ella.model.io.MyParameters;

public class DP_Matrix_last {

	private int MIN_VALUE = Integer.MIN_VALUE / 2;

	private MyDiagonalStream D = new MyDiagonalStream(), Y = new MyDiagonalStream(), X = new MyDiagonalStream(), pointers = new MyDiagonalStream();
	private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

	public void init(int n, int m) {
		int capacity = n * m + (n + m) * 2;
		D.ensureCapacity(capacity);
		X.ensureCapacity(capacity);
		Y.ensureCapacity(capacity);
		pointers.ensureCapacity(n + m);
		D.reset();
		X.reset();
		Y.reset();
		pointers.reset();
	}

	public void addDiagonal(int[] d, int[] y, int[] x) {
		pointers.add(X.size());
		X.add(x);
		Y.add(y);
		D.add(d);
	}

	private MyDiagonalStream getStream(int type) {
		return type == 0 ? D : type == 1 ? Y : X;
	}

	public int D(int i, int j) {
		return D(i, j, null);
	}

	public int D(int i, int j, int[] diagonal) {

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
		if (diagonal != null)
			return diagonal(diagonal, i, j);
		return diagonal(0, i, j);
	}

	public int Y(int i, int j) {
		return Y(i, j, null);
	}

	public int Y(int i, int j, int[] diagonal) {

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
		if (diagonal != null)
			return diagonal(diagonal, i, j);
		return diagonal(1, i, j);

	}

	public int X(int i, int j) {
		return X(i, j, null);
	}

	public int X(int i, int j, int[] diagonal) {

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
		if (diagonal != null)
			return diagonal(diagonal, i, j);
		return diagonal(2, i, j);

	}

	public int diagonal(int[] diagonal, int i, int j) {
		int k = j - diagonal[0] + 2;
		if (k >= 2 && k < diagonal[1])
			return diagonal[k];
		return MIN_VALUE;
	}

	public int diagonal(int type, int i, int j) {
		int d = i + j - 5;
		if (d < pointers.size()) {
			MyDiagonalStream a = getStream(type);
			int p = pointers.get(d);
			int k = j - a.get(p) + 2;
			if (k >= 2 && k < a.get(p + 1))
				return a.get(p + k);
		}
		return MIN_VALUE;
	}

	public void freeMemory() {
		D = Y = X = pointers = null;
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
				int val = X(i, j, null) == MIN_VALUE ? -1000 : X(i, j, null);
				System.out.print(val + "\t");
			}
			System.out.println();
		}

	}

}
