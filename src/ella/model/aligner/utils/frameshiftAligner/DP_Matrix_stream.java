package ella.model.aligner.utils.frameshiftAligner;

import ella.model.aligner.utils.streams.MyDiagonalStream;
import ella.model.io.MyParameters;

public class DP_Matrix_stream {

	private int MIN_VALUE = Integer.MIN_VALUE / 2;

	private MyDiagonalStream positions = new MyDiagonalStream();
	private MyDiagonalStream[] diagonals = new MyDiagonalStream[3];
	private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

	public void init(int n, int m) {
		int capacity = (int) Math.round((double) n * (double) m * (1.));
		for (int i = 0; i < 3; i++)
			diagonals[i] = new MyDiagonalStream(capacity);
	}

	public void newDiagonal() {
		positions.add(diagonals[0].size());
	}

	public void addValue(int type, int value) {
		diagonals[type].add(value);
	}

	public int getPosition(int i, int j) {
		int d = i + j - 5;
		return positions.get(d);
	}

	public int getValue(int type, int i, int j, int index) {
		int d = i + j - 5;
		return diagonals[type].get(positions.get(d) + index);
	}

	public int getValue(int type, int index) {
		return diagonals[type].get(index);
	}

	public void setValue(int type, int index, int value) {
		diagonals[type].set(index, value);
	}

	public void replace(int type, int index, int value) {
		if (value > diagonals[type].get(index))
			diagonals[type].set(index, value);
	}

	public void fill(int type, int fromIndex, int toIndex, int value) {
		for (int i = fromIndex; i < toIndex; i++)
			diagonals[type].add(value);
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

	public int diagonal(int type, int i, int j) {
		int d = i + j - 5;
		if (d < positions.size()) {
			int k = j - getValue(type, i, j, 0) + 2;
			if (k >= 2 && k < getValue(type, i, j, 1))
				return getValue(type, i, j, k);
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
