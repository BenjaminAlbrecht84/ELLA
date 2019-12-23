package ella.model.aligner.utils.frameshiftAligner;

import java.util.ArrayList;

import ella.model.aligner.utils.bigArrays.UnsafeIntArray;
import ella.model.io.MyParameters;

public class DP_Matrix_unsafe {

	private int MIN_VALUE = Integer.MIN_VALUE / 2;

	private int pointer = 0, offset = 100;
	private ArrayList<UnsafeIntArray[]> diagonals = new ArrayList<UnsafeIntArray[]>();
	private int gop = MyParameters.SCORING_MATRIX.getGapOpen(), gep = MyParameters.SCORING_MATRIX.getGapExtend();

	public void init(int n, int m) {
		pointer = 0;
		diagonals.ensureCapacity(n * m);
	}

	public UnsafeIntArray[] nextDiagonal(int size) {
		if (diagonals.size() == pointer) {
			diagonals.add(new UnsafeIntArray[3]);
			for (int i = 0; i < 3; i++)
				diagonals.get(pointer)[i] = new UnsafeIntArray(size > 2 ? size + offset : 2);
		}
		ensureSize(size);
		return diagonals.get(pointer++);
	}

	private void ensureSize(int size) {
		if (size > diagonals.get(pointer)[0].size()) {
			for (int i = 0; i < 3; i++)
				diagonals.get(pointer)[i].free();
			for (int i = 0; i < 3; i++)
				diagonals.get(pointer)[i] = new UnsafeIntArray(size + offset);
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

	public UnsafeIntArray getDiagonals(int type, int i, int j) {
		int d = i + j - 5;
		if (d >= 0 && d < diagonals.size())
			return diagonals.get(d)[type];
		return null;
	}

	public int diagonal(int type, int i, int j) {
		int d = i + j - 5;
		if (d < diagonals.size()) {
			UnsafeIntArray diagonal = diagonals.get(d)[type];
			int k = j - diagonal.get(0) + 2;
			if (k >= 2 && k < diagonal.get(1))
				return diagonal.get(k);
		}
		return MIN_VALUE;
	}

	public void freeMemory() {
		for (int i = 0; i < diagonals.size(); i++) {
			for (int j = 0; j < 3; j++)
				diagonals.get(i)[j].free();
		}
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
