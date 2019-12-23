package ella.model.aligner.utils.frameshiftAligner;

import java.util.ArrayList;
import java.util.Arrays;

import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.bigArrays.UnsafeIntArray;
import ella.model.aligner.utils.frameshiftAligner.Zheng_XDrop_Frameshift.AliMode;
import ella.model.io.MyParameters;

public class Aligner_Mod6 {

	private final static int MIN_VALUE = Integer.MIN_VALUE / 2;
	private String protein, dna;

	private boolean verbose = false;

	private int gop, gep, F;

	public Aligner_Mod6() {
		this.gop = MyParameters.SCORING_MATRIX.getGapOpen();
		this.gep = MyParameters.SCORING_MATRIX.getGapExtend();
		this.F = MyParameters.FRAMESHIFT_PENALTY;
	}

	public Object[] run(String s1, String s2, AliMode mode, int xDiagonalDrop, DP_Matrix_unsafe matrix) {

		this.dna = s1.replace('T', 'U');
		this.protein = s2;

		if (protein.length() == 0 || dna.length() <= 2)
			return cmpTrivialAlignment(mode);

		if (mode == AliMode.LEFT) {
			protein = new StringBuilder(protein).reverse().toString();
			dna = new StringBuilder(dna).reverse().toString();
		}

		int n = protein.length() + 1;
		int m = dna.length() + 4;

		// initializing matrices
		matrix.init(n, m);

		// filling Matrix
		int[] startingCell = fillingMatrix(n, m, matrix, xDiagonalDrop, mode);
		if (startingCell == null || (startingCell[0] == -1 && startingCell[1] == -1)) {
			Object[] emptyResult = { "", "", "", 0, dna.length(), protein.length() };
			return emptyResult;
		}

		// starting traceback
		Object[] tracebackResult = traceback(matrix, mode, startingCell, n, m);

		return tracebackResult;

	}

	private Object[] cmpTrivialAlignment(AliMode mode) {
		if (dna.length() < 3) {
			int score = protein.length() == 0 ? 0 : -gop - protein.length() * gep;
			Object[] res = { gapString(protein.length()), protein, frameString(protein.length(), 1), score, dna.length(), 0 };
			return res;
		} else {
			StringBuffer p = new StringBuffer();
			for (int i = 0; i < dna.length() - 2; i += 3)
				p.append(CodonTranslator.translateCodon(dna.charAt(i), dna.charAt(i + 1), dna.charAt(i + 2)));
			Object[] res = { p, gapString(p.length()), frameString(protein.length(), 1), -gop - p.length() * gep, 0, 0 };
			return res;
		}
	}

	private int[] fillingMatrix(int n, int m, DP_Matrix_unsafe matrix, int xDiagonalDrop, AliMode mode) {

		int[] T_Primes = { MIN_VALUE, MIN_VALUE, MIN_VALUE, MIN_VALUE };
		int[] tracebackCell = { -1, -1, 0 };
		int[] T = { 0, 0, 0, 0 };
		int[][] B = { { 4, 4 }, { 4, 5 }, { 4, 6 }, { 4, 7 } };
		int k = 4, L = 4, U = 4;
		int dropOutCounter = 0;
		int maxScore = MIN_VALUE;

		while (k < m + n - 4) {

			// printMatrix(X, protein, dna);

			k = k + 1;
			int L_Prime = -1, U_Prime = -2;
			UnsafeIntArray[] diagonals = matrix.nextDiagonal(U - L + 3);
			for (int i = 0; i < 3; i++) {
				diagonals[i].set(0, U >= L ? L : m + 1);
				diagonals[i].set(1, U >= L ? U - L + 3 : 2);
			}

			if (L > U) {
				if (dropOutCounter++ == 4)
					return tracebackCell;
			} else {

				if (L > 3 && k - U - 1 > 0) {
					UnsafeIntArray X_diagonal = matrix.getDiagonals(2, k - L - 1, L);
					int x0 = X_diagonal.get(0), x1 = X_diagonal.get(1);
					UnsafeIntArray Y_diagonal = matrix.getDiagonals(1, k - L - 1, L);
					int y0 = Y_diagonal.get(0);
					diagonals[1].fill(2, Math.min(Math.max(L, x0), U + 1) - L + 2, MIN_VALUE);
					int fromIndex = Math.max(L, x0), toIndex = Math.min(U + 1, x0 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						diagonals[1].set(j - L + 2, X_diagonal.get(j - x0 + 2) - gop);
					fromIndex = Math.max(L, x0);
					toIndex = Math.min(U + 1, x0 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						diagonals[1].set(j - L + 2, Math.max(diagonals[1].get(j - L + 2), Y_diagonal.get(j - y0 + 2) - gep));
					diagonals[1].fill(Math.max(Math.min(U + 1, x0 + x1 - 2), L) - L + 2, U + 1 - L + 2, MIN_VALUE);
				} else
					for (int j = L; j < U + 1; j++)
						diagonals[1].set(j - L + 2, Math.max(matrix.X(k - j - 1, j) - gop, matrix.Y(k - j - 1, j) - gep));

				if (L > 6 && k - U > 0) {
					UnsafeIntArray X_diagonal = matrix.getDiagonals(2, k - L, L - 3);
					int x0 = X_diagonal.get(0), x1 = X_diagonal.get(1);
					UnsafeIntArray D_diagonal = matrix.getDiagonals(0, k - L, L - 3);
					int d0 = D_diagonal.get(0);
					diagonals[0].fill(2, Math.min(Math.max(L, d0 + 3), U + 1) - L + 2, MIN_VALUE);
					int fromIndex = Math.max(L, x0 + 3), toIndex = Math.min(U + 1, x0 + 3 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						diagonals[0].set(j - L + 2, X_diagonal.get(j - x0 - 3 + 2) - gop);
					fromIndex = Math.max(L, x0 + 3);
					toIndex = Math.min(U + 1, x0 + 3 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						diagonals[0].set(j - L + 2, Math.max(diagonals[0].get(j - L + 2), D_diagonal.get(j - d0 - 3 + 2) - gep));
					diagonals[0].fill(Math.max(Math.min(U + 1, x0 + 3 + x1 - 2), L) - L + 2, U + 1 - L + 2, MIN_VALUE);
				} else
					for (int j = L; j < U + 1; j++)
						diagonals[0].set(j - L + 2, Math.max(matrix.X(k - j, j - 3) - gop, matrix.D(k - j, j - 3) - gep));

				for (int j = L; j < U + 1; j++)
					diagonals[2].set(j - L + 2, diagonals[1].get(j - L + 2));
				for (int j = L; j < U + 1; j++)
					diagonals[2].set(j - L + 2, Math.max(diagonals[2].get(j - L + 2), diagonals[0].get(j - L + 2)));

				int[] a = new int[U - L + 1];
				if (L > 6 && k - U - 1 > 0) {
					UnsafeIntArray X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 3);
					int x0 = X_diagonal.get(0), x1 = X_diagonal.get(1);
					Arrays.fill(a, 0, Math.min(Math.max(L, x0 + 3), U + 1) - L, MIN_VALUE);
					int fromIndex = Math.max(L, x0 + 3), toIndex = Math.min(U + 1, x0 + 3 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						a[j - L] = X_diagonal.get(j - x0 - 3 + 2);
					Arrays.fill(a, Math.max(Math.min(U + 1, x0 + 3 + x1 - 2), L) - L, U + 1 - L, MIN_VALUE);
				} else
					for (int j = L; j < U + 1; j++)
						a[j - L] = matrix.X(k - j - 1, j - 3);
				if (L > 5 && k - U - 1 > 0) {
					UnsafeIntArray X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 2);
					int x0 = X_diagonal.get(0), x1 = X_diagonal.get(1);
					Arrays.fill(a, 0, Math.min(Math.max(L, x0 + 2), U + 1) - L, MIN_VALUE);
					int fromIndex = Math.max(L, x0 + 2), toIndex = Math.min(U + 1, x0 + 2 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						a[j - L] = Math.max(a[j - L], X_diagonal.get(j - x0 - 2 + 2) - F);
					Arrays.fill(a, Math.max(Math.min(U + 1, x0 + 2 + x1 - 2), L) - L, U + 1 - L, MIN_VALUE);
				} else
					for (int j = Math.max(L, 6); j < U + 1; j++)
						a[j - L] = Math.max(a[j - L], matrix.X(k - j - 1, j - 2) - F);
				if (L > 7 && k - U - 1 > 0) {
					UnsafeIntArray X_diagonal = matrix.getDiagonals(2, k - L - 1, L - 4);
					int x0 = X_diagonal.get(0), x1 = X_diagonal.get(1);
					Arrays.fill(a, 0, Math.min(Math.max(L, x0 + 4), U + 1) - L, MIN_VALUE);
					int fromIndex = Math.max(L, x0 + 4), toIndex = Math.min(U + 1, x0 + 4 + x1 - 2);
					for (int j = fromIndex; j < toIndex; j++)
						a[j - L] = Math.max(a[j - L], X_diagonal.get(j - x0 - 4 + 2) - F);
					Arrays.fill(a, Math.max(Math.min(U + 1, x0 + 4 + x1 - 2), L) - L, U + 1 - L, MIN_VALUE);
				} else
					for (int j = Math.max(L, 8); j < U + 1; j++)
						a[j - L] = Math.max(a[j - L], matrix.X(k - j - 1, j - 4) - F);
				for (int j = L; j < U + 1; j++)
					diagonals[2].set(j - L + 2, Math.max(diagonals[2].get(j - L + 2), a[j - L] + matchScore(k - j, j, false, mode)));

				for (int j = L; j < U + 1; j++) {

					int xMax = diagonals[2].get(j - L + 2);
					maxScore = maxScore < xMax ? xMax : maxScore;

					if (xMax > maxScore - MyParameters.X_BESTSCORE_DROP) {

						if (xMax > T_Primes[k & 3]) {
							T_Primes[k & 3] = xMax;
							if (mode == AliMode.LEFT || mode == AliMode.RIGHT) {
								if (xMax > tracebackCell[2]) {
									int[] cell = { k - j, j, xMax };
									tracebackCell = cell;
								}
							} else {
								if (k - j > tracebackCell[0]
										|| (k - j == tracebackCell[0] && (xMax > tracebackCell[2] || j > tracebackCell[1] + 2))) {
									int[] cell = { k - j, j, xMax };
									tracebackCell = cell;
								}
							}
						}

						L_Prime = L_Prime == -1 ? j : L_Prime;
						U_Prime = j;

					} else {
						diagonals[0].set(j - L + 2, MIN_VALUE);
						diagonals[1].set(j - L + 2, MIN_VALUE);
						diagonals[2].set(j - L + 2, MIN_VALUE);
					}

				}

			}

			int[] B_Prime = new int[2];
			B_Prime[0] = L_Prime;
			B_Prime[1] = L_Prime >= 0 ? U_Prime + 4 : U_Prime;
			B[(k - 1) & 3] = B_Prime;
			L = Math.max(B[k & 3][0], k + 2 - n);
			U = Math.min(B[k & 3][1], m - 3);

			T[k & 3] = T_Primes[k & 3];
			T_Primes[k & 3] = MIN_VALUE;

		}

		// matrix.printMatrix(protein, dna, n, m);

		return tracebackCell;

	}

	private Object[] traceback(DP_Matrix_unsafe matrix, AliMode mode, int[] startingCell, int n, int m) {

		StringBuilder ali1 = new StringBuilder();
		StringBuilder ali2 = new StringBuilder();
		StringBuilder frames = new StringBuilder();
		StringBuilder[] alignment = { ali1, ali2, frames };

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
		Object[] res = { alignment[1].toString(), alignment[0].toString(), alignment[2].toString(), aliScore, qNotAligned, rNotAligned };

		return res;

	}

	private int[] tracebackRec(DP_Matrix_unsafe matrix, int i, int j, StringBuilder[] alignment, AliMode mode) {

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
					int[] prevCell = { i - 1, j - 3 };
					alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
				}

				// checking for hit in different frames
				if ((j - 2) >= 2 && matrix.X(i - 1, j - 2) == val - score + F) {
					int[] prevCell = { i - 1, j - 2 };
					alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
				}
				if ((j - 4) >= 4 && matrix.X(i - 1, j - 4) == val - score + F) {
					int[] prevCell = { i - 1, j - 4 };
					alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
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
						int[] prevCell = { u, j };
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
						int[] prevCell = { i, t };
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
		char b = mode != AliMode.LEFT ? CodonTranslator.translateCodon(dna.charAt(j - 4), dna.charAt(j - 3), dna.charAt(j - 2))
				: CodonTranslator.translateCodon(dna.charAt(j - 2), dna.charAt(j - 3), dna.charAt(j - 4));
		return MyParameters.SCORING_MATRIX.getScore(a, b);
	}

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
