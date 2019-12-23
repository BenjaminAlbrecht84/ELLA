package ella.model.aligner.aligning.seeding.chaining;

import ella.model.aligner.aligning.seeding.Seed;

import java.util.Arrays;

public class SeedRadixSort {

	public Seed[] sort(Seed[] seeds, int referenceLength) {

		Seed[] a1 = seeds;
		Seed[] a2 = new Seed[a1.length];
		for (int k = 0; k < 4; k++) {
			Seed[] i1 = k % 2 == 0 ? a1 : a2;
			Seed[] i2 = k % 2 == 1 ? a1 : a2;
			countingSort(i1, i2, k, referenceLength);
		}
		return a1;

	}

	private final int[] pointers = new int[256];

	private Seed[] countingSort(Seed[] input, Seed[] output, int k, int referenceLength) {

		Arrays.fill(pointers, 0);
//		int[] pointers = new int[256];
		for (int i = 0; i < input.length; i++) {
			int v = input[i].getYCoordinate(referenceLength);
			int b = (v >> (k * 8)) & 255;
			pointers[b]++;
		}
		for (int i = 1; i < pointers.length; i++)
			pointers[i] += pointers[i - 1];

		for (int i = input.length - 1; i >= 0; i--) {
			int v = input[i].getYCoordinate(referenceLength);
			int b = (v >> (k * 8)) & 255;
			output[pointers[b] - 1] = input[i];
			pointers[b]--;
		}

		return output;

	}

}
