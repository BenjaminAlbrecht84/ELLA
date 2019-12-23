package ella.model.aligner.utils.suffixArray;

import java.util.HashMap;

import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CyclicSeedShape;

public class SearchTable {

	public static int depth = 3;

	private int[][] buckets;
	private int[] indexArray;
	private int sigmaSize, prefixLength;
	private CyclicSeedShape seedShape;

	private HashMap<Integer, Integer> posTranslator;

	public SearchTable(CyclicSeedShape seedShape, String prefix) {
		this.seedShape = seedShape;
		setUpLengthTranslator(seedShape);
		initAlphabet();
		int n = 0;
		for (int i = 0; i < depth; i++)
			n += (int) Math.pow(sigmaSize + 1, i) * sigmaSize;
		buckets = new int[3][n];
		prefixLength = prefix.length();
	}

	private void setUpLengthTranslator(CyclicSeedShape seedShape2) {
		posTranslator = new HashMap<Integer, Integer>();
		int counter = -1, i = 0;
		while (counter < depth) {
			if (seedShape.usePosition(i))
				counter++;
			posTranslator.put(i++, counter);
		}
	}

	private void initAlphabet() {
		String sigma = Alphabet.getReducedAminoacids();
		indexArray = new int[256];
		for (int i = 0; i < indexArray.length; i++)
			indexArray[i] = -1;
		sigmaSize = sigma.length();
		for (int i = 0; i < sigma.length(); i++)
			indexArray[sigma.charAt(i)] = i;
	}

	public int[] findBucketEntry(String s) {
		int index = 0;
		for (int i = prefixLength; i < Math.min(s.length(), depth + prefixLength); i++) {
			index = cmpBucketIndex(s.charAt(i), i, index);
			if (buckets[2][index] != 0) {
				int[] result = { buckets[0][index], buckets[1][index], buckets[2][index] };
				return result;
			}
		}
		return null;
	}

	public void setBucketEntry(String s, long[] range, int length) {
		int l = length - prefixLength;
		l = posTranslator.containsKey(l) ? posTranslator.get(l) : depth + 1;
		if (l <= depth) {
			int index = cmpBucketIndex(s, prefixLength, length);
			int[] result = { (int) range[0], (int) range[1], length };
			for (int i = 0; i < 3; i++)
				buckets[i][index] = result[i];
		}
	}

	public int cmpBucketIndex(String s, int l, int r) {
		int index = 0;
		for (int i = l; i < Math.min(s.length(), r); i++) {
			if (!seedShape.usePosition(i))
				continue;
			int aaIndex = indexArray[Alphabet.reduceCharacter(s.charAt(i))] + 1;
			int pos = posTranslator.get(i - prefixLength);
			index += Math.pow(sigmaSize + 1, depth - 1 - pos) * aaIndex;
		}
		return index;
	}

	public int cmpBucketIndex(char c, int pos, int offset) {
		if (!seedShape.usePosition(pos))
			return offset;
		pos = posTranslator.get(pos - prefixLength);
		int aaIndex = indexArray[Alphabet.reduceCharacter(c)] + 1;
		int index = offset + (int) Math.pow(sigmaSize + 1, depth - 1 - pos) * aaIndex;
		return index;
	}

}
