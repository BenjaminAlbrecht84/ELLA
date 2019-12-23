package ella.model.aligner.utils.suffixArray;

import java.util.ArrayList;
import java.util.HashMap;

import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;

public class BucketTable {

	public static int depth = 3;

	private UnsignedIntArray buckets;
	private int[] indexArray;
	private int sigmaSize;
	private long n;
	private CyclicSeedShape seedShape;

	private HashMap<Integer, Integer> posTranslator;

	public BucketTable(long n, int depth, UnsignedIntArray buckets, CyclicSeedShape seedShape) {
		this.n = n;
		this.depth = depth;
		this.buckets = buckets;
		this.seedShape = seedShape;
		setUpLengthTranslator(seedShape);
		initAlphabet();
	}

	public BucketTable(CyclicSeedShape seedShape, long n) {
		this.n = n;
		this.seedShape = seedShape;
		setUpLengthTranslator(seedShape);
		initAlphabet();
		int bucketSize = 0;
		for (int i = 0; i < depth; i++)
			bucketSize += (int) Math.pow(sigmaSize + 1, i) * sigmaSize;
		buckets = new UnsignedIntArray(bucketSize + 1);
	}

	private void setUpLengthTranslator(CyclicSeedShape seedShape2) {
		posTranslator = new HashMap<Integer, Integer>();
		int counter = 0;
		for (int i = 0; i < depth; i++) {
			if (seedShape.usePosition(i)) {
				posTranslator.put(i, counter);
				counter++;
			}
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

	public long setBucketEntry(String s, int entry) {
		long index = cmpBucketIndex(s);
		if (buckets.get(index) == 0 && index > -1 && index < buckets.size())
			buckets.set(index, entry + 1);
		return index;
	}

	public long setBucketEntry(char c, int len, long offset, long entry) {
		long index = cmpBucketIndex(c, len, offset);
		if (index > -1 && buckets.get(index) == 0 && index < buckets.size())
			buckets.set(index, entry + 1);
		return index;
	}

	public Object[] getRightMostPrefixRange(String s, int m) {
		long index = 0;
		String prefix = "";
		Object[] resBucket = null;
		for (int i = 0; i < Math.min(s.length(), depth); i++) {
			prefix += s.charAt(i);
			index = cmpBucketIndex(s.charAt(i), i, index);
			long start = buckets.get(index);
			if (start != 0) {
				long end = findNextBucket(index, prefix);
				Object[] b = { start - 1, end - 1, i + 1, end - start };
				resBucket = b;
				if (end - start <= m)
					return resBucket;
			}
		}
		return resBucket;
	}

	public long[] getBucketRange(String s) {
		long index = cmpBucketIndex(s);
		if (index > -1 && index < buckets.size()) {
			long start = buckets.get(index);
			if (start != 0) {
				long end = findNextBucket(index, s);
				long[] b = { start - 1, end - 1 };
				return b;
			}
		}
		return null;
	}

	private long findNextBucket(long index, String s) {

		// happens only if bucket is the last bucket
		if (s.isEmpty())
			return n + 1;

		// position suppressed by seed shape
		int pos = s.length() - 1;
		if (!seedShape.usePosition(pos))
			return findNextBucket(index, s.substring(0, s.length() - 1));

		// checking buckets of lexicographically higher sequences that are of equal size
		int aaIndex = indexArray[Alphabet.reduceCharacter(s.charAt(s.length() - 1))];
		pos = posTranslator.get(pos);
		for (int i = aaIndex + 1; i < sigmaSize; i++) {
			index += Math.pow(sigmaSize + 1, depth - 1 - pos);
			if (buckets.get(index) != 0)
				return buckets.get(index);
		}
		index -= Math.pow(sigmaSize + 1, depth - 1 - pos) * (sigmaSize);

		// checking buckets of lexicographically higher sequences that are of one shorter length
		return findNextBucket(index, s.substring(0, s.length() - 1));

	}

	public long cmpBucketIndex(String s) {
		if (s.length() > depth)
			return -1;
		long index = 0;
		for (int i = 0; i < Math.min(s.length(), depth); i++) {
			if (!seedShape.usePosition(i))
				continue;
			int aaIndex = indexArray[Alphabet.reduceCharacter(s.charAt(i))] + 1;
			int pos = posTranslator.get(i);
			index += Math.pow(sigmaSize + 1, depth - 1 - pos) * aaIndex;
		}
		return index;
	}

	public long cmpBucketIndex(char c, int pos, long offset) {
		if (pos >= depth)
			return -1;
		if (!seedShape.usePosition(pos))
			return offset;
		pos = posTranslator.get(pos);
		int aaIndex = indexArray[Alphabet.reduceCharacter(c)] + 1;
		long index = offset + (long) Math.pow(sigmaSize + 1, depth - 1 - pos) * aaIndex;
		return index;
	}

	public int getDepth() {
		return depth;
	}

	public long getN() {
		return n;
	}

	public CyclicSeedShape getSeedShape() {
		return seedShape;
	}

	public void freeMemory() {
		buckets = null;
	}

	public long getIntSize() {
		return buckets.size();
	}

	public UnsignedIntArray toArray() {
		return buckets;
	}

	// ONLY FOR DEBUGGING -------------------------------------

	public void printTable() {
		ArrayList<String> allKeys = new ArrayList<String>();
		cmpAllKeysRec(allKeys, "");
		for (String key : allKeys) {
			long index = cmpBucketIndex(key);
			long[] bucket = getBucketRange(key);
			if (bucket != null)
				System.out.println(key + "\t" + index + "\t" + bucket[0] + " " + bucket[1]);
		}
	}

	private void cmpAllKeysRec(ArrayList<String> allKeys, String key) {
		if (key.length() > depth)
			return;
		if (!key.isEmpty())
			allKeys.add(key);
		String sigma = Alphabet.getReducedAminoacids();
		for (int i = 0; i < sigma.length(); i++) {
			char c = sigma.charAt(i);
			if (c >= 65 && c <= 90) {
				String nextKey = key.concat(c + "");
				cmpAllKeysRec(allKeys, nextKey);
			}
		}
	}

}
