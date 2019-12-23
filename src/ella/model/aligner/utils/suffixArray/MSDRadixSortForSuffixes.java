package ella.model.aligner.utils.suffixArray;

import java.util.ArrayDeque;
import java.util.ArrayList;

import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class MSDRadixSortForSuffixes {

	private final static int DELIMITER_BUCKET_INDEX = Alphabet.getAaString().length();
	private final static int BUCKET_SIZE_TRESHOLD = 10;
	private long[] pointers = new long[Alphabet.getAaString().length() + 1];
	private long[] counters = new long[Alphabet.getAaString().length() + 1];

	public Object[] run(IndexText text, UnsignedIntArray suffixes, CyclicSeedShape seedShape, ByteArrayWrapper prefix) {

		// initializing depth and helper arrays
		long n = suffixes.size();
		int initialDepth = prefix.getData().length;
		UnsignedIntArray sa1 = suffixes;
		UnsignedIntArray sa2 = new UnsignedIntArray(suffixes);
		UnsignedIntArray oracle = new UnsignedIntArray(n);

		// initializing sorted suffix array and lcp array
		LCPTreeArray lcpTreeArray = new LCPTreeArray(n);
		lcpTreeArray.setValue(0, 0);

		BucketTable bucketTable = null; // new BucketTable(seedShape.trim(prefix.getData().length), n);

		// initialize stack implementation
		ArrayDeque<Object[]> deque = new ArrayDeque<Object[]>();
		Object[] firstBucketInfo = { 0L, n, initialDepth, 0L };
		deque.add(firstBucketInfo);

		// iteratively conducting bucket sort in combination with insertion sort
		long beg, end, prefixBucketIndex;
		int depth;
		while (!deque.isEmpty()) {

			// retrieving bucket info
			Object[] bucketInfo = deque.pop();
			beg = (long) bucketInfo[0];
			end = (long) bucketInfo[1];
			depth = (int) bucketInfo[2];
			prefixBucketIndex = (long) bucketInfo[3];

			// setting up help arrays
			int b = depth & 1;
			UnsignedIntArray a = b == 0 ? sa1 : sa2;
			UnsignedIntArray aCopy = b == 0 ? sa2 : sa1;

			// checking if bucket size is small enough for doing a naive insertion sort
			// continuing with bucket sort otherwise
			if (end - beg <= BUCKET_SIZE_TRESHOLD) {
				insertionSort(beg, end, depth, initialDepth, a, sa2, lcpTreeArray, seedShape, text, prefixBucketIndex, bucketTable);
				continue;
			}

			// initializing counters
			for (int i = 0; i < counters.length; i++) {
				counters[i] = 0;
			}

			// setting up oracle array, technical speed-up trick for limiting cache and TLB misses
			// adopted from the paper "Engineering Radix Sort for Strings" of Kärkkäinen & Rantala (2008)
			for (long i = beg; i < end; i++) {
				int pos = readTextPosition(text, a.get(i) + depth);
				if (pos != DELIMITER_BUCKET_INDEX && !seedShape.usePosition(depth)) {
					oracle.set(i - beg, 0);
				} else {
					oracle.set(i - beg, pos);
				}
			}

			// setting up counters
			for (long i = beg; i < end; i++) {
				counters[(int) oracle.get(i - beg)]++;
			}

			// setting up pointers and new bucket info
			pointers[0] = beg;
			for (int i = 1; i < counters.length; i++) {
				pointers[i] = pointers[i - 1] + counters[i - 1];
				if (pointers[i - 1] < pointers[i]) {
					if (i - 1 != DELIMITER_BUCKET_INDEX) {
						long bucketIndex = bucketTable != null
								? bucketTable.setBucketEntry(Alphabet.getCharacter(i - 1), depth - initialDepth, prefixBucketIndex, pointers[i - 1])
								: 0;
						Object[] newBucketInfo = { pointers[i - 1], pointers[i], depth + 1, bucketIndex };
						deque.push(newBucketInfo);
					}
					if (pointers[i - 1] > 0 && lcpTreeArray.getEntryAtPosition(pointers[i - 1]) == 0)
						lcpTreeArray.setValue(depth, pointers[i - 1]);
				}
			}

			// sorting buckets into helper array
			ArrayList<long[]> delimiterList = new ArrayList<long[]>();
			for (long i = beg; i < end; i++) {
				long j = pointers[(int) oracle.get(i - beg)]++;
				if (oracle.get(i - beg) == DELIMITER_BUCKET_INDEX) {
					long[] tuple = { a.get(i), j };
					delimiterList.add(tuple);
					if (j > 0 && lcpTreeArray.getEntryAtPosition(j) == 0)
						lcpTreeArray.setValue(depth, j);
				} else {
					aCopy.set(j, a.get(i));
					// aCopy[j] = a[i];
				}
			}

			// reporting those suffixes that have reached a delimiter symbol
			for (long[] tuple : delimiterList) {
				sa2.set(tuple[1], tuple[0]);
				// sa2[tuple[1]] = tuple[0];
			}

		}

		// setting up arrays
		lcpTreeArray.cmpTreeArray(n);
		UnsignedIntArray suffixArray = new UnsignedIntArray(sa2);
		// this.sanityCheck(suffixArray, lcpTreeArray, text, seedShape);

		// freeing memory
		sa1 = sa2 = oracle = null;

		Object[] result = { suffixArray, lcpTreeArray, bucketTable };
		return result;

	}

	private void insertionSort(long beg, long end, int depth, int initialDepth, UnsignedIntArray sa, UnsignedIntArray sortedSuffixArray,
			LCPTreeArray lcpTreeArray, CyclicSeedShape seedShape, IndexText text, long prefixBucketIndex, BucketTable bucketTable) {

		// sorting suffixes into list
		ArrayList<Object[]> list = new ArrayList<Object[]>((int) (end - beg));
		Object[] firstTuple = { sa.get(beg), 0 };
		list.add(firstTuple);
		for (long i = beg + 1; i < end; i++) {

			// looking for a suffix lexicographically larger than or equal to s1
			// additionally using conducted character comparisons to asses lcp values
			long s1 = sa.get(i);
			boolean added = false;
			int[] comp = null, lastComp = null;
			for (int j = 0; j < list.size(); j++) {
				long s2 = (long) list.get(j)[0];
				comp = compareSuffixes(s1, s2, depth, text, seedShape);
				if (comp[0] <= 0) {
					Object[] tuple = { s1, 0 };
					list.add(j, tuple);
					list.get(j)[1] = lastComp != null ? lastComp[1] : 0;
					list.get(j + 1)[1] = comp[1];
					added = true;
					break;
				}
				lastComp = comp;
			}

			// could not find any suffix lexicographically larger than s1, consequently s2 is added to the end
			if (!added) {
				Object[] tupel = { s1, 0 };
				list.add(tupel);
				list.get(list.size() - 1)[1] = lastComp != null ? lastComp[1] : 0;
			}

		}

		// inserting sorted list of suffixes into sorted suffix array
		// setting lcp values for sorted suffixes
		long j = beg;
		for (int i = 0; i < list.size(); i++) {

			// setting suffix array
			sortedSuffixArray.set(j, (long) list.get(i)[0]);
			// sortedSuffixArray[j] = list.get(i)[0];

			// setting lcp array
			if (i > 0)
				lcpTreeArray.setValue((int) list.get(i)[1], j);

			// setting bucket table
			if (bucketTable != null) {
				long bucketIndex = prefixBucketIndex;
				for (int d = depth; d < bucketTable.getDepth() + initialDepth; d++) {
					int p = readTextPosition(text, sortedSuffixArray.get(j) + d);
					if (p == DELIMITER_BUCKET_INDEX)
						break;
					bucketIndex = bucketTable.setBucketEntry(Alphabet.getCharacter(p), d - initialDepth, bucketIndex, j);
				}
			}

			j++;
		}

	}

	private int[] compareSuffixes(long s1, long s2, int depth, IndexText text, CyclicSeedShape seedShape) {
		int len = depth;
		int p1 = readTextPosition(text, s1 + len);
		int p2 = readTextPosition(text, s2 + len);
		while (p1 != DELIMITER_BUCKET_INDEX && p2 != DELIMITER_BUCKET_INDEX && (p1 == p2 || !seedShape.usePosition(len))) {
			len++;
			p1 = readTextPosition(text, s1 + len);
			p2 = readTextPosition(text, s2 + len);
		}
		int[] result = { Integer.compare(p1, p2), len };
		return result;
	}

	private int readTextPosition(IndexText text, long pos) {
		int p = pos < text.length() ? text.readPosition(pos) : 127;
		if (p == 127)
			return DELIMITER_BUCKET_INDEX;
		return Alphabet.reducePosition(p);
	}

	// ONLY FOR DEBUGGING ----------------------------------------------------

	private void sanityCheck(int[] sa, LCPTreeArray lcpTreeArray, IndexText text, CyclicSeedShape seedShape) {
		for (int i = 1; i < sa.length; i++) {
			int s1 = sa[i - 1];
			int s2 = sa[i];
			int[] comp = compareSuffixes(s1, s2, 0, text, seedShape);
			if (comp[0] > 0) {
				System.out.println("ERROR: Wrong lexicographical ordering... " + (i - 1) + " " + i);
				System.exit(0);
			}
			if (comp[1] != lcpTreeArray.getValueAtPosition(i)) {
				System.out.println("ERROR: Wrong lcp value... " + i + " -> " + comp[1] + " != " + lcpTreeArray.getValueAtPosition(i));
				System.exit(0);
			}
		}
	}

}
