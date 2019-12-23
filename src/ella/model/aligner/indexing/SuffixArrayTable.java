package ella.model.aligner.indexing;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.Minimizer;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.streams.MyStringStream;
import ella.model.aligner.utils.suffixArray.EnhancedSuffixArray;
import ella.model.aligner.utils.suffixArray.ImprovedBinarySearch;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class SuffixArrayTable {

	private int separationDepth, p, q;
	private IndexText text;
	private ConcurrentHashMap<ByteArrayWrapper, ConcurrentHashMap<CyclicSeedShape, EnhancedSuffixArray>> suffixArrayTable;

	public SuffixArrayTable(IndexText text, int separationDepth, int p, int q) {
		this.text = text;
		this.separationDepth = separationDepth;
		this.p = p;
		this.q = q;
		this.suffixArrayTable = new ConcurrentHashMap<>();
	}

	public ArrayList<Object[]> findPattern(String pattern, int minLen, int m, boolean useBucketTable, QueryContainer qC, Aligner aligner,
			ImprovedBinarySearch search) {

		ArrayList<Object[]> allMatches = new ArrayList<Object[]>();
		if (pattern.length() < q)
			return allMatches;

		// splitting pattern on left-most minimizer
		String[] splitPattern = splitPattern(pattern, aligner);

		// searching split-pattern
		ByteArrayWrapper prefix = getPrefix(splitPattern[1]);
		if (!suffixArrayTable.containsKey(prefix) || splitPattern[1].length() <= separationDepth)
			return allMatches;

		for (EnhancedSuffixArray sa : suffixArrayTable.get(prefix).values()) {
			ArrayList<Object[]> matchResult = sa.findPattern(text, splitPattern, minLen, m, useBucketTable, qC, aligner, search);
			if (matchResult != null) {
				allMatches.ensureCapacity(allMatches.size() + matchResult.size());
				allMatches.addAll(matchResult);
			}
		}

		return allMatches;

	}

	private String[] splitPattern(String pattern, Aligner aligner) {

		if (p == q) {
			String[] splitPattern = { "", pattern };
			return splitPattern;
		}

		Minimizer minimizing = new Minimizer(p, q);
		for (int i = 0; i < q - p + 1; i++)
			minimizing.addPMer(i, pattern.substring(i, i + p));

		int minPos = (int) minimizing.getMinimium();
		String[] splitPattern = { pattern.substring(0, minPos), pattern.substring(minPos) };
		return splitPattern;

	}

	public String[] getPrefixes() {
		MyStringStream stream = new MyStringStream(suffixArrayTable.keySet().size());
		for (ByteArrayWrapper key : suffixArrayTable.keySet())
			stream.add(key.getAAPrefix());
		return stream.toArray();
	}

	private ByteArrayWrapper getPrefix(String sequence) {
		byte[] kmer = new byte[separationDepth];
		for (int i = 0; i < Math.min(separationDepth, sequence.length()); i++) {
			byte b = 0;
			b |= Alphabet.getIndex(Alphabet.reduceCharacter(sequence.charAt(i)));
			kmer[i] = b;
		}
		return new ByteArrayWrapper(kmer);
	}

	public void addSuffixArray(EnhancedSuffixArray esa) {
		suffixArrayTable.putIfAbsent(esa.getCommonPrefix(), new ConcurrentHashMap<CyclicSeedShape, EnhancedSuffixArray>());
		suffixArrayTable.get(esa.getCommonPrefix()).put(esa.getSeedShape(), esa);
	}

	public long getByteSize() {
		long countSA = 0, countLCP = 0, countBT = 0;
		for (ByteArrayWrapper key : suffixArrayTable.keySet()) {
			for (EnhancedSuffixArray esa : suffixArrayTable.get(key).values()) {
				countSA += esa.getSAByteSize();
				countLCP += esa.getLCPByteSize();
			}
		}
		return countSA + countLCP + countBT;
	}

	public int getSeparationDepth() {
		return separationDepth;
	}

	public int getP() {
		return p;
	}

	public int getQ() {
		return q;
	}

	public ConcurrentHashMap<ByteArrayWrapper, ConcurrentHashMap<CyclicSeedShape, EnhancedSuffixArray>> getSuffixArrayTable() {
		return suffixArrayTable;
	}

	public int getNumberOfPrefixes() {
		return suffixArrayTable.keySet().size();
	}

	public void clearTable() {
		for (ByteArrayWrapper key : suffixArrayTable.keySet()) {
			for (EnhancedSuffixArray esa : suffixArrayTable.get(key).values())
				esa.freeMemory();
		}
		suffixArrayTable.clear();
	}

}
