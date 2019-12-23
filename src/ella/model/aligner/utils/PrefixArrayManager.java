package ella.model.aligner.utils;

import java.util.HashMap;
import java.util.Set;

import ella.model.aligner.utils.streams.MyIntegerStream;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class PrefixArrayManager {

	private HashMap<ByteArrayWrapper, Integer> prefixToIndex = new HashMap<ByteArrayWrapper, Integer>();
	private MyIntegerStream[] streams;

	public PrefixArrayManager(int sepDepth) {
		cmpPrefixesRec(new byte[sepDepth], 0, Alphabet.getReducedAminoacids());
		streams = new MyIntegerStream[prefixToIndex.keySet().size()];
		for (ByteArrayWrapper prefix : prefixToIndex.keySet())
			streams[prefixToIndex.get(prefix)] = new MyIntegerStream();
	}

	private void cmpPrefixesRec(byte[] pref, int i, String alphabet) {
		if (i == pref.length) {
			prefixToIndex.put(new ByteArrayWrapper(pref), prefixToIndex.keySet().size());
		} else {
			for (int j = 0; j < alphabet.length(); j++) {
				byte[] prefCopy = pref.clone();
				prefCopy[i] = (byte) Alphabet.getIndex(alphabet.charAt(j));
				cmpPrefixesRec(prefCopy, i + 1, alphabet);
			}
		}
	}

	public void addValue(ByteArrayWrapper prefix, int[] values, int size) {
		streams[prefixToIndex.get(prefix)].ensureCapacity(size);
		streams[prefixToIndex.get(prefix)].add(values);
	}

	public void addValue(ByteArrayWrapper prefix, MyIntegerStream stream, Integer size) {
		streams[prefixToIndex.get(prefix)].ensureCapacity(size);
		streams[prefixToIndex.get(prefix)].add(stream);
	}

	public void addValue(ByteArrayWrapper prefix, long value) {
		streams[prefixToIndex.get(prefix)].add(value);
	}

	public void addValue(ByteArrayWrapper prefix, int[] values) {
		streams[prefixToIndex.get(prefix)].add(values);
	}

	public void clear() {
		prefixToIndex = null;
		streams = null;
	}

	public Set<ByteArrayWrapper> getPrefixes() {
		return prefixToIndex.keySet();
	}

	public MyIntegerStream getStream(ByteArrayWrapper prefix) {
		return streams[prefixToIndex.get(prefix)];
	}

	public void freeMemory(ByteArrayWrapper prefix) {
		streams[prefixToIndex.get(prefix)] = null;
	}

	public MyIntegerStream get(ByteArrayWrapper prefix) {
		return streams[prefixToIndex.get(prefix)];
	}

}
