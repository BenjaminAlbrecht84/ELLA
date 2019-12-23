package ella.model.aligner.utils;

import java.util.ArrayList;
import java.util.HashMap;

public class DNAPacker {

	private static final HashMap<Character, Integer> nucToIndex;
	static {
		nucToIndex = new HashMap<Character, Integer>();
		nucToIndex.put('A', 0);
		nucToIndex.put('C', 1);
		nucToIndex.put('G', 2);
		nucToIndex.put('T', 3);
		nucToIndex.put('a', 0);
		nucToIndex.put('c', 1);
		nucToIndex.put('g', 2);
		nucToIndex.put('t', 3);
	}

	public static byte[] unpackSequence(byte[] packed, int length, int bits) {

		byte[] result = new byte[length];
		long x = 0;
		int n = 0, pos = 0;
		int mask = (1 << bits) - 1;

		for (int i = 0; i < packed.length; i++) {
			x |= (packed[i] & 0xFF) << n;
			n += 8;

			while (n >= bits && pos < length) {
				result[pos++] = (byte) (x & mask);
				n -= bits;
				x >>>= bits;
			}
		}

		return result;

	}

	public static byte[] packSequence(String dna) {

		ArrayList<Byte> packed = new ArrayList<Byte>();
		byte p = 0;
		for (int i = 0; i < dna.length(); i++) {
			char c = dna.charAt(i);
			byte b = 0;
			int dnaIndex = nucToIndex.containsKey(c) ? nucToIndex.get(c) : 0;
			b |= dnaIndex << (i * 2) % 8;
			p |= b;
			if (i == dna.length() - 1 || (((i + 1) * 2) % 8 == 0 && i != 0)) {
				packed.add(p);
				p = 0;
			}
		}

		byte[] a = new byte[packed.size()];
		for (int i = 0; i < packed.size(); i++)
			a[i] = packed.get(i);

		return a;

	}
}
