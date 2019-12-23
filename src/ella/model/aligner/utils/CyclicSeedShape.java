package ella.model.aligner.utils;

import java.util.BitSet;

public class CyclicSeedShape {

	private BitSet seedShape;
	private int length;

	public CyclicSeedShape(String s) {
		this.length = s.length();
		this.seedShape = new BitSet();
		for (int i = 0; i < length; i++) {
			int b = Integer.parseInt(s.charAt(i) + "");
			if (b != 1 && b != 0)
				throw new IllegalArgumentException("Only 1s and 0s in seed-shapes allowed: " + s);
			seedShape.set(i, b != 0);
		}
	}

	public boolean usePosition(int pos) {
		if (pos >= length)
			pos = pos % length;
		return seedShape.get(pos);
	}

	public String toString() {
		StringBuilder builder = new StringBuilder(length);
		for (int i = 0; i < length; i++) {
			char c = seedShape.get(i) ? '1' : '0';
			builder.append(c);
		}
		return builder.toString();
	}

	public CyclicSeedShape trim(int len) {
		return new CyclicSeedShape(toString().substring(len));
	}

	public int getLength() {
		return length;
	}

}
