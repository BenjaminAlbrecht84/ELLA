package ella.model.aligner.utils.suffixArray;

import java.util.Comparator;

import ella.model.aligner.utils.bigArrays.UnsignedByteArray;

public class LCPTreeArray {

	private UnsignedByteArray lcpTreeArray;

	public LCPTreeArray(UnsignedByteArray lcpTreeArray, String prefix) {
		this.lcpTreeArray = lcpTreeArray;
	}

	public LCPTreeArray(long n) {
		lcpTreeArray = new UnsignedByteArray(2 * n);
	}

	public void setValue(int value, long pos) {
		storeValue(pos, value);
	}

	public void cmpTreeArray(long n) {
		computeLCPTreeRec(0, n - 1, n);
	}

	public int[] byteToIntArray(byte[] byteArray) {
		int intArray[] = new int[byteArray.length / 4];
		int offset = 0;
		for (int i = 0; i < intArray.length; i++) {
			intArray[i] = (byteArray[3 + offset] & 0xFF) | ((byteArray[2 + offset] & 0xFF) << 8) | ((byteArray[1 + offset] & 0xFF) << 16)
					| ((byteArray[0 + offset] & 0xFF) << 24);
			offset += 4;
		}
		return intArray;
	}

	public class BigValueComparator implements Comparator<int[]> {
		@Override
		public int compare(int[] v1, int[] v2) {
			int p1 = v1[0];
			int p2 = v2[0];
			return Integer.compare(p1, p2);
		}
	}

	private int computeLCPTreeRec(long i, long j, long n) {

		if (i == j)
			return getValueAtPosition(i);
		if (i + 1 == j)
			return getValueAtPosition(j);

		int mid = divide(i + j, 2);
		int v1 = computeLCPTreeRec(i, mid, n);
		int v2 = computeLCPTreeRec(mid, j, n);

		int v = Math.min(v1, v2);

		long pos = n + divide((i + j), 2);
		storeValue(pos, v);
		return v;
	}

	private int divide(long num, long denum) {
//		double v1 = Double.valueOf(num);
//		double v2 = Double.valueOf(denum);
		double v1 = num;
		double v2 = denum;
		return (int) Math.floor(v1 / v2);
	}

	private void storeValue(long pos, int value) {
		value = value < 128 ? value : 127;
		byte b = 0;
		b |= value;
		lcpTreeArray.set(pos, b);
	}

	private int getValue(long pos) {
		return lcpTreeArray.get(pos);
	}

	private int querySmallValues(long pos) {
		return (int) lcpTreeArray.get(pos);
	}

	public int getTreeValue(long L, long R) {
		long pos;
		if (R - L < 2) {
			pos = R;
		} else {
			pos = lcpTreeArray.size() / 2 + divide(L + R, 2);
		}
		return getValue(pos);
	}

	public int getValueAtPosition(long pos) {
		return getValue(pos);
	}

	public int getUpperBoundValueAtPosition(long pos) {
		int value = querySmallValues(pos);
		value = value >= 0 ? value : 127;
		return value;
	}

	public int getEntryAtPosition(long pos) {
		return lcpTreeArray.get(pos);
	}

	public long getIntSize() {
		return lcpTreeArray.size() / 4;
	}

	public UnsignedByteArray getLcpTreeArray() {
		return lcpTreeArray;
	}

	public void freeMemory() {
		lcpTreeArray = null;
	}

}
