package ella.model.aligner.utils.bigArrays;

import java.util.Arrays;

public class UnsignedByteArray {

	// TODO: apply http://www.omsn.de/blog/big-arrays-in-java

	private final static int MAX_SIZE = Integer.MAX_VALUE - 8;
	private final byte[] a1, a2;
	private int size;

	public UnsignedByteArray(long size) {
		this.size = (int) size;
		int[] split = splitLong(size);
		a1 = new byte[split[0]];
		a2 = new byte[split[1]];
	}

	public UnsignedByteArray(byte[] a1, byte[] a2) {
		this.size = (int) (a1.length + a2.length);
		this.a1 = a1;
		this.a2 = a2;
	}

	public UnsignedByteArray(UnsignedByteArray a) {
		this.size = (int) a.size();
		byte[] a1Clone = (byte[]) a.getArrays()[0];
		a1 = Arrays.copyOf(a1Clone, a1Clone.length);
		byte[] a2Clone = (byte[]) a.getArrays()[1];
		a2 = Arrays.copyOf(a2Clone, a2Clone.length);
	}

	public UnsignedByteArray(UnsignedByteArray a, long size) {
		this.size = (int) size;
		int[] split = splitLong(size);
		if (size > a.size()) {
			byte[] a1Clone = (byte[]) a.getArrays()[0];
			a1 = new byte[split[0]];
			System.arraycopy(a1Clone, 0, a1, 0, a1Clone.length);
			byte[] a2Clone = (byte[]) a.getArrays()[1];
			a2 = new byte[split[1]];
			System.arraycopy(a2Clone, 0, a2, 0, a2Clone.length);
		} else {
			byte[] a1Clone = (byte[]) a.getArrays()[0];
			a1 = new byte[split[0]];
			System.arraycopy(a1Clone, 0, a1, 0, (int) (size < MAX_SIZE ? size : MAX_SIZE));
			byte[] a2Clone = (byte[]) a.getArrays()[1];
			a2 = new byte[split[1]];
			System.arraycopy(a2Clone, 0, a2, 0, (int) (size > MAX_SIZE ? size - MAX_SIZE : 0));
		}
	}

	public byte get(long index) {
		int[] split = splitLong(index);
		byte value = split[0] != MAX_SIZE ? a1[split[0]] : a2[split[1]];
		return value;
	}

	public void set(long index, byte value) {
		int[] split = splitLong(index);
		if (split[0] != MAX_SIZE)
			a1[split[0]] = value;
		else
			a2[split[1]] = value;
	}

	public long size() {
		return Integer.toUnsignedLong(size);
	}

	private int[] splitLong(long value) {
		int v1, v2;
		if (value < MAX_SIZE) {
			v1 = (int) value;
			v2 = 0;
		} else {
			v1 = MAX_SIZE;
			v2 = (int) (value - MAX_SIZE);
		}
		int[] result = { v1, v2 };
		return result;
	}

	public Object[] getArrays() {
		Object[] arrays = { a1, a2 };
		return arrays;
	}

}
