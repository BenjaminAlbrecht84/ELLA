package ella.model.aligner.utils.bigArrays;

import java.lang.reflect.Field;

public class UnsafeIntArray {

	private long address;
	private final static int BYTE_SIZE = 4;
	private long size;
	private sun.misc.Unsafe unsafe;

	public UnsafeIntArray(int size) {
		this.size = size;
		this.unsafe = getUnsafe();
		address = unsafe.allocateMemory(size * BYTE_SIZE);
	}

	public long size() {
		return size;
	}

	public void set(long idx, int val) {
		unsafe.putInt(address + idx * BYTE_SIZE, val);
	}

	public void fill(long fromIDX, long toIDX, int val) {
		long startAddress = address + fromIDX * BYTE_SIZE;
		long endAddress = address + toIDX * BYTE_SIZE;
		for (long pos = startAddress; pos < endAddress; pos += BYTE_SIZE)
			unsafe.putInt(pos, val);
//		for (long idx = fromIDX; idx < toIDX; idx++)
//			set(idx, val);
	}

	public int get(long idx) {
		return unsafe.getInt(address + idx * BYTE_SIZE);
	}

	public void free() {
		unsafe.freeMemory(address);
	}

	private sun.misc.Unsafe getUnsafe() {
		try {
			Field f = sun.misc.Unsafe.class.getDeclaredField("theUnsafe");
			f.setAccessible(true);
			return (sun.misc.Unsafe) f.get(null);
		} catch (Exception e) {
			System.err.println("Error obtaining unsafe: " + e.getMessage());
			return null;
		}
	}

	public Object[] getArrays() {
		UnsignedIntArray a = new UnsignedIntArray(size);
		for (long i = 0; i < size; i++)
			a.set(i, get(i));
		return a.getArrays();
	}

}
