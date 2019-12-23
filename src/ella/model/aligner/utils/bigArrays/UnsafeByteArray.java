package ella.model.aligner.utils.bigArrays;

import java.lang.reflect.Field;

public class UnsafeByteArray {

	private long address;
	private final static int BYTE_SIZE = 1;
	private long size;

	public UnsafeByteArray(UnsignedByteArray a) {
		address = getUnsafe().allocateMemory(a.size() * BYTE_SIZE);
		for (long i = 0; i < a.size(); i++)
			set(i, a.get(i));
		size = a.size();
	}

	public long size() {
		return size;
	}

	public void set(long idx, int val) {
		getUnsafe().putInt(address + idx * BYTE_SIZE, val);
	}

	public byte get(long idx) {
		return getUnsafe().getByte(address + idx * BYTE_SIZE);
	}

	public void free() {
		getUnsafe().freeMemory(address);
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
		UnsignedByteArray a = new UnsignedByteArray(size);
		for (long i = 0; i < size; i++)
			a.set(i, get(i));
		return a.getArrays();
	}

}
