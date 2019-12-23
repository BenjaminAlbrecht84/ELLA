package ella.model.aligner.utils.streams;

import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.io.MyParameters;

public class MyIntegerStream {

	private final static int CAPACITY = MyParameters.INITIAL_STREAM_CAPACITIY;
	private UnsignedIntArray stream;
	private long pointer;

	public MyIntegerStream(int size) {
		stream = new UnsignedIntArray(size);
		pointer = 0;
	}

	public MyIntegerStream() {
		stream = new UnsignedIntArray(CAPACITY);
		pointer = 0;
	}

	public void add(MyIntegerStream s) {
		if (pointer + s.size() > stream.size()) {
			long diff = pointer + s.size() - stream.size() + 1;
			UnsignedIntArray streamCopy = new UnsignedIntArray(stream, stream.size() + diff + 2);
			stream = streamCopy;
		}
		for (long i = 0; i < s.size(); i++)
			add(s.getStream().get(i));
	}

	public void add(int[] a) {
		if (pointer + a.length > stream.size()) {
			UnsignedIntArray streamCopy = new UnsignedIntArray(stream, stream.size() + a.length + 2);
			stream = streamCopy;
		}
		for (int i : a)
			add(i);
	}

	public void add(long i) {

		// doubling the size of the stream
		if (pointer == stream.size()) {
			UnsignedIntArray streamCopy = new UnsignedIntArray(stream, expandCapacity());
			stream = streamCopy;
		}
		stream.set(pointer++, i);

	}

	private long expandCapacity() {
		long newCapacity = stream.size() * 2 + 2;
		if (newCapacity < 0)
			newCapacity = 2 * Integer.MAX_VALUE - 8;
		return newCapacity;
	}

	public UnsignedIntArray toArray() {
		return new UnsignedIntArray(stream, pointer);
	}

	public void reset() {
		stream = new UnsignedIntArray(CAPACITY);
		pointer = 0;
	}

	public long size() {
		return pointer;
	}

	public UnsignedIntArray getStream() {
		return stream;
	}

	public void freeMemory() {
		stream = null;
	}

	public void addOffset(long j, long offset) {
		stream.set(j, stream.get(j) + offset);
	}

	public void ensureCapacity(long capacity) {
		if (stream.size() < capacity) {
			UnsignedIntArray streamCopy = new UnsignedIntArray(stream, capacity);
			// System.arraycopy(stream, 0, streamCopy, 0, stream.size());
			stream = streamCopy;
		}
	}

}
