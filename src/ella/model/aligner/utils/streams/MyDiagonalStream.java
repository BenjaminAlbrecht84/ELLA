package ella.model.aligner.utils.streams;

import ella.model.io.MyParameters;

public class MyDiagonalStream {

	private final static int CAPACITY = MyParameters.INITIAL_STREAM_CAPACITIY;
	private int[] stream;
	private int pointer;

	public MyDiagonalStream(int size) {
		stream = new int[size];
		pointer = 0;
	}

	public MyDiagonalStream() {
		stream = new int[CAPACITY];
		pointer = 0;
	}

	public void add(int[] a) {
		if (pointer + a.length > stream.length) {
			int[] aCopy = new int[expandCapacity(pointer + a.length)];
			System.arraycopy(stream, 0, aCopy, 0, pointer);
			stream = aCopy;
		}
		System.arraycopy(a, 0, stream, pointer, a.length);
		pointer += a.length;
	}

	public void add(int value) {
		// doubling the size of the stream
		if (pointer == stream.length) {
			int[] streamCopy = new int[expandCapacity(pointer + 1)];
			System.arraycopy(stream, 0, streamCopy, 0, stream.length);
			stream = streamCopy;
		}
		stream[pointer++] = value;

	}

	private int expandCapacity(int min) {
		int newCapacity = min;
		if (newCapacity < 0)
			newCapacity = Integer.MAX_VALUE - 8;
		return newCapacity;
	}

	public int get(int index) {
		return stream[index];
	}

	public void reset() {
		pointer = 0;
	}

	public int size() {
		return pointer;
	}

	public void freeMemory() {
		stream = null;
	}

	public void ensureCapacity(int capacity) {
		if (stream.length < capacity) {
			int[] streamCopy = new int[capacity];
			System.arraycopy(stream, 0, streamCopy, 0, stream.length);
			stream = streamCopy;
		}
	}

	public void set(int index, int value) {
		stream[index] = value;
	}

}
