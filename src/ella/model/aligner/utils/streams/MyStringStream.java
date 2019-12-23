package ella.model.aligner.utils.streams;

import ella.model.io.MyParameters;

public class MyStringStream {

	private final static int CAPACITIY = MyParameters.INITIAL_STREAM_CAPACITIY;
	private String[] stream;
	private int pointer;

	public MyStringStream(int size) {
		stream = new String[size];
		pointer = 0;
	}

	public MyStringStream() {
		stream = new String[CAPACITIY];
		pointer = 0;
	}

	public void add(String s) {

		// doubling the size of the stream
		if (pointer == stream.length) {
			String[] streamCopy = new String[expandCapacity()];
			System.arraycopy(stream, 0, streamCopy, 0, stream.length);
			stream = streamCopy;
		}
		stream[pointer++] = s;
	}

	private int expandCapacity() {
		int newCapacity = stream.length * 2 + 2;
		if (newCapacity < 0)
			newCapacity = Integer.MAX_VALUE;
		return newCapacity;
	}

	public String[] toArray() {
		String[] out = new String[pointer];
		System.arraycopy(stream, 0, out, 0, pointer);
		return out;
	}

}
