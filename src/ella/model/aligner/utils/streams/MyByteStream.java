package ella.model.aligner.utils.streams;

import ella.model.aligner.utils.bigArrays.UnsignedByteArray;
import ella.model.io.MyParameters;

public class MyByteStream implements ByteStream {

	private final static int MAX_RESIZE = 1000000;
	private final static int CAPACITY = MyParameters.INITIAL_STREAM_CAPACITIY;
	private UnsignedByteArray stream;
	private long pointer;

	public MyByteStream(int size) {
		stream = new UnsignedByteArray(size);
		pointer = 0;
	}

	public MyByteStream() {
		stream = new UnsignedByteArray(CAPACITY);
		pointer = 0;
	}

	public void add(MyByteStream s) {
		if (pointer + s.size() > stream.size()) {
			long diff = pointer + s.size() - stream.size() + 1;
			UnsignedByteArray streamCopy = new UnsignedByteArray(stream, stream.size() + diff + 2);
			stream = streamCopy;
		}
		for (int i = 0; i < s.size(); i++)
			add(s.getStream().get(i));
	}

	public UnsignedByteArray getStream() {
		return stream;
	}

	public void add(byte b) {

		// doubling the size of the stream
		if (pointer == stream.size()) {
			long resize = Math.min(MAX_RESIZE, stream.size());
			UnsignedByteArray streamCopy = new UnsignedByteArray(stream, stream.size() + resize);
			stream = streamCopy;
		}
		stream.set(pointer++, b);

	}

	public UnsignedByteArray toArray() {
		return new UnsignedByteArray(stream, pointer);
	}

	public long size() {
		return pointer;
	}

	public void reset() {
		stream = new UnsignedByteArray(CAPACITY);
		pointer = 0;
	}

}
