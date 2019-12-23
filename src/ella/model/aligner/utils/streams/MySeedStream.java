package ella.model.aligner.utils.streams;

import ella.model.aligner.aligning.seeding.Seed;
import ella.model.io.MyParameters;

public class MySeedStream {

	private final static int CAPACITY = MyParameters.INITIAL_STREAM_CAPACITIY;
	private Seed[] stream;
	private int pointer;

	public MySeedStream(int size) {
		stream = new Seed[size];
		pointer = 0;
	}

	public MySeedStream() {
		stream = new Seed[CAPACITY];
		pointer = 0;
	}

	public void add(Seed s) {

		// doubling the size of the stream
		if (pointer == stream.length) {
			Seed[] streamCopy = new Seed[expandCapacity()];
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

	public Seed[] toArray() {
		Seed[] out = new Seed[pointer];
		System.arraycopy(stream, 0, out, 0, pointer);
		return (Seed[]) out;
	}

	public int size() {
		return pointer;
	}

	public void reset() {
		stream = new Seed[CAPACITY];
		pointer = 0;
	}

}
