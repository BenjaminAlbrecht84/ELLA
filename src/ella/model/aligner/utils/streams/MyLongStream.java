package ella.model.aligner.utils.streams;

import ella.model.io.MyParameters;

public class MyLongStream {

    private final static int CAPACITY = MyParameters.INITIAL_STREAM_CAPACITIY;
    private long[] stream;
    private int pointer;

    public MyLongStream(int size) {
        stream = new long[size];
        pointer = 0;
    }

    public MyLongStream() {
        stream = new long[CAPACITY];
        pointer = 0;
    }

    public void add(long l) {

        // doubling the size of the stream
        if (pointer == stream.length) {
            long[] streamCopy = new long[expandCapacity()];
            System.arraycopy(stream, 0, streamCopy, 0, stream.length);
            stream = streamCopy;
        }
        stream[pointer++] = l;

    }

    private int expandCapacity() {
        int newCapacity = stream.length * 2 + 2;
        if (newCapacity < 0)
            newCapacity = Integer.MAX_VALUE;
        return newCapacity;
    }

    public long[] toArray() {
        long[] out = new long[pointer];
        System.arraycopy(stream, 0, out, 0, pointer);
        return (long[]) out;
    }

    public Long get(int i) {
        if (i < pointer)
            return stream[i];
        return null;
    }

    public void reset() {
        stream = new long[CAPACITY];
        pointer = 0;
    }

    public void reset(int size) {
        stream = stream.length < size ? new long[size] : stream;
        pointer = 0;
    }

    public int size() {
        return pointer;
    }

    public long[] getStream() {
        return stream;
    }

    public void freeMemory() {
        stream = null;
    }

    public void addOffset(int j, int offset) {
        stream[j] += offset;
    }

}
