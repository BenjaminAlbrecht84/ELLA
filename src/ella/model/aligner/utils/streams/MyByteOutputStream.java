package ella.model.aligner.utils.streams;

import java.io.File;
import java.io.FileOutputStream;

import ella.model.aligner.utils.bigArrays.UnsignedByteArray;

public class MyByteOutputStream implements ByteStream {

	private UnsignedByteArray stream;
	private int size, pointer;
	private long writtenBytes = 0;
	private File out;
	private boolean append;

	public MyByteOutputStream(File out, int size, boolean append) {
		this.size = size;
		this.out = out;
		this.append = append;
		stream = new UnsignedByteArray(size);
		pointer = 0;
	}

	public void add(byte b) {
		// flushing content to output file if necessary
		if (pointer == size)
			writeToFile();
		stream.set(pointer++, b);
		writtenBytes++;
	}

	public void writeToFile() {
		FileOutputStream fos;
		try {
			fos = new FileOutputStream(out.getAbsolutePath(), append);
			UnsignedByteArray a = toArray();
			fos.write((byte[]) a.getArrays()[0]);
			fos.write((byte[]) a.getArrays()[1]);
			fos.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		stream = new UnsignedByteArray(size);
		pointer = 0;
	}

	public UnsignedByteArray toArray() {
		return new UnsignedByteArray(stream, pointer);
	}

	public long getWrittenBytes() {
		return writtenBytes;
	}

	public int size() {
		return pointer;
	}

}
