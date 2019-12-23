package ella.model.io;

import java.nio.ByteBuffer;

import ella.model.aligner.utils.streams.ByteStream;

public class ByteWriter {

	public static void write(ByteStream byteStream, byte b) {
		byteStream.add(b);
	}

	public static void write(ByteStream byteStream, int[] a) {
		for (int i : a)
			write(byteStream, i);
	}

	public static void write(ByteStream byteStream, int i) {
		byte[] bytes = ByteBuffer.allocate(4).putInt(i).array();
		write(byteStream, bytes);
	}

	public static void write(ByteStream byteStream, long l) {
		byte[] bytes = ByteBuffer.allocate(8).putLong(l).array();
		write(byteStream, bytes);
	}

	public static void write(ByteStream byteStream, byte[] bytes) {
		for (byte b : bytes)
			write(byteStream, b);
	}

	public static void write(ByteStream byteStream, String s) {
		for (int i = 0; i < s.length(); i++)
			byteStream.add((byte) s.charAt(i));
	}

}
