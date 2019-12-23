package ella.model.aligner.utils.wrapper;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

import ella.model.aligner.utils.Alphabet;

public class ByteArrayWrapper {
	
	private final byte[] data;

	public ByteArrayWrapper(byte[] data) {
		if (data == null) {
			throw new NullPointerException();
		}
		this.data = data;
	}

	@Override
	public boolean equals(Object other) {
		if (!(other instanceof ByteArrayWrapper)) {
			return false;
		}
		return Arrays.equals(data, ((ByteArrayWrapper) other).data);
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(data);
	}

	public byte[] getData() {
		return data;
	}
	
	public int getLength() {
		return data.length;
	}
	
	public String getAAPrefix() {
		StringBuilder builder = new StringBuilder(data.length);
		for (byte b : data) {
			byte[] bytes = { b };
			ByteBuffer buffer = ByteBuffer.wrap(bytes);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			builder.append(Alphabet.getCharacter((int) buffer.get()));
		}
		return builder.toString();
	}

}

