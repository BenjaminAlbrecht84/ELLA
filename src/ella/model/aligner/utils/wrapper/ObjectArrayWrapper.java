package ella.model.aligner.utils.wrapper;

import java.util.Arrays;

public class ObjectArrayWrapper {
	
	private final Object[] data;

	public ObjectArrayWrapper(Object[] data) {
		if (data == null) {
			throw new NullPointerException();
		}
		this.data = data;
	}

	@Override
	public boolean equals(Object other) {
		if (!(other instanceof ObjectArrayWrapper)) {
			return false;
		}
		return Arrays.equals(data, ((ObjectArrayWrapper) other).data);
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(data);
	}

	public Object[] getData() {
		return data;
	}
	
	public int getLength() {
		return data.length;
	}

}

